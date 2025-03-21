"""
Written by Yannick Schubert Jul 2022
"""
import os
from pathlib import Path
from typing import Any, Dict, Generator, List, Tuple
from pydantic import ConfigDict

import numpy as np
from sklearn.metrics import mean_absolute_error as mae

from koopmans import calculators, ml, utils
from koopmans.files import File
from koopmans.process_io import IOModel
from koopmans.processes.power_spectrum import (
    ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)
from koopmans.status import Status

from ._workflow import Workflow


class SelfHartreeOutput(IOModel):
    descriptors: List[float]


class SelfHartreeWorkflow(Workflow[SelfHartreeOutput]):

    output_model = SelfHartreeOutput

    def _run(self) -> None:
        assert self.bands
        self.outputs = self.output_model(descriptors=[b.self_hartree for b in self.bands.to_solve])
        self.status = Status.COMPLETED
        return


class PowerSpectrumDecompositionOutput(IOModel):
    descriptors: List[File]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class PowerSpectrumDecompositionWorkflow(Workflow[PowerSpectrumDecompositionOutput]):

    output_model = PowerSpectrumDecompositionOutput  # type: ignore

    def __init__(self, calc_that_produced_orbital_densities: calculators.Calc, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Specify from which calculation we extract the real space densities. Currently only
        # the KI-calculation at the beginning of the alpha calculations is a valid option
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities

    def _run(self) -> None:
        """
        Runs the PowerSpectrumDecomposition workflow.

        If the input data are orbital densities, this workflow consists of three steps:
        1) converting the binary files containing real space densities to xml-files
        2) reading the real space densities from the xml files and calculating the decomposition into spherical
           harmonics and radial basis functions
        3) computing the power spectra of the resulting coefficient vectors
        """

        if self.ml.estimator == 'mean':
            raise ValueError('A mean estimator does not require input data; this `PowerSpectrumDecompositionWorkflow` '
                             'should not have been called')

        # Specify for which bands we want to compute the decomposition
        assert self.bands
        self.num_bands_occ = [len([band for band in self.bands if (band.filled and band.spin == spin)])
                              for spin in [0, 1]]
        self.num_bands_to_extract = [len([band for band in self.bands.to_solve if band.filled == filled])
                                     for filled in [True, False]]

        if self.bands.spin_polarized:
            self.nspin_to_extract = 2
        else:
            self.nspin_to_extract = 1

        assert self.ml.descriptor == 'orbital_density'
        self.extract_input_vector_from_orbital_densities()

        self.status = Status.COMPLETED

        return

    def extract_input_vector_from_orbital_densities(self):
        """
        Performs the three steps of the PowerSpectrumDecomposition workflow
        """

        # Convert the binary files to xml format
        bin2xml_workflow = ConvertOrbitalFilesToXMLWorkflow.fromparent(
            self, calc_that_produced_orbital_densities=self.calc_that_produced_orbital_densities)
        bin2xml_workflow.run()

        # Extract the coefficients from the xml files
        if self.parameters.spin_polarized:
            wannier_centers = [b.center for b in self.bands]
        else:
            wannier_centers = [b.center for b in self.bands if b.spin == 0]

        decomposition_process = ExtractCoefficientsFromXMLProcess(n_max=self.ml.n_max,
                                                                  l_max=self.ml.l_max,
                                                                  r_min=self.ml.r_min,
                                                                  r_max=self.ml.r_max,
                                                                  r_cut=min(
                                                                      self.atoms.get_cell_lengths_and_angles()[:3]),
                                                                  wannier_centers=wannier_centers,
                                                                  bands=self.bands.to_solve,
                                                                  cell=self.atoms.cell,
                                                                  total_density_xml=bin2xml_workflow.outputs.total_density,
                                                                  orbital_densities_xml=bin2xml_workflow.outputs.orbital_densities)
        status = self.run_steps(decomposition_process)
        if status != Status.COMPLETED:
            return

        orb_coeffs = decomposition_process.outputs.orbital_coefficients
        tot_coeffs = decomposition_process.outputs.total_coefficients
        descriptors = []
        for i, (orb_coeff, tot_coeff) in enumerate(zip(orb_coeffs, tot_coeffs)):
            power_spectrum_process = ComputePowerSpectrumProcess(
                n_max=self.ml.n_max,
                l_max=self.ml.l_max,
                orbital_coefficients=orb_coeff,
                total_coefficients=tot_coeff,
            )
            power_spectrum_process.name += '_orbital_' + str(i + 1)
            status = self.run_steps(power_spectrum_process)
            if status != Status.COMPLETED:
                return
            descriptors.append(power_spectrum_process.outputs.power_spectrum)

        self.outputs = self.output_model(descriptors=descriptors)

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        calc_that_produced_orbital_densities = dct.pop('calc_that_produced_orbital_densities')
        return super(PowerSpectrumDecompositionWorkflow, cls).fromdict(
            dct, calc_that_produced_orbital_densities=calc_that_produced_orbital_densities, **kwargs)


class ConvertOrbitalFilesToXMLOutput(IOModel):
    total_density: File
    orbital_densities: List[File]
    model_config = ConfigDict(arbitrary_types_allowed=True)


class ConvertOrbitalFilesToXMLWorkflow(Workflow[ConvertOrbitalFilesToXMLOutput]):

    output_model = ConvertOrbitalFilesToXMLOutput  # type: ignore

    def __init__(self, calc_that_produced_orbital_densities, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities

    def _run(self) -> None:
        """
        Converts the binary files produced by a previous calculation to python-readable xml files.
        """

        from koopmans.processes.bin2xml import Bin2XMLProcess

        # Convert total density to XML
        binary = self.calc_that_produced_orbital_densities.write_directory / 'charge-density.dat'
        bin2xml_total_density = Bin2XMLProcess(name='bin2xml_total_density', binary=binary)
        status = self.run_steps(bin2xml_total_density)
        if status != Status.COMPLETED:
            return

        # Convert orbital densities to XML
        orbital_densities: List[File] = []
        assert self.bands
        for band in self.bands.to_solve:
            if band.filled:
                occ_id = 'occ'
            else:
                occ_id = 'emp'
            binary = self.calc_that_produced_orbital_densities.write_directory / f'real_space_orb_density.{occ_id}.{band.spin}.{band.index:05d}.dat'

            bin2xml_orbital_density = Bin2XMLProcess(
                name=f'bin2xml_{occ_id}_spin_{band.spin}_orb_{band.index}_density', binary=binary)
            status = self.run_steps(bin2xml_orbital_density)
            if status != Status.COMPLETED:
                return
            orbital_densities.append(bin2xml_orbital_density.outputs.xml)

        self.outputs = self.output_model(total_density=bin2xml_total_density.outputs.xml,
                                         orbital_densities=orbital_densities)

        self.status = Status.COMPLETED
