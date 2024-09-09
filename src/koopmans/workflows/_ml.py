"""
Written by Yannick Schubert Jul 2022
"""
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from sklearn.metrics import mean_absolute_error as mae

from koopmans import calculators, ml, utils
from koopmans.bands import Band
from koopmans.files import FilePointer
from koopmans.outputs import OutputModel
from koopmans.processes.power_spectrum import (
    ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)
from koopmans.settings import KoopmansCPSettingsDict

from ._workflow import Workflow


class SelfHartreeOutput(OutputModel):
    descriptors: List[float]


class SelfHartreeWorkflow(Workflow):

    output_model = SelfHartreeOutput  # type: ignore

    def _run(self):
        self.outputs = self.output_model(descriptors=[b.self_hartree for b in self.bands.to_solve])


class PowerSpectrumDecompositionOutput(OutputModel):
    descriptors: List[FilePointer]

    class Config:
        arbitrary_types_allowed = True


class PowerSpectrumDecompositionWorkflow(Workflow):

    output_model = PowerSpectrumDecompositionOutput  # type: ignore

    def __init__(self, calc_that_produced_orbital_densities: calculators.Calc, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Specify from which calculation we extract the real space densities. Currently only
        # the KI-calculation at the beginning of the alpha calculations is a valid option
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities

    def _run(self):
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
        self.run_process(decomposition_process)

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
            self.run_process(power_spectrum_process)
            descriptors.append(power_spectrum_process.outputs.power_spectrum)

        self.outputs = self.output_model(descriptors=descriptors)

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        calc_that_produced_orbital_densities = dct.pop('calc_that_produced_orbital_densities')
        return super(PowerSpectrumDecompositionWorkflow, cls).fromdict(
            dct, calc_that_produced_orbital_densities=calc_that_produced_orbital_densities, **kwargs)


class ConvertOrbitalFilesToXMLOutput(OutputModel):
    total_density: FilePointer
    orbital_densities: List[FilePointer]

    class Config:
        arbitrary_types_allowed = True


class ConvertOrbitalFilesToXMLWorkflow(Workflow):

    output_model = ConvertOrbitalFilesToXMLOutput  # type: ignore

    def __init__(self, calc_that_produced_orbital_densities, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities

    def _run(self):
        """
        Converts the binary files produced by a previous calculation to python-readable xml files.
        """

        from koopmans.processes.bin2xml import Bin2XMLProcess

        # Convert total density to XML
        binary = FilePointer(self.calc_that_produced_orbital_densities,
                             self.calc_that_produced_orbital_densities.write_directory / 'charge-density.dat')
        bin2xml_total_density = Bin2XMLProcess(name='bin2xml_total_density', binary=binary)
        self.run_process(bin2xml_total_density)

        # Convert orbital densities to XML
        orbital_densities: List[FilePointer] = []
        for band in self.bands.to_solve:
            if band.filled:
                occ_id = 'occ'
            else:
                occ_id = 'emp'
            binary = FilePointer(self.calc_that_produced_orbital_densities,
                                 self.calc_that_produced_orbital_densities.write_directory / f'real_space_orb_density.{occ_id}.{band.spin}.{band.index:05d}.dat')

            bin2xml_orbital_density = Bin2XMLProcess(
                name=f'bin2xml_{occ_id}_spin_{band.spin}_orb_{band.index}_density', binary=binary)
            self.run_process(bin2xml_orbital_density)
            orbital_densities.append(bin2xml_orbital_density.outputs.xml)

        self.outputs = self.output_model(total_density=bin2xml_total_density.outputs.xml,
                                         orbital_densities=orbital_densities)
