"""
Written by Yannick Schubert Jul 2022
"""
import copy
import os
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
from deepdiff import DeepDiff
from sklearn.metrics import mean_absolute_error as mae

from koopmans import calculators, ml, utils
from koopmans.bands import Band
from koopmans.settings import KoopmansCPSettingsDict

from ._workflow import Workflow


class MLFittingWorkflow(Workflow):

    def __init__(self, calc_that_produced_orbital_densities=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Specify from which calculation we extract the real space densities. Currently only
        # the ki-calculation at the beginning of the alpha calculations is a valid option
        if calc_that_produced_orbital_densities is None:
            raise ValueError(
                "please provide this workflow with a calculation that produced the real space orbital densities.")
        else:
            self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities

        # Define and create the different subdirectories for the outputs of the ML-workflow
        ml_dir = self.calc_that_produced_orbital_densities.directory / 'ml' / 'TMP'
        dir_suffix = '_'.join(str(x) for x in [self.ml.n_max,
                              self.ml.l_max, self.ml.r_min, self.ml.r_max])
        if self.ml.input_data_for_ml_model == 'orbital_density':
            self.dirs = {
                'ml': ml_dir,
                'xml': ml_dir / 'xml',
                'alphas': ml_dir / 'alphas',
                'betas': ml_dir / 'betas',
                'coeff': ml_dir / ('coefficients_' + dir_suffix),
                'coeff_orb': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_orb',
                'coeff_tot': ml_dir / ('coefficients_' + dir_suffix) / 'coeff_tot',
                'power': ml_dir / ('power_spectra_' + dir_suffix)
            }
        elif self.ml.input_data_for_ml_model == 'self_hartree':
            self.dirs = {
                'ml': ml_dir,
                'SH': ml_dir / 'SH'
            }
        else:
            raise ValueError(
                f"{self.ml.input_data_for_ml_model} is currently not implemented as a valid input for the ml model.")

        for dir in self.dirs.values():
            dir.mkdir(parents=True, exist_ok=True)

        # normally this dictonary contains the power spectrum
        self.input_vectors_for_ml: Dict[str, np.ndarray] = {}

        # initialize lists for the results of the prediction and the calculations
        self.predicted_alphas = []
        self.calculated_alphas = []
        self.fillings_of_predicted_alphas = []
        self.use_predictions = []

    def _run(self):
        """
        Runs the MLFitting workflow.

        If the input data are orbital densities, this workflow consists of three steps:
        1) converting the binary files containing real space densities to xml-files
        2) reading the real space densities from the xml files and calculating the decomposition into spherical
           harmonics and radial basis functions
        3) computing the power spectra of the resulting coefficient vectors
        """

        # Specify for which bands we want to compute the decomposition
        self.bands_to_extract = self.bands.to_solve
        self.num_bands_occ = [len([band for band in self.bands if (band.filled == True and band.spin == spin)])
                              for spin in [0, 1]]
        self.num_bands_to_extract = [len([band for band in self.bands_to_extract if band.filled == filled])
                                     for filled in [True, False]]

        if self.bands.spin_polarized:
            self.nspin_to_extract = 2
        else:
            self.nspin_to_extract = 1

        if self.ml.type_of_ml_model == 'mean':
            return  # this model needs no X-data
        else:
            if self.ml.input_data_for_ml_model == 'orbital_density':
                # start the actual three steps
                self.extract_input_vector_from_orbital_densities()
            elif self.ml.input_data_for_ml_model == 'self_hartree':
                # get the self-Hartree energies
                self.extract_input_vector_from_self_hartrees()
            else:
                raise ValueError(
                    f"{self.ml.input_data_for_ml_model} is currently not implemented as a valid input for the ml model.")

    def extract_input_vector_from_self_hartrees(self):
        """
        Extracts the self-Hartree energies from the corresponding calculation.
        """
        SH = self.calc_that_produced_orbital_densities.results['orbital_data']['self-Hartree']
        for band in self.bands_to_extract:
            if band.filled:
                filled_str = 'occ'
            else:
                filled_str = 'emp'
            np.savetxt(self.dirs['SH'] / f"SH.orbital.{filled_str}.{band.index}.txt",
                       np.array([SH[band.spin][band.index-1]]))
            self.input_vectors_for_ml[f"SH.orbital.{filled_str}.{band.index}"] = SH

        return

    def extract_input_vector_from_orbital_densities(self):
        """
        Performs the three steps of the MLFitting workflow
        """

        self.convert_bin2xml()
        self.compute_decomposition()
        self.compute_power_spectrum()

    def convert_bin2xml(self):
        """
        Converts the binary files produced by a previous calculation to python-readable xml files.
        """
        orbital_densities_bin_dir = self.calc_that_produced_orbital_densities.parameters.outdir / \
            f'kc_{self.calc_that_produced_orbital_densities.parameters.ndw}.save'

        calculation_title = 'Converting the real-space densities from binary to xml format'
        is_complete = self.check_if_bin2xml_is_complete()

        if not self.parameters.from_scratch and is_complete:
            self.print(f'Not {calculation_title.lower()} as it is already complete')
        else:
            self.print(f'{calculation_title}...', end='', flush=True)

            # Convert total density to XML
            command = f'bin2xml.x {orbital_densities_bin_dir}/charge-density.dat {self.dirs["xml"]}/charge-density.xml'
            utils.system_call(command)

            # Convert orbital densities to XML
            for band in self.bands_to_extract:
                if band.filled:
                    occ_id = 'occ'
                else:
                    occ_id = 'emp'
                dat_seed = f'real_space_orb_density.{occ_id}.{band.spin}.{band.index:05d}'
                xml_seed = f'orbital.{occ_id}.{band.spin}.{band.index:05d}'

                command = f'bin2xml.x {orbital_densities_bin_dir}/{dat_seed}.dat {self.dirs["xml"]}/{xml_seed}.xml'
                utils.system_call(command)

            self.print(f' done')

    def check_if_bin2xml_is_complete(self) -> bool:
        """
        Checks if the bin2xml conversion was already performed
        """

        # If there are as many xml files as there are bands to solve, the calculation was already completed
        return len([name for name in os.listdir(self.dirs['xml']) if os.path.isfile(self.dirs['xml'] / name)]) \
            == (self.num_bands_to_extract[0] + self.num_bands_to_extract[1] + 1)

    def compute_decomposition(self):
        """
        Performs the decomposition into radial basis functions and spherical harmonics
        """

        calculation_title = 'the decomposition of the real-space density'
        is_complete = self.check_if_compute_decomposition_is_complete()
        if not self.parameters.from_scratch and is_complete:
            self.print(f'Not calculating {calculation_title} as it is already complete')
        else:
            self.print(f'Calculating {calculation_title}...', end='', flush=True)

            # the maximum radius will be set to the minimum of self.r_cut and half of the cell-size later on
            self.r_cut = min(self.atoms.get_cell_lengths_and_angles()[:3])

            # Extract the wannier-centres
            w90_calcs = [c for c in self.calculations if isinstance(
                c, calculators.Wannier90Calculator) and c.command.flags == ''][-len(self.projections):]

            if self.parameters.spin_polarized:
                spins = ['up', 'down']
            else:
                spins = [None]

            centers_list = []
            for spin in spins:
                for filling in ['occ', 'emp']:
                    calc_presets = filling
                    if spin:
                        calc_presets += '_' + spin
                    centers_list.append(np.array([center for c in w90_calcs for center in c.results['centers']
                                                  if calc_presets in c.directory.name]))
            centers = np.concatenate(centers_list)

            ml.precompute_parameters_of_radial_basis(self.ml.n_max, self.ml.l_max,
                                                     self.ml.r_min, self.ml.r_max, self.dirs)
            ml.compute_decomposition(self.ml.n_max, self.ml.l_max, self.ml.r_min,
                                     self.ml.r_max, self.r_cut, self.dirs, self.bands_to_extract, self.atoms, centers)
            self.print(f' done')

    def check_if_compute_decomposition_is_complete(self) -> bool:
        """
        Checks if the expansion coefficients were already calculated
        """

        # If there are as many coefficient files as there are bands to solve, the calculation was already completed
        return len([name for name in os.listdir(self.dirs['coeff'] / 'coeff_tot') if
                    os.path.isfile(self.dirs['coeff'] / 'coeff_tot' / name)]) == sum(self.num_bands_to_extract)

    def compute_power_spectrum(self):
        """
        Performs the computation of the power spectra.
        """

        self.print('Calculating the power spectrum...', end='', flush=True)
        ml.compute_power(self.ml.n_max, self.ml.l_max, self.dirs,
                         self.bands_to_extract, self.input_vectors_for_ml)
        self.print(' done')

    def get_input_data(self, band) -> np.ndarray:
        """
        Loads the input data depending on the ML model
        """
        if self.ml.type_of_ml_model == 'mean':
            input_data = np.array([1.0])  # dummy value
        else:
            if self.ml.input_data_for_ml_model == 'orbital_density':
                input_data = self.load_power_spectrum(band)
            elif self.ml.input_data_for_ml_model == 'self_hartree':
                input_data = self.load_SH(band)
            else:
                raise ValueError(f"{self.ml.input_data_for_ml_model} is currently not implemented as a valid input "
                                 "for the ML model")

        return input_data

    def predict(self, band: Band) -> float:
        """
        Make the prediction for one band.
        """

        input_data = self.get_input_data(band)

        if self.ml.occ_and_emp_together:
            y_predict = self.ml.ml_model.predict(input_data)[0]
        else:
            if band.filled:
                y_predict = self.ml.ml_model_occ.predict(input_data)[0]
            else:
                y_predict = self.ml.ml_model_emp.predict(input_data)[0]

        self.predicted_alphas.append(y_predict)
        self.fillings_of_predicted_alphas.append(band.filled)

        return y_predict

    def train(self):
        """
        Reset the model and train the ML-model (including the StandardScaler) with all training data added so far
        """

        self.print('Training the ML model')
        if self.ml.occ_and_emp_together:
            self.ml.ml_model.train()
        else:
            self.ml.ml_model_occ.train()
            self.ml.ml_model_emp.train()

    def add_training_data(self, band: Band):
        """
        Add training data to the ML-model.
        """

        self.print('Adding this orbital to the training data')
        input_data = self.get_input_data(band)
        alpha = band.alpha
        assert isinstance(alpha, float)
        self.calculated_alphas.append(alpha)
        if self.ml.occ_and_emp_together:
            self.ml.ml_model.add_training_data(input_data, alpha)
        else:
            if band.filled:
                self.ml.ml_model_occ.add_training_data(input_data, alpha)
            else:
                self.ml.ml_model_emp.add_training_data(input_data, alpha)

    def load_power_spectrum(self, band: Band) -> np.ndarray:
        """
        Reads the power spectrum of one band from a file
        """

        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'
        return self.input_vectors_for_ml[f"power_spectrum.orbital.{filled_str}.{band.index}"]

    def load_SH(self, band: Band) -> np.ndarray:
        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'
        return self.input_vectors_for_ml[f"SH.orbital.{filled_str}.{band.index}"]
        # np.loadtxt(self.dirs['SH'] / f"SH.orbital.{filled_str}.{band.index}.txt")

    def use_prediction(self) -> bool:
        """
        Check if the prediction criterium specified by the user is satisfied.

        If True: use the ML-prediction for the alpha value
        If False: compute this alpha value ab-initio
        """

        # Default is to not use the prediction
        use_prediction = False
        if self.ml.criterium == 'after_fixed_num_of_snapshots':
            if self.ml.current_snapshot < self.ml.number_of_training_snapshots:
                use_prediction = False
            else:
                use_prediction = True
        else:
            raise NotImplementedError(f'criterium = {self.ml.criterium} is currently not implemented')
        if use_prediction:
            self.print('Predicting the screening parameter with the ML model')

        # Store whether the prediction was used or not
        self.use_predictions.append(use_prediction)

        return use_prediction

    def write_predicted_alphas(self):
        """
        Wrapper to write the predicted alpha values to a file
        """

        if self.nspin_to_extract == 1:
            duplicate = 2
        else:
            duplicate = 1
        alphas = [a for _ in range(duplicate) for a in self.predicted_alphas]
        fillings = [f for _ in range(duplicate) for f in self.fillings_of_predicted_alphas]
        utils.write_alpha_file(self.dirs['ml'], alphas, fillings)

    def print_error_of_single_orbital(self, alpha_predicted: float, alpha_calculated: float, indent: int = 0):
        """
        Prints a summary of the prediction of the alphas values of one band
        """

        utils.indented_print(f'\npredicted screening parameter:  {alpha_predicted:.5f}', indent=indent)
        utils.indented_print(f'calculated screening parameter: {alpha_calculated:.5f}', indent=indent)
        utils.indented_print(f'absolute error:                 {np.abs(alpha_predicted-alpha_calculated):.5f}\n',
                             indent=indent)

    def print_error_of_all_orbitals(self, indent: int = 0):
        """
        Prints a summary of the prediction of the alphas values of all bands
        """
        utils.indented_print('\nThe mean absolute error of the predicted screening parameters of this snapshot is '
                             f'{mae(self.predicted_alphas, self.calculated_alphas):.4f}\n', indent=indent)

    def get_alpha_from_file_for_debugging(self, band: Band) -> Tuple[float, float]:
        """
        Auxillary function to test the ML workflow with screening parameters from a file instead of computing them
        ab initio
        """

        flat_alphas = utils.read_alpha_file(Path())
        params = self.calculator_parameters['kcp']
        assert isinstance(params, KoopmansCPSettingsDict)
        alphas = calculators.convert_flat_alphas_for_kcp(flat_alphas, params)
        if self.parameters.spin_polarized:
            raise NotImplementedError('Need to check implementation')
        assert isinstance(band.index, int)
        return alphas[0][band.index-1], 0.0

    @classmethod
    def fromdict(cls, dct: Dict[str, Any], **kwargs) -> Workflow:
        calc_that_produced_orbital_densities = dct.pop('calc_that_produced_orbital_densities')
        return super(MLFittingWorkflow, cls).fromdict(
            dct, calc_that_produced_orbital_densities=calc_that_produced_orbital_densities, **kwargs)

    def __eq__(self, other):
        return DeepDiff(self, other) == {}
