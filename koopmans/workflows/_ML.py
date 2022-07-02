
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from koopmans.bands import Band
from ._workflow import Workflow
from koopmans import utils
from koopmans.settings import KoopmansCPSettingsDict

from koopmans import ML_utils
import numpy as np

from typing import List, Dict, Any, Tuple
from pathlib import Path

from sklearn.metrics import mean_absolute_error as mae

class MLFiitingWorkflow(Workflow):

    def __init__(self, calc_that_produced_orbital_densities, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.method_to_extract_from_binary        = 'from_ki'
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities
        self.ML_dir                               = self.calc_that_produced_orbital_densities.directory / 'ML' / 'TMP'
        self.predicted_alphas                     = []
        self.calculated_alphas                    = []
        self.fillings_of_predicted_alphas         = []
        self.use_predictions                      = []
        


    def _run(self):
        
        self.bands_to_extract     = self.bands.to_solve
        self.n_bands_to_extract   = [len([band for band in self.bands_to_extract if band.filled==filled]) for filled in [True, False]]
        if self.bands.spin_polarised:
            self.nspin_to_extract = 2
        else:
            self.nspin_to_extract = 1
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.print(f'Computing the power spectra from the real space densities...', end='', flush=True)
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()
        self.print(f' done')
            

    # TODO: currently it extracts all orbitals in [1,..,self.n_bands_to_extract[0]] instead of the indices given by self.bands_to_extract
    def convert_binary_to_xml(self):
        orbital_densities_bin_dir            = self.calc_that_produced_orbital_densities.parameters.outdir/ f'kc_{self.calc_that_produced_orbital_densities.parameters.ndw}.save'
        

        if self.method_to_extract_from_binary == 'from_ki':
            utils.system_call(f'mkdir -p {self.ML_dir}')
            command  = str(calculators.bin_directory / 'bin2xml_real_space_density.x ') + ' '.join(str(x) for x in [orbital_densities_bin_dir, self.ML_dir, self.n_bands_to_extract[0], self.n_bands_to_extract[1], self.nspin_to_extract])
            utils.system_call(command)

    def compute_decomposition(self):
        self.r_cut = min(self.atoms.get_cell_lengths_and_angles()[:3])#/2.5 #the maximum radius will be set to the minimum of self.r_cut and half of the cell-size later on
        print("r_cut = ", self.r_cut)
        if self.method_to_extract_from_binary == 'from_ki':
            centers_occ = np.array(self.calculations[-11].results['centers'])
            centers_emp = np.array(self.calculations[-8].results['centers'])
            centers     = np.concatenate([centers_occ, centers_emp])
        else: 
             raise ValueError(f'Currently it is only implemented to extract the real space orbital densities from the ki-trial calculation after the initial wannier-calculation')
        
        ML_utils.precompute_radial_basis(self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, self.ML_dir)
        ML_utils.func_compute_decomposition(self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, self.r_cut, self.ML_dir, self.bands_to_extract, self.atoms, centers)
    
    def compute_power_spectrum(self):
        self.dir_power = self.ML_dir / f'power_spectra_{self.parameters.n_max}_{self.parameters.l_max}_{self.parameters.r_min}_{self.parameters.r_max}'
        ML_utils.main_compute_power(self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, self.ML_dir, self.dir_power, self.bands_to_extract)

    def predict(self, band: Band) -> float:
        self.print('Predicting screening parameter')
        # TODO: currently only capable of making one prediction at a time
        power_spectrum = self.load_power_spectrum(band)
        y_predict      = self.ml_model.predict(power_spectrum)[0]
        self.predicted_alphas.append(y_predict)
        self.fillings_of_predicted_alphas.append(band.filled)

        return y_predict
    
    def train(self):
        self.print('Training the ML model')
        self.ml_model.train()
    
    def add_training_data(self, band: Band):
        self.print('Adding this orbital to the training data')
        power_spectrum = self.load_power_spectrum(band)
        alpha          = band.alpha
        assert isinstance(alpha, float)
        self.calculated_alphas.append(alpha)
        self.ml_model.add_training_data(power_spectrum, alpha)

    def load_power_spectrum(self, band: Band) -> np.ndarray:
        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'
        return np.loadtxt(self.dir_power / f"power_spectrum.orbital.{filled_str}.{band.index}.txt")


    
    
    def use_prediction(self) -> bool:
        # defualt is to not use the prediction
        use_prediction = False
        if self.parameters.criterium == 'after_fixed_num_of_snapshots':
            if self.parameters.current_snapshot < self.parameters.number_of_snapshots:
                use_prediction = False
            else:
                use_prediction = True
        else: 
             raise ValueError(f'criterium = {self.parameters.criterium} is currently not implemented')
        if use_prediction:
            self.print('The prediction-criterium is satisfied -> Use the predicted screening parameter')
        else:
            self.print('The prediction-criterium is not yet satsified -> Compute the screening parameter ab initio')
        
        self.use_predictions.append(use_prediction)
        
        return use_prediction


    def write_predicted_alphas(self):
        if self.nspin_to_extract==1:
            duplicate = 2
        else:
            duplicate = 1
        alphas   = [a for _ in range(duplicate) for a in self.predicted_alphas]
        fillings = [f for _ in range(duplicate) for f in self.fillings_of_predicted_alphas]
        utils.write_alpha_file(self.ML_dir, alphas, fillings)

    
    def print_error_of_single_orbital(self, alpha_predicted:float, alpha_calculated:float, indent: int = 0):
        # Printing out the error of the predicted alpha
        utils.indented_print('\npredicted  screening parameter: {0:.5f}'.format(alpha_predicted), indent=indent)
        utils.indented_print(  'calculated screening parameter: {0:.5f}'.format(alpha_calculated), indent=indent)
        utils.indented_print(  'absoulute error               : {0:.5f}'.format(np.abs(alpha_predicted-alpha_calculated)), indent=indent)
        utils.indented_print('')

    def print_error_of_all_orbitals(self, indent: int = 0):
        # Printing out a progress summary
        # utils.indented_print(f'\nerror of predictions', indent=indent)
        # utils.indented_print('calculated ' + str(self.calculated_alphas), indent=indent)
        # utils.indented_print('predicted  ' + str(self.predicted_alphas), indent=indent)
        # utils.indented_print(str(self.bands.alpha_history), indent=indent)
        # utils.indented_print('')
        print("Yannick Debug: len(predicted alphas) = ", len(self.predicted_alphas))
        print("Yannick Debug: len(calculated alphas) = ", len(self.calculated_alphas))
        utils.indented_print('\nThe mean absolut error of the predicted screening parameters of this snapshot is {0:.5f}'.format(mae(self.predicted_alphas, self.calculated_alphas)), indent=indent)
        utils.indented_print('')
    
    def get_alpha_from_file_for_debugging(self, band:Band) -> Tuple[float, float]:
        flat_alphas = utils.read_alpha_file(Path())
        params = self.master_calc_params['kcp']
        assert isinstance(params, KoopmansCPSettingsDict)
        alphas = calculators.convert_flat_alphas_for_kcp(flat_alphas, params)
        assert isinstance(band.index, float)
        return alphas[0][band.index-1], 0.0

