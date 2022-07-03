
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from koopmans.bands import Band
from ._workflow import Workflow
from koopmans import utils
from koopmans.settings import KoopmansCPSettingsDict

from koopmans import ml_utils
import numpy as np

from typing import List, Dict, Any, Tuple
from pathlib import Path

import os

from sklearn.metrics import mean_absolute_error as mae

class MLFiitingWorkflow(Workflow):

    def __init__(self, calc_that_produced_orbital_densities, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.method_to_extract_from_binary        = 'from_ki'
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities
        
        
        ml_dir                                    = self.calc_that_produced_orbital_densities.directory / 'ml' / 'TMP'
        dir_suffix                                = '_'.join(str(x) for x in [self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max])
        self.dirs                                 = {
            'ml'        : ml_dir,
            'xml'       : ml_dir / 'xml', 
            'alphas'    : ml_dir / 'alphas', 
            'betas'     : ml_dir / 'betas', 
            'coeff'     : ml_dir / ('coefficients_'  + dir_suffix),
            'coeff_orb' : ml_dir / ('coefficients_'  + dir_suffix) / 'coeff_orb',
            'coeff_tot' : ml_dir / ('coefficients_'  + dir_suffix) / 'coeff_tot',
            'power'     : ml_dir / ('power_spectra_' + dir_suffix)
        }


        for dir in self.dirs.values():
            utils.system_call(f'mkdir -p {dir}')



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
        self.extract_input_vector_for_ml_model()


    def extract_input_vector_for_ml_model(self):
        self.convert_bin2xml()
        self.compute_decomposition()
        self.compute_power_spectrum()
            

    # TODO: currently it extracts all orbitals in [1,..,self.n_bands_to_extract[0]] instead of the indices given by self.bands_to_extract
    def convert_bin2xml(self):
        if self.method_to_extract_from_binary == 'from_ki':
            orbital_densities_bin_dir            = self.calc_that_produced_orbital_densities.parameters.outdir/ f'kc_{self.calc_that_produced_orbital_densities.parameters.ndw}.save'
        else:
            raise NotImplementedError(f'Currently it is only implemented to extract the real space orbital densities from the ki-trial calculation after the initial wannier-calculation')

        calculation_title = 'conversion binary->xml of real-space-densities'
        is_complete = self.check_if_bin2xml_is_complete()
        
        if not self.parameters.from_scratch and is_complete:
            self.print(f'Not running {calculation_title} as it is already complete')
        else:
            self.print(f'Running {calculation_title}...', end='', flush=True)
            command  = str(calculators.bin_directory / 'bin2xml_real_space_density.x ') + ' '.join(str(x) for x in [orbital_densities_bin_dir, self.dirs['xml'], self.n_bands_to_extract[0], self.n_bands_to_extract[1], self.nspin_to_extract])
            utils.system_call(command)
            self.print(f' done')
            

    def check_if_bin2xml_is_complete(self) -> bool:
        if (len([name for name in os.listdir(self.dirs['xml']) if os.path.isfile(self.dirs['xml'] / name)])==(self.n_bands_to_extract[0] + self.n_bands_to_extract[1]+1)):
            return True
        else:
            return False

    def compute_decomposition(self):
        calculation_title = 'computation of decomposition of real-space-density'
        is_complete       = self.check_if_compute_decomposition_is_complete()
        if not self.parameters.from_scratch and is_complete:
            self.print(f'Not running {calculation_title} as it is already complete')
        else:
            self.print(f'Running {calculation_title}...', end='', flush=True)

            self.r_cut = min(self.atoms.get_cell_lengths_and_angles()[:3])#/2.5 #the maximum radius will be set to the minimum of self.r_cut and half of the cell-size later on
            if self.method_to_extract_from_binary == 'from_ki':
                # Store the original w90 calculations
                w90_calcs = [c for c in self.calculations if isinstance(c, calculators.Wannier90Calculator) and c.command.flags == ''][-len(self.projections):]
                
                # TODO: implement also the spin-unpolarized case?
                calc_presets_occ = 'occ'
                calc_presets_emp = 'emp'
                centers_occ      = np.array([center for c in w90_calcs for center in c.results['centers'] if calc_presets_occ in c.directory.name])
                centers_emp      = np.array([center for c in w90_calcs for center in c.results['centers'] if calc_presets_emp in c.directory.name]) 
                centers = np.concatenate([centers_occ, centers_emp])

            else: 
                raise NotImplementedError(f'Currently it is only implemented to extract the real space orbital densities from the ki-trial calculation after the initial wannier-calculation')
            
            ml_utils.precompute_radial_basis(self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, self.dirs)
            ml_utils.func_compute_decomposition(self.parameters.n_max, self.parameters.l_max, self.parameters.r_min, self.parameters.r_max, self.r_cut, self.dirs, self.bands_to_extract, self.atoms, centers)
            self.print(f' done')
    
    def check_if_compute_decomposition_is_complete(self) -> bool:
        if (len([name for name in os.listdir(self.dirs['coeff'] / 'coeff_tot') if os.path.isfile(self.dirs['coeff'] / 'coeff_tot' / name)])==(self.n_bands_to_extract[0] + self.n_bands_to_extract[1])):
            return True
        else:
            return False
    
    def compute_power_spectrum(self):
        calculation_title = 'computation of power spectrum'
        self.print(f'Running {calculation_title}...', end='', flush=True)
        ml_utils.main_compute_power(self.parameters.n_max, self.parameters.l_max, self.dirs, self.bands_to_extract)
        self.print(f' done')

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
        return np.loadtxt(self.dirs['power'] / f"power_spectrum.orbital.{filled_str}.{band.index}.txt")


    
    
    def use_prediction(self) -> bool:
        # defualt is to not use the prediction
        use_prediction = False
        if self.parameters.criterium == 'after_fixed_num_of_snapshots':
            if self.parameters.current_snapshot < self.parameters.number_of_snapshots:
                use_prediction = False
            else:
                use_prediction = True
        else: 
             raise NotImplementedError(f'criterium = {self.parameters.criterium} is currently not implemented')
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
        utils.write_alpha_file(self.dirs['ml'], alphas, fillings)

    
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
        assert isinstance(band.index, int)
        return alphas[0][band.index-1], 0.0

