
from abc import ABC, abstractmethod
from xmlrpc.client import boolean

from koopmans import calculators
from ._workflow import Workflow
from koopmans import utils

from koopmans import ML_utils
import numpy as np

from typing import List, Dict, Any

import copy


class MLFiitingWorkflow(Workflow):

    def __init__(self, calc_that_produced_orbital_densities, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.method_to_extract_from_binary        = 'from_ki'
        self.calc_that_produced_orbital_densities = calc_that_produced_orbital_densities
        self.ML_dir                               = self.calc_that_produced_orbital_densities.directory / 'ML' / 'TMP'
        ML_params = self.master_calc_params['ML']
        self.n_max = ML_params.n_max
        self.l_max = ML_params.l_max
        self.r_min = ML_params.r_min
        self.r_max = ML_params.r_max


    def _run(self):
        
        self.bands_to_extract     = self.bands.to_solve
        self.n_bands_to_extract   = [len([band for band in self.bands_to_extract if band.filled==filled]) for filled in [True, False]]
        if self.bands.spin_polarised:
            self.nspin_to_extract = 2
        else:
            self.nspin_to_extract = 1
        self.extract_input_vector_for_ML_model()


    def extract_input_vector_for_ML_model(self):
        self.convert_binary_to_xml()
        self.compute_decomposition()
        self.compute_power_spectrum()
            

    # TODO: currently it extracts all orbitals in [1,..,self.n_bands_to_extract[0]] instead of the indices given by self.bands_to_extract
    def convert_binary_to_xml(self):
        print("Convert binary to xml")
        orbital_densities_bin_dir            = self.calc_that_produced_orbital_densities.parameters.outdir/ f'kc_{self.calc_that_produced_orbital_densities.parameters.ndw}.save'
        

        if self.method_to_extract_from_binary == 'from_ki':
            utils.system_call(f'mkdir -p {self.ML_dir}')
            command  = str(calculators.bin_directory / 'bin2xml_real_space_density.x ') + ' '.join(str(x) for x in [orbital_densities_bin_dir, self.ML_dir, self.n_bands_to_extract[0], self.n_bands_to_extract[1], self.nspin_to_extract])
            utils.system_call(command)

    def compute_decomposition(self):
        print("compute decomposition")
        self.r_cut = min(self.atoms.get_cell_lengths_and_angles()[:3])/2.5 #the maximum radius has to be smaller than half of the cell-size
        print(self.r_cut)
        if self.method_to_extract_from_binary == 'from_ki':
            centers_occ = np.array(self.calculations[4].results['centers'])
            centers_emp = np.array(self.calculations[7].results['centers'])
            centers     = np.concatenate([centers_occ, centers_emp])
        
        
        ML_utils.precompute_radial_basis(self.n_max, self.l_max, self.r_min, self.r_max, self.ML_dir)
        ML_utils.func_compute_decomposition(self.n_max, self.l_max, self.r_min, self.r_max, self.r_cut, self.ML_dir, self.bands_to_extract, self.atoms, centers)
    
    def compute_power_spectrum(self):
        print("compute power spectrum")
        self.dir_power = self.ML_dir / f'power_spectra_{self.n_max}_{self.l_max}_{self.r_min}_{self.r_max}'
        ML_utils.main_compute_power(self.n_max, self.l_max, self.r_min, self.r_max, self.ML_dir, self.dir_power, self.bands_to_extract)

    def predict(self, band):
        power_spectrum = self.load_power_spectrum(band)
        y_predict = self.ml_model.predict(power_spectrum)
        print("y_predict = ", y_predict)
        return y_predict
    
    def train(self):
        print("Now training")
        self.ml_model.train()
    
    def add_training_data(self, band):
        power_spectrum = self.load_power_spectrum(band)
        alpha          = band.alpha
        print("adding orbital ", band.index)
        print("filling orbital ", band.filled)
        print("alpha orbital ", band.alpha)
        self.ml_model.add_training_data(power_spectrum, alpha)

    def load_power_spectrum(self, band):
        if band.filled:
            filled_str = 'occ'
        else:
            filled_str = 'emp'
        return np.loadtxt(self.dir_power / f"power_spectrum.orbital.{filled_str}.{band.index}.txt")


    
    
    def use_prediction(self):
        print("Use prediction -> False")
        return False