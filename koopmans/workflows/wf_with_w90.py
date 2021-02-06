"""

Workflow module for python_KI, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import os
import numpy as np
from koopmans import utils
from koopmans.calculators.pw2wannier import PW2Wannier_calc
from ase.dft.kpoints import bandpath
from koopmans.workflows.generic import Workflow


class WannierizeWorkflow(Workflow):

    def __init__(self, workflow_settings, calcs_dct):
        super().__init__(workflow_settings, calcs_dct)

        if 'pw' not in self.master_calcs:
            raise ValueError(
                'You need to provide a pw block in your input when init_variational_orbitals = "mlwfs" or "projw"')
        if 'w90_occ' not in self.master_calcs and 'w90_emp' not in self.master_calcs:
            raise ValueError(
                'You need to provide a w90 block in your input when init_variational_orbitals = "mlwfs" or "projw"')

        if 'pw2wannier' not in self.master_calcs:
            # Generate a pw2wannier calculator with default settings
            p2w_calc = PW2Wannier_calc()
            # Copying over the atoms object from the PW calculator
            p2w_calc.calc.atoms = self.master_calcs['pw'].calc.atoms
            p2w_calc.calc.atoms.calc = p2w_calc.calc
            self.master_calcs['pw2wannier'] = p2w_calc

        # Make sure num_wann (occ/empty), num_bands (occ/empty), and nbnd are present and consistent
        pw_calc = self.master_calcs['pw']
        w90_occ_calc = self.master_calcs['w90_occ']
        w90_emp_calc = self.master_calcs['w90_emp']

        # Filling out missing fields
        if w90_emp_calc.num_bands is None:
            w90_emp_calc.num_bands = pw_calc.nbnd - w90_occ_calc.num_bands
        if w90_occ_calc.exclude_bands is None:
            w90_occ_calc.exclude_bands = f'{w90_occ_calc.num_bands + 1}-{pw_calc.nbnd}'

        # Sanity checking
        assert pw_calc.nbnd == w90_occ_calc.num_bands + w90_emp_calc.num_bands
        if w90_emp_calc.num_wann == 0:
            raise ValueError('Cannot run a wannier90 calculation with num_wann = 0. Please set empty_states_nbnd > 0 '
                             'in the setup block, or num_wann > 0 in the wannier90 empty subblock')

    def run(self, skip_wann2odd=False):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        print('\nWANNIZERIZATION')

        if self.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc = self.new_calculator('pw', nbnd=None, nspin=1)
        calc.directory = 'wannier'
        calc.name = 'scf'
        self.run_calculator(calc)

        calc = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True, nspin=1)
        calc.directory = 'wannier'
        calc.name = 'nscf'
        self.run_calculator(calc)

        calc_p2w = self.new_calculator('pw2wannier')
        for typ in ['occ', 'emp']:
            calc_w90 = self.new_calculator('w90_' + typ)

            calc_w90.directory = 'wannier/' + typ
            calc_p2w.directory = calc_w90.directory
            if calc_w90.num_bands != calc_w90.num_wann and calc_w90.dis_num_iter is None:
                calc_w90.dis_num_iter = 5000
            if self.init_variational_orbitals == 'projw':
                calc_w90.num_iter = 0

            # 1) pre-processing Wannier90 calculation
            calc_w90.name = 'wann'
            calc_w90.preprocessing_flags = '-pp'
            self.run_calculator(calc_w90)
            # 2) standard pw2wannier90 calculation
            calc_p2w.name = 'pw2wan'
            self.run_calculator(calc_p2w)
            # 3) Wannier90 calculation
            calc_w90.preprocessing_flags = ''
            self.run_calculator(calc_w90)
            # 4) pw2wannier90 calc (wannier2odd mode for WFs conversion to supercell)
            if not skip_wann2odd:
                calc_w2o = self.new_calculator('pw2wannier', wan_mode='wannier2odd')
                calc_w2o.directory = calc_w90.directory
                calc_w2o.name = 'wan2odd'
                if typ == 'emp':
                    calc_w2o.split_evc_file = True
                self.run_calculator(calc_w2o)

        print()

        return
