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
                'You need to provide a pw block in your input when init_variational_orbitals = "mlwfs"')
        if 'w90_occ' not in self.master_calcs and 'w90_emp' not in self.master_calcs:
            raise ValueError(
                'You need to provide a w90 block in your input when init_variational_orbitals = "mlwfs"')

        if 'pw2wannier' not in self.master_calcs:
            # Generate a pw2wannier calculator with default settings
            p2w_calc = PW2Wannier_calc()
            # Copying over the atoms object from the PW calculator
            p2w_calc.calc.atoms = self.master_calcs['pw'].calc.atoms
            p2w_calc.calc.atoms.calc = p2w_calc.calc
            self.master_calcs['pw2wannier'] = p2w_calc

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        print('\nWANNIZERIZATION')

        if self.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc = self.new_calculator('pw', nbnd=None)
        calc.directory = 'wannier'
        calc.name = 'scf'
        self.run_calculator(calc)

        calc = self.new_calculator('pw', calculation='nscf')
        calc.directory = 'wannier'
        calc.name = 'nscf'
        # Workaround to generate the right explicit MP mesh >>>
        kpts = np.indices(calc.calc.parameters['kpts']).transpose(
            1, 2, 3, 0).reshape(-1, 3)
        kpts = kpts / calc.calc.parameters['kpts']
        kpts = bandpath(kpts, calc.calc.atoms.cell, npoints=len(kpts) - 1)
        calc.calc.parameters['kpts'] = kpts
        # <<<
        self.run_calculator(calc)

        calc_p2w = self.new_calculator('pw2wannier')
        for typ in ['occ', 'emp']:
            calc_w90 = self.new_calculator('w90_' + typ)
            calc_w90.directory = 'wannier/' + typ
            calc_p2w.directory = calc_w90.directory

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
            calc_w2o = self.new_calculator('pw2wannier', wan_mode='wannier2odd')
            calc_w2o.directory = calc_w90.directory
            calc_w2o.name = 'wan2odd'
            if typ == 'emp':
                calc_w2o.empty_states = True
            # Workaround to run wannier2odd in serial >>>
            # if calc_w2o.calc.command[:6] == 'mpirun':
            #     calc_w2o.calc.command = calc_w2o.calc.command[13:]
            # elif calc_w2o.calc.command[:4] == 'srun':
            #     calc_w2o.calc.command = calc_w2o.calc.command[5:]
            # <<<
            self.run_calculator(calc_w2o)

        print()

        return
