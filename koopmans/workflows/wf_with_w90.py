"""

Workflow module for python_KI, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import copy
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
                'You need to provide a pw block in your input when init_orbitals = "mlwfs" or "projwfs"')
        if 'w90_occ' not in self.master_calcs and 'w90_emp' not in self.master_calcs:
            raise ValueError(
                'You need to provide a w90 block in your input when init_orbitals = "mlwfs" or "projwfs"')

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
        calc_pw = self.new_calculator('pw', nbnd=None, nspin=1)
        calc_pw.directory = 'wannier'
        calc_pw.name = 'scf'
        self.run_calculator(calc_pw)

        calc_pw = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True, nspin=1)
        calc_pw.directory = 'wannier'
        calc_pw.name = 'nscf'
        self.run_calculator(calc_pw)

        for typ in ['occ', 'emp']:
            # 1) pre-processing Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ, name='wann_preproc',
                                           preprocessing_flags='-pp')
            self.run_calculator(calc_w90)
            utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

            # 2) standard pw2wannier90 calculation
            calc_p2w = self.new_calculator('pw2wannier', directory=calc_w90.directory, outdir=calc_pw.outdir,
                                           name='pw2wan')
            self.run_calculator(calc_p2w)

            # 3) Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ, name='wann')
            self.run_calculator(calc_w90)

            # 4) pw2wannier90 calc (wannier2odd mode for WFs conversion to supercell)
            if not skip_wann2odd:
                calc_w2o = self.new_calculator('pw2wannier', wan_mode='wannier2odd', directory=calc_w90.directory,
                                               name='wan2odd')
                if typ == 'emp':
                    calc_w2o.split_evc_file = True
                self.run_calculator(calc_w2o)

        print()

        return

    def new_calculator(self, calc_type, *args, **kwargs):
        calc = super().new_calculator(calc_type, *args, **kwargs)

        # Extra tweaks for Wannier90 calculations
        if calc_type.startswith('w90'):
            if calc.num_bands != calc.num_wann and calc.dis_num_iter is None:
                calc.dis_num_iter = 5000
            if self.init_orbitals == 'projwfs':
                calc.num_iter = 0

        # Checking that gamma_trick is consistent with do_wf_cmplx
        if calc_type == 'pw2wannier' and calc.wan_mode == 'wannier2odd':
            kcp_calc = self.master_calcs['kcp']
            if calc.gamma_trick == kcp_calc.do_wf_cmplx:
                utils.warn(
                    f'if do_wf_cmplx is {kcp_calc.do_wf_cmplx}, gamma_trick cannot be {calc.gamma_trick}. '
                    f'Changing gamma_trick to {not kcp_calc.do_wf_cmplx}')
                calc.gamma_trick = not kcp_calc.do_wf_cmplx
            elif calc.gamma_trick is None and not kcp_calc.do_wf_cmplx:
                calc.gamma_trick = True
            else:
                pass

        # Use a unified tmp directory
        if hasattr(calc, 'outdir'):
            calc.outdir = os.path.abspath('wannier/TMP')

        return calc
