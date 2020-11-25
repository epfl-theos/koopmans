"""

Workflow module for python_KI, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import numpy as np
import os
import koopmans.utils as utils
from koopmans.calculators.calculator import run_qe
from koopmans.calculators.pw2wannier import PW2Wannier_calc
from koopmans.ase import set_up_pw2wannier
from ase.dft.kpoints import bandpath


def run(workflow_settings, calcs_dct):
    '''

    Wrapper for the calculation of (maximally localized) Wannier functions
    using PW and Wannier90

    '''

    print('\nWANNIZERIZATION')

    if workflow_settings["from_scratch"]:
        utils.system_call("rm -r wannier", False)

    if 'pw' not in calcs_dct:
        raise ValueError(
            'You need to provide a pw block in your input when init_variational_orbitals = "mlwfs"')
    if 'w90_occ' not in calcs_dct and 'w90_emp' not in calcs_dct:
        raise ValueError(
            'You need to provide a w90 block in your input when init_variational_orbitals = "mlwfs"')

    # Run PW scf and nscf calculations
    calc = calcs_dct['pw']
    calc.directory = 'wannier'
    calc.name = 'scf'
    # PWscf needs only the valence bands
    nbnd = calc.nbnd
    calc.nbnd = None
    run_qe(calc, silent=False)

    calc.name = 'nscf'
    calc.calculation = 'nscf'
    calc.nbnd = nbnd
    # Workaround to generate the right explicit MP mesh >>>
    kpts = np.indices(calc._ase_calc.parameters['kpts']).transpose(
        1, 2, 3, 0).reshape(-1, 3)
    kpts = kpts / calc._ase_calc.parameters['kpts']
    kpts = bandpath(kpts, calc._ase_calc.atoms.cell, npoints=len(kpts) - 1)
    calc._ase_calc.parameters['kpts'] = kpts
    # <<<
    run_qe(calc, silent=False)

    # Run Wannier90 and pw2wannier90 for occupied and empty states separately
    if 'pw2wannier' in calcs_dct:
        calc_p2w = calcs_dct['pw2wannier']
    else:
        calc_p2w = set_up_pw2wannier()

    w90_dir = {}
    for typ in ['occ', 'emp']:
        calc_w90 = calcs_dct['w90_' + typ]
        calc_w90.directory = 'wannier/' + typ
        calc_p2w.directory = calc_w90.directory
        w90_dir[typ] = calc_w90.directory

        # 1) pre-processing Wannier90 calculation
        calc_w90.name = 'wann'
        calc_w90.preprocessing_flags = '-pp'
        run_qe(calc_w90, silent=False)
        # 2) standard pw2wannier90 calculation
        calc_p2w.name = 'pw2wan'
        run_qe(calc_p2w, silent=False)
        # 3) Wannier90 calculation
        calc_w90.preprocessing_flags = ''
        run_qe(calc_w90, silent=False)
        # 4) pw2wannier90 calc (wannier2odd mode for WFs conversion to supercell)
        calc_w2o = set_up_pw2wannier()
        calc_w2o.directory = calc_w90.directory
        calc_w2o.name = 'wan2odd'
        calc_w2o.wan_mode = 'wannier2odd'
        # Workaround to run wannier2odd in serial >>>
#        if calc_w2o._ase_calc.command[:6] == 'mpirun':
#            calc_w2o._ase_calc.command = calc_w2o._ase_calc.command[13:]
#        elif calc_w2o._ase_calc.command[:4] == 'srun':
#            calc_w2o._ase_calc.command = calc_w2o._ase_calc.command[5:]
        # <<<
        run_qe(calc_w2o, silent=False)

    print()

    return w90_dir
