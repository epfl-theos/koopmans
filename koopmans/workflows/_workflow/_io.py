"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

import os
import json as json_ext
from ase.atoms import Atoms
from ase.io.espresso.utils import ibrav_to_cell
from ase.calculators.espresso import Espresso_kcp
from ase.io.espresso.koopmans_cp import KEYS as kcp_keys
from koopmans.pseudopotentials import set_up_pseudos, nelec_from_pseudos
from koopmans.utils import read_atomic_species, read_atomic_positions, read_cell_parameters, read_kpoints_block


def read_setup_dict(dct, task):
    '''

    Reads the "setup" block. This block uses the same syntax as kcp

    '''

    calc = Espresso_kcp(atoms=Atoms())

    compulsory_block_readers = {'atomic_species': read_atomic_species,
                                'atomic_positions': read_atomic_positions}

    for block, subdct in dct.items():
        if block in compulsory_block_readers or block in ['cell_parameters', 'k_points']:
            # We will read these afterwards
            continue
        elif block in kcp_keys:
            for key, value in subdct.items():
                if value == "":
                    continue

                # Force pseudo_dir to be an absolute path
                if key == 'pseudo_dir' and value[0] != '/':
                    value = os.path.abspath(value) + '/'

                try:
                    value = json_ext.loads(value)
                except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                    pass
                calc.parameters[key] = value
        else:
            raise ValueError(f'Unrecognised block "setup:{block}" in the input file')

    # Calculating the simulation cell
    cell = None
    if 'cell_parameters' in dct:
        subdct = dct['cell_parameters']
        cell = read_cell_parameters(calc, subdct)

    # Generating cell if it is missing
    if cell is None:
        _, cell = ibrav_to_cell(calc.parameters)

    # Attaching the cell to the calculator
    calc.atoms = Atoms(cell=cell)

    # Calculating kpoints
    if 'k_points' in dct:
        read_kpoints_block(calc, dct['k_points'])

    if task != 'ui':
        def read_compulsory_block(block_name, extract_function):
            if block_name in dct:
                subdct = dct[block_name]
                extract_function(calc, subdct)
                del dct[block_name]
            else:
                raise ValueError(f'{block_name} not found in "setup" block')

        for block_name, extract_function in compulsory_block_readers.items():
            read_compulsory_block(block_name, extract_function)

        set_up_pseudos(calc)

        # If they are missing, fill the nelec/nelup/neldw field using information contained
        # in the pseudopotential files
        if 'nelec' not in calc.parameters:
            nelec = nelec_from_pseudos(calc.atoms, calc.parameters.pseudopotentials, calc.parameters.pseudo_dir)
            calc.parameters.nelec = nelec
        else:
            nelec = calc.parameters.nelec

        if calc.parameters.get('nspin', 2) == 2:
            if 'tot_charge' in calc.parameters:
                tot_charge = calc.parameters.tot_charge
                calc.parameters.nelec -= tot_charge
                nelec -= tot_charge
            if 'tot_magnetization' in calc.parameters:
                tot_mag = calc.parameters.tot_magnetization
            else:
                tot_mag = nelec % 2
                calc.parameters.tot_magnetization = tot_mag
            if 'nelup' not in calc.parameters:
                calc.parameters.nelup = int(nelec / 2 + tot_mag / 2)
            if 'neldw' not in calc.parameters:
                calc.parameters.neldw = int(nelec / 2 - tot_mag / 2)
            if tot_mag != 0:
                if 'starting_magnetization(1)' not in calc.parameters:
                    calc.atoms.set_initial_magnetic_moments([tot_mag / len(calc.atoms) for _ in calc.atoms])

        # Work out the number of filled and empty bands
        n_filled = nelec // 2 + nelec % 2
        n_empty = calc.parameters.get('empty_states_nbnd', 0)
    else:
        # these parameters require a definition but they are not necessary for a UI workflow
        calc.parameters.nelec = 0
        n_filled = 0
        n_empty = 0

    # Separamting the output into atoms, parameters, and psp+kpoint information
    atoms = calc.atoms
    atoms.calc = None
    parameters = calc.parameters
    psps_and_kpts = {}
    for key in ['pseudopotentials', 'gamma_only', 'kgrid', 'koffset', 'kpath']:
        if key in parameters:
            psps_and_kpts[key] = parameters.pop(key)

    return atoms, parameters, psps_and_kpts, n_filled, n_empty
