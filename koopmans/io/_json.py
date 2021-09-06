"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

import os
import copy
import json as json_ext
import xml.etree.ElementTree as ET
from ase.atoms import Atoms
from ase.dft.kpoints import BandPath
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.espresso import Espresso, Espresso_kcp, KoopmansHam, KoopmansScreen, Wann2KC
from ase.calculators.espresso import PW2Wannier as ASEPW2Wannier
from ase.io.espresso.utils import ibrav_to_cell
from ase.io.espresso.pw import KEYS as pw_keys
from ase.io.espresso.koopmans_cp import KEYS as kcp_keys
from ase.io.wannier90 import Wannier90 as ASEWannier90
from koopmans import utils
from ._utils import read_atomic_species, read_atomic_positions, read_cell_parameters, read_kpoints_block, read_kpath, \
    set_up_pseudos, nelec_from_pseudos, construct_cell_parameters_block


def read_w90_dict(dct, generic_atoms):

    from koopmans import calculators

    # Setting up ASE calc object, copying over the generic atoms object and non-kcp-specific settings
    calc = ASEWannier90()
    atoms = copy.deepcopy(generic_atoms)
    atoms.calc.parameters.pop('input_data')
    atoms.calc.parameters.pop('pseudopotentials')
    calc.parameters = atoms.calc.parameters
    calc.atoms = atoms
    calc.atoms.calc = calc

    dct = {k.lower(): v for k, v in dct.items()}

    for k in dct.keys():
        if k in ['atoms_frac', 'unit_cell_cart', 'mp_grid']:
            raise ValueError('Please specify {k} in the "setup" block using kcp syntax rather than in the "w90" block')

    # Read in parameters
    calc.parameters.update(read_dict(dct))

    # Return koopmans-type calculator object rather than ASE calculator
    return calculators.Wannier90Calculator(calc)


def read_w90_occ_dict(dct, generic_atoms):

    calc = read_w90_dict(dct, generic_atoms)

    # Auto-generate values if they have not been provided
    n_filled = generic_atoms.calc.parameters['input_data']['system']['nelec'] // 2
    if calc.num_bands is None:
        calc.num_bands = n_filled
    if calc.num_wann is None:
        calc.num_wann = n_filled

    return calc


def read_w90_empty_dict(dct, generic_atoms):

    calc = read_w90_dict(dct, generic_atoms)

    # Auto-generate values if they have not been provided
    n_filled = generic_atoms.calc.parameters['input_data']['system']['nelec'] // 2
    n_empty = generic_atoms.calc.parameters['input_data']['electrons'].get('empty_states_nbnd', 0)

    if calc.num_wann is None:
        calc.num_wann = n_empty
    if calc.exclude_bands is None:
        calc.exclude_bands = f'1-{n_filled}'

    return calc


def read_pw2wannier_dict(dct, generic_atoms):

    from koopmans import calculators

    # Setting up ASE calc object
    calc = ASEPW2Wannier()

    calc.parameters['inputpp'] = read_dict(dct)

    # Attaching an atoms object (though pw2wannier doesn't need this information,
    # ASE will complain if it's missing)
    calc.atoms = copy.deepcopy(generic_atoms)
    calc.atoms.calc = calc

    # Giving the pw2wannier calculator access to the kpts info for reference
    if 'kpts' in generic_atoms.calc.parameters:
        calc.parameters['kpts'] = generic_atoms.calc.parameters['kpts']

    # Return koopmans-type calculator object rather than ASE calculator
    return calculators.PW2WannierCalculator(calc)


def read_setup_dict(dct):
    '''

    Reads the "setup" block. This block uses the same syntax as kcp

    '''

    calc = Espresso_kcp(atoms=Atoms())

    calc.parameters['input_data'] = {k: {} for k in pw_keys.keys()}

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
                calc.parameters['input_data'][block][key] = value
        else:
            raise ValueError(f'Unrecognised block "setup:{block}" in the input file')

    # Calculating the simulation cell
    cell = None
    if 'cell_parameters' in dct:
        subdct = dct['cell_parameters']
        cell = read_cell_parameters(calc, subdct)

    # Generating cell if it is missing
    if cell is None:
        _, cell = ibrav_to_cell(calc.parameters['input_data']['system'])

    # Attaching the cell to the calculator
    calc.atoms = Atoms(cell=cell)

    # Calculating kpoints
    if 'k_points' in dct:
        read_kpoints_block(calc, dct['k_points'])

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
    if 'nelec' not in calc.parameters['input_data']['system']:
        nelec = nelec_from_pseudos(calc)
        calc.parameters['input_data']['system']['nelec'] = nelec
    else:
        nelec = calc.parameters['input_data']['system']['nelec']
    if 'tot_charge' in calc.parameters['input_data']['system']:
        tot_charge = calc.parameters['input_data']['system']['tot_charge']
        calc.parameters['input_data']['system']['nelec'] -= tot_charge
        nelec -= tot_charge
    if 'tot_magnetization' in calc.parameters['input_data']['system']:
        tot_mag = calc.parameters['input_data']['system']['tot_magnetization']
    else:
        tot_mag = nelec % 2
        calc.parameters['input_data']['system']['tot_magnetization'] = tot_mag
    if 'nelup' not in calc.parameters['input_data']['system']:
        calc.parameters['input_data']['system']['nelup'] = int(
            nelec / 2 + tot_mag / 2)
    if 'neldw' not in calc.parameters['input_data']['system']:
        calc.parameters['input_data']['system']['neldw'] = int(
            nelec / 2 - tot_mag / 2)

    return calc.atoms


def read_kcp_or_pw_dict(dct, calc):
    '''

    Reads in dict of kcp/pw input file and adds the settings to the provided calc object

    '''

    if isinstance(calc, Espresso_kcp):
        qe_keys = kcp_keys
    elif isinstance(calc, Espresso):
        qe_keys = pw_keys
    else:
        raise TypeError(f'io.read_kcp_or_pw_dict() does not accept "calc" with class {calc.__class__}')

    skipped_blocks = []
    for block, subdct in dct.items():
        if block in qe_keys:
            if block not in calc.parameters['input_data']:
                calc.parameters['input_data'][block] = {}

            for key, value in subdct.items():
                if value == "":
                    continue

                if key == 'pseudo_dir':
                    raise ValueError('Please specify "pseudo_dir" in the "setup" control block')

                if key in calc.parameters['input_data'][block]:
                    utils.warn(f'Overwriting value for {block}:{key} provided in the "setup" block')

                try:
                    value = json_ext.loads(value)
                except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                    pass
                calc.parameters['input_data'][block][key] = value
        else:
            skipped_blocks.append(block)

    return skipped_blocks


def read_kcp_dict(dct, generic_atoms):
    '''

    Reads in dict of kcp input file and returns a KoopmansCPCalculator object

    '''

    from koopmans import calculators

    # Copy over the settings from the "setup" block
    generic_atoms_copy = copy.deepcopy(generic_atoms)
    calc = generic_atoms_copy.calc
    calc.atoms.calc = calc

    # Read in any settings provided in the "kcp" block
    skipped_blocks = read_kcp_or_pw_dict(dct, calc)

    for block in skipped_blocks:
        utils.warn(f'The {block} block is not yet implemented and will be ignored')

    return calculators.KoopmansCPCalculator(calc)


def read_pw_dict(dct, generic_atoms):
    '''

    Reads in dict of pw input file and returns a PWCalculator object

    '''

    from koopmans import calculators

    # Initialising a pw calculator tethered to a copy of the atoms object
    calc = Espresso()
    calc.atoms = copy.deepcopy(generic_atoms)
    calc.atoms.calc = calc

    # Remove settings using kcp-syntax
    kcp_settings = generic_atoms.calc.parameters.pop('input_data')

    # Copy over parameters that use generic syntax (pseudos, k_points)
    calc.parameters = generic_atoms.calc.parameters

    # Initialise the pw-specific settings
    calc.parameters['input_data'] = {key: {} for key in pw_keys}

    # Convert kcp-syntax settings to pw-syntax settings
    for block, subdct in kcp_settings.items():
        for key, val in subdct.items():
            if key in ['nelec', 'nelup', 'neldw', 'empty_states_nbnd', 'tot_magnetization']:
                continue
            elif block in pw_keys and key in pw_keys[block]:
                # PW and KCP share this keyword so we can copy it over directly
                calc.parameters['input_data'][block][key] = val
            elif key == 'conv_thr':
                # Pw uses Ry, KCP uses Ha = 2 Ry
                calc.parameters['input_data'][block][key] = val * 2
            else:
                raise ValueError(f'Could not convert {block}:{key} to a pw keyword')

    # Read the content of the pw block
    skipped_blocks = read_kcp_or_pw_dict(dct, calc)
    for block in skipped_blocks:
        if block == 'k_points':
            read_kpoints_block(calc, dct[block])
        else:
            utils.warn(f'The {block} block is not yet implemented and will be ignored')

    # If no nbnd is provided, auto-generate it
    if 'nbnd' not in calc.parameters['input_data']['system']:
        n_elec = kcp_settings['system']['nelec']
        n_empty = kcp_settings['electrons'].get('empty_states_nbnd', 0)
        calc.parameters['input_data']['system']['nbnd'] = n_elec // 2 + n_empty

    return calculators.PWCalculator(calc)


def read_ui_dict(dct, generic_atoms):

    from koopmans import calculators

    # For UI, just use a generic calculator with no command
    atoms = copy.deepcopy(generic_atoms)
    calc = atoms.calc
    calc.atoms.calc = calc

    calc.command = ''

    # Overwrite the parameters with the provided JSON dict
    calc.parameters = read_dict(dct)

    # Use kpath specified in the setup block if it is not present in the ui dict
    setup_kpath = generic_atoms.calc.parameters.get('kpath', None)
    if 'kpath' not in dct and setup_kpath:
        calc.parameters['kpath'] = setup_kpath

    # Convert units of alat_sc
    if 'alat_sc' in calc.parameters:
        calc.parameters['alat_sc'] *= utils.units.Bohr

    return calculators.UnfoldAndInterpolateCalculator(calc)


def read_kc_wann_dict(dct, generic_atoms):

    from koopmans.calculators import KoopmansHamCalculator, KoopmansScreenCalculator, Wann2KCCalculator

    calcs = {}
    for key, ase_calc_class, calc_class in (('kc_ham', KoopmansHam, KoopmansHamCalculator),
                                            ('kc_screen', KoopmansScreen, KoopmansScreenCalculator),
                                            ('wann2kc', Wann2KC, Wann2KCCalculator)):
        # Create a copy of generic_atoms
        atoms = copy.deepcopy(generic_atoms)

        # Create an ASE calculator object and copy over the settings
        calc = ase_calc_class(atoms=atoms)
        calc.parameters = copy.deepcopy(generic_atoms.calc.parameters)

        # Read in the parameters
        relevant_subblocks = ['control', 'system']
        if key == 'kc_ham':
            relevant_subblocks.append('ham')
            if 'ham' not in dct.keys():
                if 'kpts' not in calc.parameters:
                    # Populating kpts if absent
                    calc.parameters['kpts'] = [1, 1, 1]
                    calc.atoms.cell.pbc = True
                    generic_atoms.calc.parameters['kpath'] = BandPath(calc.atoms.cell, [[0, 0, 0]])
                    calc.atoms.cell.pbc = False
                dct['ham'] = {}

            # kc_ham stores the kpoint properties slightly differently
            # kpts -> mp1-3
            for i, nk in enumerate(calc.parameters['kpts']):
                dct['ham'][f'mp{i+1}'] = nk

            # kpath -> kpts
            if 'kpath' in dct['ham']:
                kpath = dct['ham'].pop('kpath')
                read_kpath(calc, kpath)
                calc.parameters['kpts'] = calc.parameters['kpath']
            elif 'kpath' in generic_atoms.calc.parameters:
                calc.parameters['kpts'] = generic_atoms.calc.parameters['kpath']
            else:
                continue

        elif key == 'kc_screen':
            relevant_subblocks.append('screen')

        flattened_settings = {k: v for block_name, block in dct.items() for k, v in block.items()
                              if block_name in relevant_subblocks}

        # Convert to a ExtendedCalculator class, loading the settings
        calcs[key] = calc_class(calc, **flattened_settings)

    # Unlike most read_*_dict functions, this returns three calculators in a dict
    return calcs


def read_dict(dct):
    '''

    Reads in a dict, formatting the values appropriately if they are not already

    '''
    settings = {}
    for k, v in dct.items():
        # Deal with bools separately since JSON strictly only will interpret
        # 'false' as False, while 'False' will be left as a string and
        # any statement to the effect of 'if param' will evaluate to True if
        # param = 'False'
        if isinstance(v, str) and v.lower() in ['f', 'false']:
            settings[k] = False
        elif isinstance(v, str) and v.lower() in ['t', 'true']:
            settings[k] = True
        else:
            try:
                settings[k] = json_ext.loads(v)
            except (TypeError, json_ext.decoder.JSONDecodeError) as e:
                settings[k] = v
    return settings


def update_nested_dict(dct_to_update, second_dct):
    for k, v in second_dct.items():
        if k in dct_to_update and isinstance(v, dict):
            update_nested_dict(dct_to_update[k], second_dct[k])
        else:
            dct_to_update[k] = v


def read_json(fd, override={}):
    '''

    Reads in settings listed in JSON file

    Values in the JSON file can be overridden by values provided in the override argument

    '''

    from koopmans.calculators import KoopmansCPCalculator
    from koopmans.workflows.singlepoint import SinglepointWorkflow
    from koopmans.workflows.convergence import ConvergenceWorkflow
    from koopmans.workflows.pbe_dscf_with_pw import DeltaSCFWorkflow
    from koopmans.workflows.ui import UnfoldAndInterpolateWorkflow
    from koopmans.workflows.wf_with_w90 import WannierizeWorkflow

    if isinstance(fd, str):
        fd = open(fd, 'r')

    bigdct = json_ext.loads(fd.read())

    # Override all keywords provided explicitly
    update_nested_dict(bigdct, override)

    # Deal with w90 subdicts
    if 'w90' in bigdct:
        bigdct['w90_occ'] = bigdct['w90'].pop('occ', {})
        bigdct['w90_emp'] = bigdct['w90'].pop('emp', {})
        for k, v in bigdct['w90'].items():
            bigdct['w90_occ'][k] = v
            bigdct['w90_emp'][k] = v
        del bigdct['w90']

    # Deal with UI subdicts
    if 'ui' in bigdct:
        subdcts = {}
        keys = ['occ', 'emp']
        for key in keys:
            # First, we must remove the occ and emp subdicts from the UI dict
            if key in bigdct['ui']:
                subdcts[key] = bigdct['ui'].pop(key)

        # Now, we add the ui_occ and ui_emp calculators to master_calcs
        for key in keys:
            if key in subdcts:
                # Add the corresponding subdict to the rest of the UI block
                bigdct[f'ui_{key}'] = dict(bigdct['ui'], **subdcts[key])
            else:
                # Duplicate the UI block
                bigdct[f'ui_{key}'] = bigdct['ui']

    # Define which function to use to read each block
    readers = {'kcp': read_kcp_dict, 'w90_occ': read_w90_occ_dict, 'w90_emp': read_w90_empty_dict,
               'pw2wannier': read_pw2wannier_dict, 'pw': read_pw_dict, 'ui': read_ui_dict,
               'ui_occ': read_ui_dict, 'ui_emp': read_ui_dict, 'kc_wann': read_kc_wann_dict}

    # Check for unexpected blocks
    for block in bigdct:
        if block not in list(readers.keys()) + ['workflow', 'setup']:
            raise ValueError(f'Unrecognised block "{block}" in json input file; '
                             'valid options are workflow/' + '/'.join(readers.keys()))

    # Loading workflow settings
    workflow_settings = read_dict(bigdct.get('workflow', {}))
    task_name = workflow_settings.pop('task', 'singlepoint')

    # Load default values
    if 'setup' in bigdct:
        generic_atoms = read_setup_dict(bigdct['setup'])
        del bigdct['setup']
    elif task_name != 'ui':
        raise ValueError('You must provide a "setup" block in the input file, specifying atomic positions, atomic '
                         'species, etc.')
    else:
        # Create an empty dummy ASE calculator to attach settings to
        generic_atoms = Atoms(calculator=FileIOCalculator())
        generic_atoms.calc.atoms = generic_atoms

    # Loading calculator-specific settings
    calcs_dct = {}

    # Generate a master calculator for every single kind of calculator, regardless of whether or not there was a
    # corresponding block in the json file
    for block, reader in readers.items():
        # For the UI task, the input file won't contain enough for us to generate all the other kinds of calculators
        if task_name == 'ui' and block != 'ui':
            if block in bigdct and not block.startswith('ui'):
                utils.warn(f'Ignoring the {block} block since "task": "ui"')
            continue

        # Read the block and add the resulting calculator to the calcs_dct
        dct = bigdct.get(block, {})
        if block == 'kc_wann':
            calcs_dct.update(readers[block](dct, generic_atoms))
        else:
            calcs_dct[block] = readers[block](dct, generic_atoms)

    name = fd.name.replace('.json', '')
    if task_name == 'singlepoint':
        workflow = SinglepointWorkflow(workflow_settings, calcs_dct, name)
    elif task_name == 'convergence':
        workflow = ConvergenceWorkflow(workflow_settings, calcs_dct, name)
    elif task_name in ['wannierize', 'wannierise']:
        workflow = WannierizeWorkflow(workflow_settings, calcs_dct, name, check_wannierisation=True)
    elif task_name == 'environ_dscf':
        workflow = DeltaSCFWorkflow(workflow_settings, calcs_dct, name)
    elif task_name == 'ui':
        workflow = UnfoldAndInterpolateWorkflow(workflow_settings, calcs_dct, name)
    else:
        raise ValueError('Invalid task name "{task_name}"')

    return workflow


def write_json(workflow, filename):
    '''

    Writes out settings to a JSON file

    '''

    from koopmans.workflows.generic import valid_settings
    from koopmans.calculators import PWCalculator, KoopmansCPCalculator

    fd = open(filename, 'w')

    calcs = workflow.master_calcs

    bigdct = {}

    # "workflow" block (not printing any values that match the defaults)
    bigdct['workflow'] = {}
    for k, v in workflow.settings.items():
        if v is None:
            continue
        setting = [s for s in valid_settings if s.name == k][0]
        if v != setting.default:
            bigdct['workflow'][k] = v

    # "setup" block
    # Working out ibrav
    kcp_calc = calcs.get('kcp', None)
    pw_calc = calcs.get('pw', None)
    if kcp_calc:
        calc = kcp_calc.calc
        ibrav = calc.parameters['input_data']['system'].get('ibrav', 0)
    elif pw_calc:
        calc = pw_calc
        ibrav = calc.parameters['input_data']['system'].get('ibrav', 0)
    else:
        calc = calcs.values().calc
        ibrav = 0

    bigdct['setup'] = {}

    # cell parameters
    if ibrav == 0:
        bigdct['setup']['cell_parameters'] = construct_cell_parameters_block(calc)

    # atomic positions
    if calc.atoms.has('labels'):
        labels = calc.atoms.get_array('labels')
    else:
        labels = calc.atoms.get_chemical_symbols()
    if ibrav == 0:
        bigdct['setup']['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_positions())],
            'units': 'angstrom'}
    else:
        bigdct['setup']['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, calc.atoms.get_scaled_positions())],
            'units': 'alat'}

    # atomic species
    bigdct['setup']['atomic_species'] = {'species': [[key, 1.0, val]
                                                     for key, val in calc.parameters.pseudopotentials.items()]}

    for code, calc in calcs.items():
        if isinstance(calc, (KoopmansCPCalculator, PWCalculator)):
            bigdct[code] = {}
            calc = calc.calc

            # pseudo directory
            pseudo_dir = calc.parameters['input_data']['control'].pop('pseudo_dir', None)
            if pseudo_dir is not None:
                bigdct['setup']['control'] = {'pseudo_dir': pseudo_dir}

            # Populate bigdct with the settings
            for key, block in calc.parameters['input_data'].items():
                if len(block) > 0:
                    bigdct[code][key] = {k: v for k, v in dict(
                        block).items() if v is not None}

            if code == 'pw':
                # Adding kpoints to "setup"
                if 'kpts' in calc.parameters:
                    kpts = {'kind': 'automatic',
                            'kpts': calc.parameters['kpts'],
                            'koffset': calc.parameters['koffset']}
                else:
                    kpts = {'kind': 'gamma'}
                bigdct['setup']['k_points'] = kpts
        else:
            raise ValueError(
                f'Writing of {calc.__class__} with write_json is not yet implemented')

    json_ext.dump(bigdct, fd, indent=2)

    fd.close()

    return
