"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

from koopmans import pseudopotentials
from koopmans.utils import read_atomic_species, read_atomic_positions, read_cell_parameters, read_kpoints_block, read_kpath, \
    construct_cell_parameters_block
from koopmans import workflows
from koopmans import calculators
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
from koopmans.pseudopotentials import set_up_pseudos, nelec_from_pseudos
from koopmans.settings import KoopmansCPSettingsDict, KoopmansHamSettingsDict, KoopmansScreenSettingsDict, PWSettingsDict, PW2WannierSettingsDict, UnfoldAndInterpolateSettingsDict, Wann2KCSettingsDict, Wannier90SettingsDict


def read_setup_dict(dct):
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
        nelec = nelec_from_pseudos(calc)
        calc.parameters.nelec = nelec
    else:
        nelec = calc.parameters.nelec
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

    # Work out the number of filled and empty bands
    n_filled = nelec // 2
    n_empty = calc.parameters.get('empty_states_nbnd', 0)

    # Separamting the output into atoms, parameters, and psp+kpoint information
    atoms = calc.atoms
    atoms.calc = None
    parameters = calc.parameters
    psps_and_kpts = {}
    for key in ['pseudopotentials', 'kpts', 'koffset', 'kpath']:
        if key in parameters:
            psps_and_kpts[key] = parameters.pop(key)

    return atoms, parameters, psps_and_kpts, n_filled, n_empty


# def read_kcp_or_pw_dict(dct, calc):
#     '''
#
#     Reads in dict of kcp/pw input file and adds the settings to the provided calc object
#
#     '''
#
#     if isinstance(calc, Espresso_kcp):
#         qe_keys = kcp_keys
#     elif isinstance(calc, Espresso):
#         qe_keys = pw_keys
#     else:
#         raise TypeError(f'io.read_kcp_or_pw_dict() does not accept "calc" with class {calc.__class__}')
#
#     skipped_blocks = []
#     for block, subdct in dct.items():
#         if block in qe_keys:
#             if block not in calc.parameters['input_data']:
#                 calc.parameters['input_data'][block] = {}
#
#             for key, value in subdct.items():
#                 if value == "":
#                     continue
#
#                 if key == 'pseudo_dir':
#                     raise ValueError('Please specify "pseudo_dir" in the "setup" control block')
#
#                 if key in calc.parameters['input_data'][block]:
#                     utils.warn(f'Overwriting value for {block}:{key} provided in the "setup" block')
#
#                 try:
#                     value = json_ext.loads(value)
#                 except (TypeError, json_ext.decoder.JSONDecodeError) as e:
#                     pass
#                 calc.parameters['input_data'][block][key] = value
#         else:
#             skipped_blocks.append(block)
#
#     return skipped_blocks
#
#
# def read_kcp_dict(dct, generic_atoms):
#     '''
#
#     Reads in dict of kcp input file and returns a KoopmansCPCalculator object
#
#     '''
#
#     from koopmans import calculators
#
#     # Copy over the settings from the "setup" block
#     generic_atoms_copy = copy.deepcopy(generic_atoms)
#     calc = generic_atoms_copy.calc
#     calc.atoms.calc = calc
#
#     # Read in any settings provided in the "kcp" block
#     skipped_blocks = read_kcp_or_pw_dict(dct, calc)
#
#     for block in skipped_blocks:
#         utils.warn(f'The {block} block is not yet implemented and will be ignored')
#
#     return calculators.KoopmansCPCalculator(calc)
#
#
# def read_pw_dict(dct, generic_atoms):
#     '''
#
#     Reads in dict of pw input file and returns a PWCalculator object
#
#     '''
#
#     from koopmans import calculators
#
#     # Initialising a pw calculator tethered to a copy of the atoms object
#     calc = Espresso()
#     calc.atoms = copy.deepcopy(generic_atoms)
#     calc.atoms.calc = calc
#
#     # Remove settings using kcp-syntax
#     kcp_settings = copy.deepcopy(generic_atoms.calc.parameters)
#
#     # Copy over parameters that use generic syntax
#     for key in ['pseudopotentials', 'kpts', 'kpath']:
#         if key in kcp_settings:
#             calc.parameters[key] = kcp_settings.pop(key)
#
#     # Convert kcp-syntax settings to pw-syntax settings
#     for key, val in kcp_settings.items():
#         if key in ['nelec', 'nelup', 'neldw', 'empty_states_nbnd', 'tot_magnetization']:
#             continue
#         elif key in [k for block in pw_keys.values() for k in block]:
#             # PW and KCP share this keyword so we can copy it over directly
#             calc.parameters[key] = val
#         elif key == 'conv_thr':
#             # PW uses Ry, KCP uses Ha = 2 Ry
#             calc.parameters[key] = val * 2
#         else:
#             raise ValueError(f'Could not convert {key} to a pw keyword')
#
#     # Read the content of the pw block
#     skipped_blocks = read_kcp_or_pw_dict(dct, calc)
#     for block in skipped_blocks:
#         if block == 'k_points':
#             read_kpoints_block(calc, dct[block])
#         else:
#             utils.warn(f'The {block} block is not yet implemented and will be ignored')
#
#     # If no nbnd is provided, auto-generate it
#     if 'nbnd' not in calc.parameters:
#         n_elec = kcp_settings.nelec
#         n_empty = kcp_settings.get('empty_states_nbnd', 0)
#         calc.parameters.nbnd = n_elec // 2 + n_empty
#
#     return calculators.PWCalculator(calc)
#
#
# def read_kc_wann_dict(dct, generic_atoms):
#
#     from koopmans.calculators import KoopmansHamCalculator, KoopmansScreenCalculator, Wann2KCCalculator
#
#     calcs = {}
#     for key, ase_calc_class, calc_class in (('kc_ham', KoopmansHam, KoopmansHamCalculator),
#                                             ('kc_screen', KoopmansScreen, KoopmansScreenCalculator),
#                                             ('wann2kc', Wann2KC, Wann2KCCalculator)):
#         # Create a copy of generic_atoms
#         atoms = copy.deepcopy(generic_atoms)
#
#         # Create an ASE calculator object and copy over the settings
#         calc = ase_calc_class(atoms=atoms)
#         calc.parameters = copy.deepcopy(generic_atoms.calc.parameters)
#
#         # Read in the parameters
#         relevant_subblocks = ['control', 'system']
#         if key == 'kc_ham':
#             relevant_subblocks.append('ham')
#             if 'ham' not in dct.keys():
#                 if 'kpts' not in calc.parameters:
#                     # Populating kpts if absent
#                     calc.parameters['kpts'] = [1, 1, 1]
#                     calc.atoms.cell.pbc = True
#                     generic_atoms.calc.parameters['kpath'] = BandPath(calc.atoms.cell, [[0, 0, 0]])
#                     calc.atoms.cell.pbc = False
#                 dct['ham'] = {}
#
#             # kc_ham stores the kpoint properties slightly differently
#             # kpts -> mp1-3
#             for i, nk in enumerate(calc.parameters['kpts']):
#                 dct['ham'][f'mp{i+1}'] = nk
#
#             # kpath -> kpts
#             if 'kpath' in dct['ham']:
#                 kpath = dct['ham'].pop('kpath')
#                 read_kpath(calc, kpath)
#                 calc.parameters['kpts'] = calc.parameters['kpath']
#             elif 'kpath' in generic_atoms.calc.parameters:
#                 calc.parameters['kpts'] = generic_atoms.calc.parameters['kpath']
#             else:
#                 continue
#
#         elif key == 'kc_screen':
#             relevant_subblocks.append('screen')
#
#         flattened_settings = {k: v for block_name, block in dct.items() for k, v in block.items()
#                               if block_name in relevant_subblocks}
#
#         # Convert to a CalculatorExt class, loading the settings
#         calcs[key] = calc_class(calc, **flattened_settings)
#
#     # Unlike most read_*_dict functions, this returns three calculators in a dict
#     return calcs


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

        # Now, we add the ui_occ and ui_emp calculators to master_calc_params
        for key in keys:
            if key in subdcts:
                # Add the corresponding subdict to the rest of the UI block
                bigdct[f'ui_{key}'] = dict(bigdct['ui'], **subdcts[key])
            else:
                # Duplicate the UI block
                bigdct[f'ui_{key}'] = bigdct['ui']

    # Deal with kc_wann subdicts
    kc_wann_blocks = bigdct.pop('kc_wann', {'kc_ham': {}, 'kc_screen': {}, 'wann2kc': {}})
    bigdct.update(**kc_wann_blocks)

    # Define which function to use to read each block
    settings_classes = {'kcp': KoopmansCPSettingsDict,
                        'kc_ham': KoopmansHamSettingsDict,
                        'kc_screen': KoopmansScreenSettingsDict,
                        'wann2kc': Wann2KCSettingsDict,
                        'pw': PWSettingsDict,
                        'pw2wannier': PW2WannierSettingsDict,
                        'ui': UnfoldAndInterpolateSettingsDict,
                        'ui_occ': UnfoldAndInterpolateSettingsDict,
                        'ui_emp': UnfoldAndInterpolateSettingsDict,
                        'w90_occ': Wannier90SettingsDict,
                        'w90_emp': Wannier90SettingsDict}

    # Check for unexpected blocks
    for block in bigdct:
        if block not in list(settings_classes.keys()) + ['workflow', 'setup']:
            raise ValueError(f'Unrecognised block "{block}" in json input file; '
                             'valid options are workflow/' + '/'.join(settings_classes.keys()))

    # Loading workflow settings
    workflow_settings = utils.parse_dict(bigdct.get('workflow', {}))
    task_name = workflow_settings.pop('task', 'singlepoint')

    # Load default values
    if 'setup' in bigdct:
        atoms, setup_parameters, psps_and_kpts, n_filled, n_empty = read_setup_dict(bigdct['setup'])
        del bigdct['setup']
    elif task_name != 'ui':
        raise ValueError('You must provide a "setup" block in the input file, specifying atomic positions, atomic '
                         'species, etc.')
    else:
        # Create dummy objects
        atoms = Atoms()
        setup_parameters = {}
        psps_and_kpts = {}
        n_filled = 0
        n_empty = 0

    # Loading calculator-specific settings
    master_calc_params = {}

    # Generate a master SettingsDict for every single kind of calculator, regardless of whether or not there was a
    # corresponding block in the json file
    for block, settings_class in settings_classes.items():
        # # For the UI task, the input file won't contain enough for us to generate all the other kinds of SettingsDicts
        # if task_name == 'ui' and block != 'ui':
        #     if block in bigdct and not block.startswith('ui'):
        #         utils.warn(f'Ignoring the {block} block since "task": "ui"')
        #     continue

        # Read the block and add the resulting calculator to the calcs_dct
        dct = bigdct.get(block, {})
        # Populating missing settings based on nelec, n_filled, n_empty etc
        if block == 'kcp':
            if 'nelec' not in dct.get('system', {}):
                dct['nelec'] = n_filled * 2
        elif block == 'pw':
            if 'nbnd' not in dct.get('system', {}):
                dct['nbnd'] = n_filled + n_empty
        if block == 'w90_occ':
            if 'num_wann' not in dct:
                dct['num_wann'] = n_filled
            if 'num_bands' not in dct:
                dct['num_bands'] = n_filled
        elif block == 'w90_emp':
            if 'num_wann' not in dct:
                dct['num_wann'] = n_empty
            if 'exclude_bands' not in dct:
                dct['exclude_bands'] = f'1-{n_filled}'
        elif block.startswith('ui'):
            # Dealing with redundancies in UI keywords
            if 'sc_dim' in dct and 'kpts' in psps_and_kpts:
                # In this case, the sc_dim keyword is redundant
                if psps_and_kpts['kpts'] != dct['sc_dim']:
                    raise ValueError('sc_dim in the UI block should match the kpoints provided in the setup block')
                dct.pop('sc_dim')
            if 'kpath' in dct and 'kpath' in psps_and_kpts:
                if psps_and_kpts['kpath'] != dct['kpath']:
                    raise ValueError('kpath in the UI block should match that provided in the setup block')
                dct.pop('kpath')

        master_calc_params[block] = settings_class(**dct)
        master_calc_params[block].update(
            **{k: v for k, v in setup_parameters.items() if k in master_calc_params[block].valid})
        master_calc_params[block].parse_algebraic_settings(nelec=n_filled * 2)

    name = fd.name.replace('.json', '')
    if task_name == 'singlepoint':
        workflow = workflows.SinglepointWorkflow(
            atoms, workflow_settings, master_calc_params, name=name, **psps_and_kpts)
    elif task_name == 'convergence':
        workflow = workflows.ConvergenceWorkflow(
            atoms, workflow_settings, master_calc_params, name=name, **psps_and_kpts)
    elif task_name in ['wannierize', 'wannierise']:
        workflow = workflows.WannierizeWorkflow(
            atoms, workflow_settings, master_calc_params, name=name, check_wannierisation=True, **psps_and_kpts)
    elif task_name == 'environ_dscf':
        workflow = workflows.DeltaSCFWorkflow(atoms, workflow_settings, master_calc_params, name=name, **psps_and_kpts)
    elif task_name == 'ui':
        workflow = workflows.UnfoldAndInterpolateWorkflow(
            atoms, workflow_settings, master_calc_params, name=name, **psps_and_kpts)
    else:
        raise ValueError('Invalid task name "{task_name}"')

    return workflow


def write_json(workflow, filename):
    '''

    Writes out settings to a JSON file

    '''

    fd = open(filename, 'w')

    calcs = workflow.master_calc_parameters

    bigdct = {}

    # "workflow" block (not printing any values that match the defaults)
    bigdct['workflow'] = {}
    for k, v in workflow.settings.items():
        if v is None:
            continue
        setting = [s for s in workflows.valid_settings if s.name == k][0]
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
        if isinstance(calc, (calculators.KoopmansCPCalculator, calculators.PWCalculator)):
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
