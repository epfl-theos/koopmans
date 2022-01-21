"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

import os
import json as json_ext
import xml.etree.ElementTree as ET
from typing import TextIO, Dict, List, Union
from pathlib import Path
from ase.atoms import Atoms
from ase.io.espresso.utils import ibrav_to_cell
from ase.calculators.espresso import Espresso_kcp
from ase.io.espresso.koopmans_cp import KEYS as kcp_keys, construct_namelist
from koopmans import utils, bands, projections
from koopmans.pseudopotentials import set_up_pseudos, nelec_from_pseudos, valence_from_pseudo
from koopmans.settings import KoopmansCPSettingsDict, KoopmansHamSettingsDict, KoopmansScreenSettingsDict, \
    PWSettingsDict, PW2WannierSettingsDict, UnfoldAndInterpolateSettingsDict, Wann2KCSettingsDict, \
    Wannier90SettingsDict, WorkflowSettingsDict
from koopmans.utils import read_atomic_species, read_atomic_positions, read_cell_parameters, read_kpoints_block, \
    construct_cell_parameters_block
from koopmans import workflows


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

        # If they are missing, fill the nelec/nelup/neldw fields using information contained
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
            if 'starting_magnetization(1)' in calc.parameters:
                labels = [s + str(t) if t != 0 else s for s, t in zip(calc.atoms.symbols, calc.atoms.get_tags())]
                starting_magmoms = {}
                for i, (l, p) in enumerate(calc.parameters.pseudopotentials.items()):
                    # ASE uses absoulte values; QE uses the fraction of the valence
                    frac_mag = calc.parameters.pop(f'starting_magnetization({i + 1})', 0.0)
                    valence = valence_from_pseudo(p, calc.parameters.pseudo_dir)
                    starting_magmoms[l] = frac_mag * valence
                calc.atoms.set_initial_magnetic_moments([starting_magmoms[l] for l in labels])
            elif tot_mag != 0:
                calc.atoms.set_initial_magnetic_moments([tot_mag / len(calc.atoms) for _ in calc.atoms])

        # Work out the number of filled and empty bands
        n_filled = nelec // 2 + nelec % 2
        if 'nbnd' in calc.parameters:
            n_empty = calc.parameters.nbnd - n_filled
        else:
            n_empty = 0
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


def update_nested_dict(dct_to_update, second_dct):
    for k, v in second_dct.items():
        if k in dct_to_update and isinstance(v, dict):
            update_nested_dict(dct_to_update[k], second_dct[k])
        else:
            dct_to_update[k] = v


def read_json(fd: TextIO, override={}):
    '''

    Reads in settings listed in JSON file

    Values in the JSON file can be overridden by values provided in the override argument

    '''

    bigdct = json_ext.loads(fd.read())

    # Override all keywords provided explicitly
    update_nested_dict(bigdct, override)

    # Deal with the nested w90 subdictionaries
    if 'w90' in bigdct:
        for filling in ['occ', 'emp']:
            for spin in ['up', 'down']:
                # Add any keywords in the filling:spin subsubdictionary
                subsubdct = bigdct['w90'].get(filling, {}).get(spin, {})
                bigdct[f'w90_{filling}_{spin}'] = subsubdct
                # Add any keywords in the filling subdictionary
                subdct = {k: v for k, v in bigdct['w90'].get(filling, {}).items() if k not in ['up', 'down']}
                bigdct[f'w90_{filling}_{spin}'].update(subdct)
                # Add any keywords in the main dictionary
                dct = {k: v for k, v in bigdct['w90'].items() if k not in ['occ', 'emp']}
                bigdct[f'w90_{filling}_{spin}'].update(dct)
            # Also create a spin-independent set of parameters
            bigdct[f'w90_{filling}'] = {}
            bigdct[f'w90_{filling}'].update(subdct)
            bigdct[f'w90_{filling}'].update(dct)
        # Finally, remove the nested w90 entry
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
                        'w90_emp': Wannier90SettingsDict,
                        'w90_occ_up': Wannier90SettingsDict,
                        'w90_emp_up': Wannier90SettingsDict,
                        'w90_occ_down': Wannier90SettingsDict,
                        'w90_emp_down': Wannier90SettingsDict}

    # Check for unexpected blocks
    for block in bigdct:
        if block not in list(settings_classes.keys()) + ['workflow', 'setup']:
            raise ValueError(f'Unrecognised block "{block}" in json input file; '
                             'valid options are workflow/' + '/'.join(settings_classes.keys()))

    # Loading workflow settings
    parameters = WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))

    # Load default values
    if 'setup' in bigdct:
        atoms, setup_parameters, workflow_kwargs, n_filled, n_empty = read_setup_dict(bigdct['setup'], parameters.task)
        del bigdct['setup']
    elif parameters.task != 'ui':
        raise ValueError('You must provide a "setup" block in the input file, specifying atomic positions, atomic '
                         'species, etc.')
    else:
        # Create dummy objects
        atoms = Atoms()
        setup_parameters = {}
        workflow_kwargs = {}
        n_filled = 0
        n_empty = 0

    # Loading calculator-specific settings
    master_calc_params = {}
    w90_block_projs: List = []
    w90_block_filling: List[bool] = []
    w90_block_spins: List[Union[str, None]] = []

    # Generate a master SettingsDict for every single kind of calculator, regardless of whether or not there was a
    # corresponding block in the json file
    for block, settings_class in settings_classes.items():
        # Read the block and add the resulting calculator to the calcs_dct
        dct = bigdct.get(block, {})
        # Populating missing settings based on nelec, n_filled, n_empty etc
        if block == 'kcp':
            if 'nelec' not in dct.get('system', {}):
                dct['nelec'] = setup_parameters['nelec']
        elif block == 'pw':
            if 'nbnd' not in dct.get('system', {}):
                dct['nbnd'] = n_filled + n_empty
        elif block.startswith('ui'):
            # Dealing with redundancies in UI keywords
            if 'sc_dim' in dct and 'kpts' in workflow_kwargs:
                # In this case, the sc_dim keyword is redundant
                if workflow_kwargs['kpts'] != dct['sc_dim']:
                    raise ValueError('sc_dim in the UI block should match the kpoints provided in the setup block')
                dct.pop('sc_dim')
            if 'kpath' in dct and 'kpath' in workflow_kwargs:
                if workflow_kwargs['kpath'] != dct['kpath']:
                    raise ValueError('kpath in the UI block should match that provided in the setup block')
                dct.pop('kpath')
        elif block.startswith('w90'):
            # If we are spin-polarised, don't store the spin-independent w90 block
            # Likewise, if we are not spin-polarised, don't store the spin-dependent w90 blocks
            if parameters.spin_polarised is not ('up' in block or 'down' in block):
                continue
            if 'projections' in dct and 'projections_blocks' in dct:
                raise ValueError('"projections" and "projections_block" are mutually exclusive')
            elif 'projections_blocks' in dct:
                projs = dct.pop('projections_blocks')
            else:
                projs = [dct.pop('projections', [])]
            w90_block_projs += projs
            w90_block_filling += ['occ' in block for _ in range(len(projs))]
            if 'up' in block:
                w90_block_spins += ['up' for _ in range(len(projs))]
            elif 'down' in block:
                w90_block_spins += ['down' for _ in range(len(projs))]
            else:
                w90_block_spins += [None for _ in range(len(projs))]
            for kw in ['exclude_bands', 'num_wann', 'num_bands', 'projections']:
                if kw in dct:
                    utils.warn(f'{kw} will be overwritten by the workflow; it is best to leave this keyword out of the '
                               'JSON input file and to then double-check this keyword in the various .win files '
                               'generated by the workflow.')

        master_calc_params[block] = settings_class(**dct)
        master_calc_params[block].update(
            **{k: v for k, v in setup_parameters.items() if k in master_calc_params[block].valid})
        master_calc_params[block].parse_algebraic_settings(nelec=setup_parameters['nelec'])

    # Adding the projections to the workflow kwargs (this is unusual in that this is an attribute of the workflow object
    # but it is provided in the w90 subdictionary)
    workflow_kwargs['projections'] = projections.ProjectionBlocks.fromprojections(
        w90_block_projs, w90_block_filling, w90_block_spins, atoms)

    name = fd.name.replace('.json', '')
    workflow: workflows.Workflow
    if parameters.task == 'singlepoint':
        workflow = workflows.SinglepointWorkflow(
            atoms, parameters, master_calc_params, name=name, **workflow_kwargs)
    elif parameters.task == 'convergence':
        workflow = workflows.ConvergenceWorkflow(
            atoms, parameters, master_calc_params, name=name, **workflow_kwargs)
    elif parameters.task == 'wannierise':
        workflow = workflows.WannierizeWorkflow(
            atoms, parameters, master_calc_params, name=name, **workflow_kwargs)
        workflow.parameters.check_wannierisation = True
    elif parameters.task == 'environ_dscf':
        workflow = workflows.DeltaSCFWorkflow(atoms, parameters, master_calc_params,
                                              name=name, **workflow_kwargs)
    elif parameters.task == 'ui':
        workflow = workflows.UnfoldAndInterpolateWorkflow(
            atoms, parameters, master_calc_params, name=name, **workflow_kwargs)
    elif parameters.task == 'dft_bands':
        workflow = workflows.PWBandStructureWorkflow(
            atoms, parameters, master_calc_params, name=name, **workflow_kwargs)
    else:
        raise ValueError('Invalid task name "{task_name}"')

    return workflow


def write_json(workflow: workflows.Workflow, filename: Path):
    '''

    Writes out settings to a JSON file

    '''

    fd = open(filename, 'w')

    bigdct: Dict[str, Dict] = {}

    # "workflow" block (not printing any values that match the defaults)
    bigdct['workflow'] = {}
    for k, v in workflow.parameters.items():
        if v is None:
            continue
        if isinstance(v, Path):
            v = str(v)
        default = workflow.parameters.defaults.get(k, None)
        if v != default:
            bigdct['workflow'][k] = v

    # "setup" block
    # Working out ibrav
    ibrav = workflow.master_calc_params['kcp'].get('ibrav', workflow.master_calc_params['pw'].get('ibrav', 0))

    bigdct['setup'] = {}

    # cell parameters
    if ibrav == 0:
        bigdct['setup']['cell_parameters'] = construct_cell_parameters_block(workflow.atoms)

    # atomic positions
    if len(set(workflow.atoms.get_tags())) > 1:
        labels = [s + str(t) if t > 0 else s for s, t in zip(workflow.atoms.symbols, workflow.atoms.get_tags())]
    else:
        labels = workflow.atoms.symbols

    if ibrav == 0:
        bigdct['setup']['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, workflow.atoms.get_positions())],
            'units': 'angstrom'}
    else:
        bigdct['setup']['atomic_positions'] = {'positions': [
            [label] + [str(x) for x in pos] for label, pos in zip(labels, workflow.atoms.get_scaled_positions())],
            'units': 'crystal'}

    # atomic species
    bigdct['setup']['atomic_species'] = {'species': [[key, 1.0, val]
                                                     for key, val in workflow.pseudopotentials.items()]}

    for code, params in workflow.master_calc_params.items():
        if isinstance(params, (KoopmansCPSettingsDict, PWSettingsDict)):
            bigdct[code] = {}

            # pseudo directory
            pseudo_dir = params.pop('pseudo_dir', None)
            if pseudo_dir is not None:
                bigdct['setup']['control'] = {'pseudo_dir': str(pseudo_dir)}

            # Remove default settings and convert Paths to strings
            params_dict = {k: v for k, v in params.items() if params.defaults.get(k, None) != v}
            for k in params_dict:
                if isinstance(params_dict[k], Path):
                    params_dict[k] = str(params_dict[k])

            # Populate bigdct with the settings
            input_data = construct_namelist(params_dict)
            for key, block in input_data.items():

                if len(block) > 0:
                    bigdct[code][key] = {k: v for k, v in dict(
                        block).items() if v is not None}

            if code == 'pw':
                # Adding kpoints to "setup"
                kpts = {'kgrid': workflow.kgrid,
                        'kpath': workflow.kpath.path}
                bigdct['setup']['k_points'] = kpts
        else:
            raise NotImplementedError(
                f'Writing of {params.__class__.__name__} with write_json is not yet implemented')

    json_ext.dump(bigdct, fd, indent=2)

    fd.close()

    return
