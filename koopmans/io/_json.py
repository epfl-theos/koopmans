"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

import json as json_ext
from typing import TextIO, Dict, Any, Type
from pathlib import Path
from ase.io.espresso.koopmans_cp import construct_namelist
from koopmans import utils
from koopmans.settings import KoopmansCPSettingsDict, PWSettingsDict, WorkflowSettingsDict
from koopmans.utils import construct_cell_parameters_block
from koopmans import workflows


def update_nested_dict(dct_to_update, second_dct):
    for k, v in second_dct.items():
        if k in dct_to_update and isinstance(v, dict):
            update_nested_dict(dct_to_update[k], second_dct[k])
        else:
            dct_to_update[k] = v


def read_json(fd: TextIO, override: Dict[str, Any] = {}):
    '''

    Reads in settings listed in JSON file

    Values in the JSON file can be overridden by values provided in the override argument

    '''

    bigdct = json_ext.loads(fd.read())

    # Override all keywords provided explicitly
    update_nested_dict(bigdct, override)

    # Loading workflow settings
    parameters = WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))

    workflow_cls: Type[workflows.Workflow]
    if parameters.task == 'singlepoint':
        workflow_cls = workflows.SinglepointWorkflow
    elif parameters.task == 'convergence':
        workflow_cls = workflows.ConvergenceWorkflow
    elif parameters.task == 'wannierise':
        workflow_cls = workflows.WannierizeWorkflow
    elif parameters.task == 'environ_dscf':
        workflow_cls = workflows.DeltaSCFWorkflow
    elif parameters.task == 'ui':
        workflow_cls = workflows.SingleUnfoldAndInterpolateWorkflow
    elif parameters.task == 'dft_bands':
        workflow_cls = workflows.PWBandStructureWorkflow
    else:
        raise ValueError('Invalid task name "{task_name}"')

    return workflow_cls.fromjson(fd.name)


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
    if workflow.atoms.has('labels'):
        labels = workflow.atoms.get_array('labels')
    else:
        labels = workflow.atoms.get_chemical_symbols()
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
