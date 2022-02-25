"""

I/O of json input file for koopmans

Written by Edward Linscott Jan 2020

"""

import json as json_ext
from typing import TextIO, Dict, Any, Type
from pathlib import Path
from koopmans import utils
from koopmans.settings import WorkflowSettingsDict
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

    json_ext.dump(workflow.toinputjson(), fd, indent=2)

    fd.close()

    return
