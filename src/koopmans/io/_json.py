"""I/O of json input file for koopmans."""

import json as json_ext
from pathlib import Path
from typing import Any, Dict, TextIO, Type

from koopmans import utils, workflows
from koopmans.settings import WorkflowSettingsDict


def read_json(fd: TextIO, override: Dict[str, Any] = {}, **kwargs):
    """Read in settings listed in JSON file.

    Values in the JSON file can be overridden by values provided in the override argument
    """
    bigdct = json_ext.loads(fd.read())

    # Override all keywords provided explicitly
    utils.update_nested_dict(bigdct, override)

    # Loading workflow settings
    parameters = WorkflowSettingsDict(**utils.parse_dict(bigdct.get('workflow', {})))

    workflow_cls: Type[workflows.Workflow]
    if parameters.task == 'singlepoint':
        workflow_cls = workflows.SinglepointWorkflow
    elif parameters.task == 'wannierize':
        workflow_cls = workflows.WannierizeWorkflow
    elif parameters.task == 'dft_bands':
        workflow_cls = workflows.DFTBandsWorkflow
    elif parameters.task == 'dft_eps':
        workflow_cls = workflows.DFTPhWorkflow
    elif parameters.task == 'trajectory':
        workflow_cls = workflows.TrajectoryWorkflow
    else:
        raise ValueError(f'Invalid task name `{parameters.task}`')

    return workflow_cls.fromjson(fd.name, override, **kwargs)


def write_json(workflow: workflows.Workflow, filename: Path):
    """Write out settings to a JSON file."""
    fd = open(filename, 'w')

    json_ext.dump(workflow.toinputjson(), fd, indent=2)

    fd.close()

    return
