"""I/O of json input file for koopmans."""

import json as json_ext
from pathlib import Path
from typing import Any, Dict, TextIO, Type

from koopmans import utils, workflows
from koopmans.workflows.config import WorkflowConfig, Task


def read_json(fd: TextIO, override: Dict[str, Any] = {}, **kwargs):
    """Read in settings listed in JSON file.

    Values in the JSON file can be overridden by values provided in the override argument
    """
    bigdct = json_ext.loads(fd.read())

    # Override all keywords provided explicitly
    utils.update_nested_dict(bigdct, override)

    # Loading workflow settings
    parameters = WorkflowConfig(**bigdct.get('workflow', {}))

    workflow_cls: Type[workflows.Workflow]
    if parameters.task == Task.SINGLEPOINT:
        workflow_cls = workflows.SinglepointWorkflow
    elif parameters.task == Task.WANNIERIZE:
        workflow_cls = workflows.WannierizeWorkflow
    elif parameters.task == Task.DFT_BANDS:
        workflow_cls = workflows.DFTBandsWorkflow
    elif parameters.task == Task.DFT_EPS:
        workflow_cls = workflows.DFTPhWorkflow
    elif parameters.task == Task.TRAJECTORY:
        workflow_cls = workflows.TrajectoryWorkflow
    else:
        raise NotImplementedError(f'Need to implement `{parameters.task.value}`')

    return workflow_cls.fromjson(fd.name, override, **kwargs)


def write_json(workflow: workflows.Workflow, filename: Path):
    """Write out settings to a JSON file."""
    fd = open(filename, 'w')

    json_ext.dump(workflow.toinputjson(), fd, indent=2)

    fd.close()

    return
