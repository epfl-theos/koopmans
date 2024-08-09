"""

kwf (Koopmans WorkFlow) I/O for koopmans

Written by Edward Linscott Mar 2021, largely modelled off ase.io.jsonio

"""

from typing import TextIO, Union

import koopmans.workflows as workflows
from koopmans.utils import serialization


def read_kwf(fd: TextIO):
    return serialization.decode(fd.read())


def write_kwf(obj: Union[workflows.Workflow, dict], fd: TextIO):
    if isinstance(obj, workflows.Workflow):
        use_relpath = obj.parameters.use_relative_paths
        obj.parameters.use_relative_paths = True
    fd.write(serialization.encode(obj))
    if isinstance(obj, workflows.Workflow):
        obj.parameters.use_relative_paths = use_relpath
