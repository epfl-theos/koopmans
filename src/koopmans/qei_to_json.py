"""Function for converting `Quantum ESPRESSO` input files to `koopmans` JSON input files."""

from pathlib import Path
from typing import Any, Dict, Union

from ase_koopmans.dft.kpoints import BandPath

from koopmans.engines import LocalhostEngine
from koopmans.io import read, write
from koopmans.kpoints import Kpoints
from koopmans.pseudopotentials import local_libraries
from koopmans.settings import WorkflowSettingsDict
from koopmans.utils.warnings import warn
from koopmans.workflows import SinglepointWorkflow


def qei_to_json(input_file: Union[str, Path], json: Union[str, Path],
                workflow_settings: Union[Dict[str, Any], WorkflowSettingsDict] = {}):
    """Convert a QE input file to a json file.

    Arguments
    ---------
        input_file: the name of the QE input file to read in
        json: the name of the .json file to write out to
        workflow_settings: the koopmans.py keywords to add to the .json file (these are not
                    included in the QE input file)

    """
    # Sanitizing input variables
    input_file = Path(input_file)
    json = Path(json)
    if not isinstance(workflow_settings, WorkflowSettingsDict):
        workflow_settings = WorkflowSettingsDict(**workflow_settings)

    # Checking we are starting from a sensible input file
    if input_file.suffix == '.cpi':
        key = 'kcp'
    elif input_file.suffix == '.pwi':
        key = 'pw'
    else:
        raise ValueError('Unrecognized input file format: must be either `.cpi` or `.pwi`')

    calc = read(input_file)
    calc.atoms.calc = None

    kwargs = {'pseudopotentials': calc.parameters.pop('pseudopotentials')}
    if key == 'pw':
        if isinstance(calc.parameters.kpts, BandPath):
            kwargs['kpoints'] = Kpoints(path=calc.parameters.kpts)
        else:
            kwargs['kpoints'] = Kpoints(grid=calc.parameters.kpts)

    # Determine the pseudopotential library
    pseudo_library = str(calc.parameters.pop('pseudo_dir'))
    default_library = 'SG15/1.2/PBE/SR'
    try:
        [matching_library] = [
            local_library for local_library in local_libraries if pseudo_library.endswith(local_library)]
    except ValueError:
        warn(f'Could not find a matching pseudopotential library; defaulting to {default_library}. If you want to use '
             'these specific pseudopotentials, install them via `koopmans pseudos install`')
        matching_library = default_library
        kwargs.pop('pseudopotentials', None)
    kwargs['pseudo_library'] = matching_library

    # Construct the workflow and write the input file to disk
    wf = SinglepointWorkflow(atoms=calc.atoms,
                             engine=LocalhostEngine(),
                             parameters=workflow_settings,
                             calculator_parameters={key: calc.parameters},
                             **kwargs)
    write(wf, json)

    return calc
