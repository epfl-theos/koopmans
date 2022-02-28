from typing import Union, Dict, Any
from pathlib import Path
from ase.dft.kpoints import BandPath
from koopmans.io import write, read
from koopmans.workflows import SinglepointWorkflow


def qei_to_json(input_file: Union[str, Path], json: Union[str, Path], workflow_settings: Dict[str, Any] = {}):
    '''

    Converts a QE input file to a json file

    Arguments
    ---------
        input_file: the name of the QE input file to read in
        json: the name of the .json file to write out to
        workflow_settings: the koopmans.py keywords to add to the .json file (these are not
                    included in the QE input file)

    '''

    # Sanitise input
    input_file = Path(input_file)
    json = Path(json)

    if input_file.suffix == '.cpi':
        key = 'kcp'
    elif input_file.suffix == '.pwi':
        key = 'pw'
    else:
        raise ValueError('Unrecognised input file format')

    calc = read(input_file)
    calc.atoms.calc = None

    kwargs = {'pseudopotentials': calc.parameters.pop('pseudopotentials')}
    if key == 'pw':
        if isinstance(calc.parameters.kpts, BandPath):
            kwargs['kpath'] = calc.parameters.kpts
        else:
            kwargs['kgrid'] = calc.parameters.kpts

    wf = SinglepointWorkflow(atoms=calc.atoms,
                             parameters=workflow_settings,
                             master_calc_params={key: calc.parameters},
                             **kwargs)

    write(wf, json)

    return calc
