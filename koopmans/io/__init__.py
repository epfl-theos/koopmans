"""

I/O module for koopmans

Written by Edward Linscott Jan 2020

"""

import os
from glob import glob
from typing import TextIO, Union, List, Type
from ._json import read_json, write_json, read_ui_dict
from ._kwf import read_kwf, write_kwf
from ._utils import indented_print, write_alpha_file, read_alpha_file, construct_cell_parameters_block, \
    nelec_from_pseudos, read_kpath, read_cell_parameters
from koopmans.calculators.generic import GenericCalc
from koopmans.workflows.generic import Workflow


def read(filename: str, **kwargs) -> Workflow:

    # Generic "read" function

    if filename.endswith('kwf'):
        with open(filename, 'w') as fd:
            out = read_kwf(fd)
        return out
    elif filename.endswith('.json'):
        return read_json(filename, **kwargs)
    else:
        raise ValueError(f'Unrecognised file type for {filename}')


def write(obj: Workflow, filename: str):

    # Generic "write" function

    if filename.endswith('kwf'):
        with open(filename, 'w') as fd:
            write_kwf(obj, fd)
    elif filename.endswith('.json'):
        write_json(obj, filename)
    else:
        raise ValueError(f'Unrecognised file type for {filename}')


def load_calculator(filenames: Union[str, List[str]]) -> GenericCalc:

    from koopmans.calculators import kcp, pw, wannier90, pw2wannier, ui, wann2kc, kc_screen, kc_ham

    if not isinstance(filenames, list):
        filenames = [filenames]

    valid_extensions = ['cpi', 'cpo', 'pwi', 'pwo', 'win', 'wout', 'p2wi',
                        'p2wo', 'uii', 'uio', 'w2ki', 'w2ko', 'ksi', 'kso', 'khi', 'kho']
    if not all([os.path.isfile(f) for f in filenames]):
        filenames = [f for prefix in filenames for f in glob(f'{prefix}.*') if f.split('.')[-1] in
                     valid_extensions]

    extensions = set([f.split('.')[-1] for f in filenames])

    calc_class: Type[GenericCalc]

    if extensions.issubset(set(['cpi', 'cpo'])):
        calc_class = kcp.KCP_calc
    elif extensions.issubset(set(['pwi', 'pwo'])):
        calc_class = pw.PW_calc
    elif extensions.issubset(set(['win', 'wout'])):
        calc_class = wannier90.W90_calc
    elif extensions.issubset(set(['p2wi', 'p2wo'])):
        calc_class = pw2wannier.PW2Wannier_calc
    elif extensions.issubset(set(['uii', 'uio'])):
        calc_class = ui.UI_calc
    elif extensions.issubset(set(['w2ki', 'w2ko'])):
        calc_class = wann2kc.Wann2KCCalc
    elif extensions.issubset(set(['ksi', 'kso'])):
        calc_class = kc_screen.KoopmansScreenCalc
    elif extensions.issubset(set(['khi', 'kho'])):
        calc_class = kc_ham.KoopmansHamCalc
    else:
        raise ValueError('Could not identify the extensions of ' + '/'.join(filenames))

    return calc_class(qe_files=filenames)
