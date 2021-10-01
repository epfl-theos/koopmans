from glob import glob
from typing import Union, List, Type
from pathlib import Path
from koopmans.calculators import ExtendedCalculator, KoopmansCPCalculator, PWCalculator, Wannier90Calculator, \
    PW2WannierCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, KoopmansScreenCalculator, \
    KoopmansHamCalculator


def read_calculator(filenames: Union[Path, List[Path]]) -> ExtendedCalculator:

    # Interpreting the filenames argument
    if not isinstance(filenames, list):
        filenames = [filenames]
    valid_extensions = ['cpi', 'cpo', 'pwi', 'pwo', 'win', 'wout', 'p2wi',
                        'p2wo', 'uii', 'uio', 'w2ki', 'w2ko', 'ksi', 'kso', 'khi', 'kho']
    if not all([f.is_file() for f in filenames]):
        raise NotImplementedError('Need to code this addition of suffixes using Path')
        filenames = [f for prefix in filenames for f in glob(f'{prefix}.*') if f.split('.')[-1] in
                     valid_extensions]
    extensions = set([f.suffix for f in filenames])

    calc_class: Type[ExtendedCalculator]

    if extensions.issubset(set(['.cpi', '.cpo'])):
        calc_class = KoopmansCPCalculator
    elif extensions.issubset(set(['.pwi', '.pwo'])):
        calc_class = PWCalculator
    elif extensions.issubset(set(['.win', '.wout'])):
        calc_class = Wannier90Calculator
    elif extensions.issubset(set(['.p2wi', '.p2wo'])):
        calc_class = PW2WannierCalculator
    elif extensions.issubset(set(['.uii', '.uio'])):
        calc_class = UnfoldAndInterpolateCalculator
    elif extensions.issubset(set(['.w2ki', '.w2ko'])):
        calc_class = Wann2KCCalculator
    elif extensions.issubset(set(['.ksi', '.kso'])):
        calc_class = KoopmansScreenCalculator
    elif extensions.issubset(set(['.khi', '.kho'])):
        calc_class = KoopmansHamCalculator
    else:
        raise ValueError('Could not identify the extensions of ' + '/'.join([str(f) for f in filenames]))

    calc = calc_class.fromfile(filenames)

    return calc
