'''

warnings functions for koopmans.utils

Written by Edward Linscott May 2020

'''


import sys
import traceback
import warnings
from typing import Optional, TextIO, Type, Union


class CalculatorNotConvergedWarning(UserWarning):
    pass


def _warning(message: Union[str, Warning], category: Type[Warning] = UserWarning, filename: str = '',
             lineno: int = -1, file: Optional[TextIO] = None, line: Optional[str] = None) -> None:
    '''
    Monkey-patching warnings.warn
    '''
    print(f'{category.__name__}: {message}')


def _warn_with_traceback(message: Union[str, Warning], category: Type[Warning] = UserWarning, filename: str = '',
                         lineno: int = -1, file: Optional[TextIO] = None, line: Optional[str] = None) -> None:
    log = file if file and hasattr(file, 'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


# warnings.showwarning = _warn_with_traceback
warnings.showwarning = _warning


def warn(message: str, cls: Type[Warning] = UserWarning) -> None:
    '''
    Allowing the monkey-patched warnings.warn to be imported as utils.warn
    '''
    warnings.warn(message, cls)
