'''

utils module for python_KI

Written by Edward Linscott May 2020

'''

import warnings


def _warning(message, category=UserWarning, filename='', lineno=-1, file=None, line=None):
    '''
    Monkey-patching warnings.warn
    '''
    print(f'{category.__name__}: {message}')


warnings.showwarning = _warning


def warn(message):
    '''
    Allowing the monkey-patched warnings.warn to be imported as utils.warn
    '''
    warnings.warn(message)
