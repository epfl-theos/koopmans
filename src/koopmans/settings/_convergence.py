from typing import Any

from ._utils import Setting, SettingsDictWithChecks


class ConvergenceSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs) -> None:
        settings = [
            Setting('observable',
                    'System observable of interest which we converge',
                    str, 'total energy', None),
            Setting('threshold',
                    'threshold for the convergence of the observable of interest',
                    (str, float), None, None),
            Setting('variables',
                    'The observable of interest will be converged with respect to this (these) '
                    'simulation variable(s)',
                    (list, str), ['ecutwfc'], None),
            Setting('max_steps',
                    'The maximum number of values of each variable to try when attempting to achieve convergence',
                    int, 10, None)]

        super().__init__(settings=settings, physicals=['threshold'], **kwargs)

    @property
    def _other_valid_keywords(self):
        return []

    def __setitem__(self, key: str, value: Any):
        # Convert parameters to a list
        if key == 'parameters' and isinstance(value, str):
            value = [value]

        return super().__setitem__(key, value)
