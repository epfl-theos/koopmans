"""
Settings module for plotting

"""

from typing import List

from ._utils import Setting, SettingsDictWithChecks

valid_settings: List[Setting] = [
    Setting('degauss',
            'gaussian broadening (in eV) for the DOS interpolation, as in QE',
            (str, float, int), 0.05, None),
    Setting('nstep',
            'number of steps for the plot of the interpolated DOS',
            int, 1000, None),
    Setting('Emin',
            'minimum energy for the plot of the interpolated DOS',
            (str, float, int), None, None),
    Setting('Emax',
            'maximum energy for the plot of the interpolated DOS',
            (str, float, int), None, None)]


class PlotSettingsDict(SettingsDictWithChecks):

    def __init__(self, **kwargs):
        super().__init__(settings=valid_settings, **kwargs)
