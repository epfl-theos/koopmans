import numpy as np
from typing import Any
from ._utils import SettingsDict
from koopmans.ML_utils.ML_models import RidgeRegression


class MLSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['use_ML', 'n_max', 'l_max', 'r_min', 'r_max', 'criterium', 'number_of_snapshots', 'current_snapshot'],
                         defaults={'use_ML': False, 'n_max': 4, 'l_max': 4,
                                   'r_min': 0.5, 'r_max': 4.0, 'criterium': 'after_fixed_num_of_snapshots', 'number_of_snapshots': 1, 'current_snapshot': 0},
                         **kwargs)


    def __setitem__(self, key: str, value: Any):
        return super().__setitem__(key, value)