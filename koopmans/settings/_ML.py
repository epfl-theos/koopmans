import numpy as np
from typing import Any
from ._utils import SettingsDict
from koopmans.ML_utils.ML_models import RidgeRegression


class MLSettingsDict(SettingsDict):
    def __init__(self, **kwargs) -> None:

        super().__init__(valid=['use_ML', 'n_max', 'l_max', 'r_min', 'r_max'],
                         defaults={'use_ML': False, 'n_max': 4, 'l_max': 4,
                                   'r_min': 0.5, 'r_max': 4.0},
                         **kwargs)

    def update(self, *args, **kwargs) -> None:
        super().update(*args, **kwargs)

    # @property
    # def _other_valid_keywords(self):
    #     return ['kgrid', 'kpath']

    def __setitem__(self, key: str, value: Any):
        if key == 'use_ML':
            assert isinstance(value, bool)
            if value == True:
                self.ml_model = RidgeRegression()

        return super().__setitem__(key, value)