from typing import Any, List
from ._utils import Setting
from ._workflow import WorkflowSettingsDict


class MLSettingsDict(WorkflowSettingsDict):
    def __init__(self, **kwargs) -> None:
        valid_settings = [
            Setting('n_max',
                    'The maximum expansion coefficient n for radial basis functions',
                    int, 4, None),
            Setting('l_max',
                    'The maximum angular expansion coefficient',
                    int, 4, None),
            Setting('r_min',
                    'The width of the narrowest radial basis function',
                    float, 0.5, None),
            Setting('r_max',
                    'The width of the broadest radial basis function',
                    float, 4.0, None),
            Setting('criterium',
                    'The criteium which has to be satisfied in order to use the ML-predicted screening coefficients instead of computing them ab-initio',
                    str, 'after_fixed_num_of_snapshots', ('after_fixed_num_of_snapshots', )),
            Setting('number_of_training_snapshots',
                    'Number of snapshots needed for the "after_fixed_num_of_snapshots"-criterium',
                    int, 1, None),
            Setting('current_snapshot',
                    'Number of snapshots already trained on',
                    int, 0, None),
            Setting('alphas_from_file_for_debugging_ml_model', 'If true, read the screening coefficients from file instead of calculating them ab-initio',
                    bool, False, (True, False)),
            Setting('train_on_the_fly', 'If true, the ML-model gets trained after the calculation of each orbital. If false, the ML-model gets trained at the end of each snapshot',
                    bool, False, (True, False)),
            Setting('occ_and_emp_together', 'If true, there will be one ML Model for both, occupied and empty states. If False, there will be one ML Model for occupied states and one for empty ones',
                    bool, True, (True, False)),
            Setting('type_ml_model', 'Which ML model to use for making the predictions',
                    str, 'Ridge Regression', ('Ridge Regression', 'Mean')),
            Setting('input_data', 'Which Data to use in case of the Ridge Regression Model',
                    str, 'Orbital Density', ('Orbital Density', 'Self-Hartree'))]

        super().__init__(settings=valid_settings, **kwargs)

        if not self['r_min'] < self['r_max']:
            raise ValueError(
                f"r_min should be smaller smaller than r_max. The provided values are r_min={self['r_min']} and r_max={self['r_max']}")
        if not self['l_max'] >= 0:
            raise ValueError(
                f"l_max has to be equal or larger than zero. The provided value is l_max={self['l_max']}")
        if not self['n_max'] > 0:
            raise ValueError(f"n_max has to be larger than zero. The provided value is n_max={self['n_max']}")

    def __setitem__(self, key: str, value: Any):
        return super().__setitem__(key, value)
