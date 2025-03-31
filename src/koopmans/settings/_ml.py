from ._utils import Setting, SettingsDictWithChecks


class MLSettingsDict(SettingsDictWithChecks):
    """Settings for the machine learning configuration."""

    def __init__(self, **kwargs) -> None:
        valid_settings = [
            Setting('train', 'train a machine learning model to predict the screening parameters',
                    bool, False, (True, False)),
            Setting('test', 'test the machine learning model',
                    bool, False, (True, False)),
            Setting('predict', 'use a machine learning model to predict the screening parameters',
                    bool, False, (True, False)),
            Setting('model_file', 'JSON file containing the ML model information (generated from a '
                    'prior training calculation)',
                    str, None, None),
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
            Setting('alphas_from_file',
                    'If True, read the screening coefficients from file instead of calculating them ab initio',
                    bool, False, (True, False)),
            Setting('train_on_the_fly',
                    'If True, the ML-model gets trained after the calculation of each orbital. If False, the ML-model '
                    'gets trained at the end of each snapshot',
                    bool, False, (True, False)),
            Setting('occ_and_emp_together',
                    'If True, there will be one ML model for both occupied and empty states. If False, there will be '
                    'one ML Model for occupied states and one for empty states',
                    bool, True, (True, False)),
            Setting('estimator',
                    'What to use as the estimator for the ML model',
                    str, 'ridge_regression', ('ridge_regression', 'linear_regression', 'mean')),
            Setting('descriptor',
                    'What to use as the descriptor for the ML model',
                    str, 'orbital_density', ('orbital_density', 'self_hartree'))
        ]

        super().__init__(settings=valid_settings, **kwargs)
