from typing import Any, List

from koopmans import utils

from ._utils import Setting, SettingsDictWithChecks


class MLSettingsDict(SettingsDictWithChecks):
    def __init__(self, **kwargs) -> None:
        valid_settings = [
            Setting('use_ml', 'use a machine learning model to predict the screening parameters',
                    bool, False, (True, False)),
            Setting('n_max',
                    'The maximum expansion coefficient n for radial basis functions. If a list is provided in the '
                    'convergence_ml-task, a grid search will be performed',
                    (int, list), 4, None),
            Setting('l_max',
                    'The maximum angular expansion coefficient. If a list is provided in the convergence_ml-task, a '
                    'grid search will be performed',
                    (int, list), 4, None),
            Setting('r_min',
                    'The width of the narrowest radial basis function. If a list is provided in the '
                    'convergence_ml-task, a grid search will be performed',
                    (float, list), 0.5, None),
            Setting('r_max',
                    'The width of the broadest radial basis function. If a list is provided in the '
                    'convergence_ml-task, a grid search will be performed',
                    (float, list), 4.0, None),
            Setting('criterium',
                    'The criterium which has to be satisfied in order to use the ML-predicted screening coefficients '
                    'instead of computing them ab-initio',
                    str, 'after_fixed_num_of_snapshots', ('after_fixed_num_of_snapshots', )),
            Setting('number_of_training_snapshots',
                    'Number of snapshots needed for the "after_fixed_num_of_snapshots"-criterium. In case of the '
                    'convergence_ml task, this number is taken to be the highest number of training samples for the '
                    'convergence analysis',
                    int, 1, None),
            Setting('current_snapshot',
                    'Number of snapshots already trained on',
                    int, 0, None),
            Setting('alphas_from_file',
                    'If true, read the screening coefficients from file instead of calculating them ab-initio. The '
                    'files have to be provided in the snapshot_ folders',
                    bool, False, (True, False)),
            Setting('train_on_the_fly',
                    'If true, the ML-model gets trained after the calculation of each orbital. If false, the ML-model '
                    'gets trained at the end of each snapshot',
                    bool, False, (True, False)),
            Setting('occ_and_emp_together',
                    'If true, there will be one ML model for both occupied and empty states. If False, there will be '
                    'one ML Model for occupied states and one for empty states',
                    bool, True, (True, False)),
            Setting('type_of_ml_model',
                    'Which ML model to use for making the predictions',
                    str, 'ridge_regression', ('ridge_regression', 'linear_regression', 'mean')),
            Setting('input_data_for_ml_model',
                    'Which data to use in case of the ridge_regression or the linear-regression Model',
                    str, 'orbital_density', ('orbital_density', 'self_hartree')),
            Setting('quantities_of_interest', 'Which quantities are we interested in the convergence_ml-task. Note '
                    'that the eigenvalues (evs) require performing the final calculation afresh for every snapshot.',
                    (str, list), ['alphas'], None)
        ]

        super().__init__(settings=valid_settings, **kwargs)
