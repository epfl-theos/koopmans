from typing import Any, List

from koopmans import utils

from ._utils import Setting, SettingsDictWithChecks


class MLSettingsDict(SettingsDictWithChecks):
    def __init__(self, **kwargs) -> None:
        valid_settings = [
            Setting('use_ml', 'wheather to use a Machine Learning model to predict the alpha-parameters or not',
                    bool, False, (True, False)),
            Setting('n_max',
                    'The maximum expansion coefficient n for radial basis functions. If a list is provided in the convergence_ml-task, a grid search will be preformed.',
                    (int, list), 4, None),
            Setting('l_max',
                    'The maximum angular expansion coefficient. If a list is provided in the convergence_ml-task, a grid search will be preformed.',
                    (int, list), 4, None),
            Setting('r_min',
                    'The width of the narrowest radial basis function. If a list is provided in the convergence_ml-task, a grid search will be preformed.',
                    (float, list), 0.5, None),
            Setting('r_max',
                    'The width of the broadest radial basis function. If a list is provided in the convergence_ml-task, a grid search will be preformed.',
                    (float, list), 4.0, None),
            Setting('criterium',
                    'The criteium which has to be satisfied in order to use the ML-predicted screening coefficients instead of computing them ab-initio.',
                    str, 'after_fixed_num_of_snapshots', ('after_fixed_num_of_snapshots', )),
            Setting('number_of_training_snapshots',
                    'Number of snapshots needed for the "after_fixed_num_of_snapshots"-criterium. In case of the convergence_ml task, this number is taken to be the highest number of training samples for the convergence analysis.',
                    int, 1, None),
            Setting('current_snapshot',
                    'Number of snapshots already trained on',
                    int, 0, None),
            Setting('alphas_from_file_for_debugging_ml_model',
                    'If true, read the screening coefficients from file instead of calculating them ab-initio. The files have to be provided in the snapshot_ folders.',
                    bool, False, (True, False)),
            Setting('train_on_the_fly',
                    'If true, the ML-model gets trained after the calculation of each orbital. If false, the ML-model gets trained at the end of each snapshot.',
                    bool, False, (True, False)),
            Setting('occ_and_emp_together',
                    'If true, there will be one ML Model for both, occupied and empty states. If False, there will be one ML Model for occupied states and one for empty ones',
                    bool, True, (True, False)),
            Setting('type_of_ml_model',
                    'Which ML model to use for making the predictions',
                    str, 'ridge_regression', ('ridge_regression', 'linear_regression', 'mean')),
            Setting('input_data_for_ml_model',
                    'Which Data to use in case of the ridge_regression or the linear-regression Model',
                    str, 'orbital_density', ('orbital_density', 'self_hartree')),
            Setting('quantities_of_interest', 'Which quantities are we interested in the convergence_ml-task. Note that the eigenvalues (evs) require performing the final calculation afresh for every snapshot.',
                    (str, list), 'alphas', None)
        ]

        task = kwargs.pop('task')

        super().__init__(settings=valid_settings, **kwargs)

        if not task == 'convergence_ml':
            assert isinstance(self['n_max'], int)
            assert isinstance(self['l_max'], int)
            assert isinstance(self['r_min'], float)
            assert isinstance(self['r_max'], float)

        # convert now each parameter to a list to be able to run the same checks irregardless of the task

        def convert_to_list(param, type):
            if isinstance(param, type):  # if param is an int or a float convert it for the checks to a list
                return [param]
            else:  # if param is not an int or a float check that it is a list of ints / floats
                assert(isinstance(param, list))
                for value in param:
                    assert(isinstance(value, type))
                return param

        n_maxs = convert_to_list(self['n_max'], int)
        l_maxs = convert_to_list(self['l_max'], int)
        r_mins = convert_to_list(self['r_min'], float)
        r_maxs = convert_to_list(self['r_max'], float)

        # check that each n_max, l_max, r_max and r_min are greater or equal to 0 and that r_min is smaller than r_max
        for n_max in n_maxs:
            if not n_max > 0:
                raise ValueError(f"n_max has to be larger than zero. The provided value is n_max={n_max}")
        for l_max in l_maxs:
            if not l_max >= 0:
                raise ValueError(f"l_max has to be equal or larger than zero. The provided value is l_max={l_max}")
        for r_min in r_mins:
            if not r_min >= 0:
                raise ValueError(f"r_min has to be equal or larger than zero. The provided value is r_min={r_min}")
            if r_min < 0.5:
                utils.warn(
                    f"Small values of r_min (<0.5) can lead to problems in the construction of the radial basis. The provided value is r_min={r_min}.")
        for r_max in r_maxs:
            if not any(r_min < r_max for r_min in r_mins):
                raise ValueError(f"All provided values of r_min are larger or equal to r_max={r_max}.")

        # for the convergence_ml task we want to have each parameter in a list form
        if task == 'convergence_ml':

            implemented_quantities_of_interest = ['alphas', 'evs']
            self['quantities_of_interest'] = convert_to_list(self['quantities_of_interest'], str)

            for qoi in self['quantities_of_interest']:
                if qoi not in implemented_quantities_of_interest:
                    raise NotImplementedError(
                        "Performing the convergence_analysis w.r.t. {qoi} has not yet been implement.")

        def __setitem__(self, key: str, value: Any):

            return super().__setitem__(key, value)
