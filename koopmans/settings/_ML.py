import numpy as np
from typing import Any, List
from ._utils import Setting
from ._workflow import WorkflowSettingsDict
from koopmans import utils
from koopmans.ML_utils.ML_models import RidgeRegression

# TODO: perform more sanity checks, like: r_max>r_min, n_max>0, l_max>0, r_max<cell_size/2, the angles of the simulation cell are all 90° 

class MLSettingsDict(WorkflowSettingsDict):
	def __init__(self, **kwargs) -> None:
		valid_settings = [
			Setting('use_ML',
					'wheather to use the ML model or not',
					bool, False, (True, False)),
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
			Setting('number_of_snapshots',
					'Number of snapshots needed for the "after_fixed_num_of_snapshots"-criterium',
					int, 1, None),
			Setting('current_snapshot',
					'Number of snapshots already trained on',
					int, 0, None)]
		super().__init__(settings=valid_settings, **kwargs)


	def __setitem__(self, key: str, value: Any):
		return super().__setitem__(key, value)