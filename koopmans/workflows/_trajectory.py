

import copy
import json as json_ext
import os
import statistics
from pathlib import Path
from tracemalloc import Snapshot
from typing import Any, Dict, List, Optional, Tuple, Union
from unittest import result

import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
from sklearn.metrics import mean_absolute_error, r2_score

from ase import Atoms, io
from ase.calculators.calculator import Calculator
from koopmans import calculators, utils
from koopmans.ml_utils import RidgeRegression

from ._workflow import Workflow

load_results_from_output = True


class TrajectoryWorkflow(Workflow):

    def __init__(self, snapshots, *args, **kwargs):
        self.snapshots: List[Atoms] = snapshots
        self.number_of_snapshots = len(self.snapshots)
        super().__init__(*args, **kwargs)

    def _run(self, indices: Optional[List[int]] = None, save_dir: Optional[Path] = None, delete_final_dir: Optional[bool] = False):
        """
        Runs the trajectory workflow.
        """

        # Import it like this so if they have been monkey-patched, we will get the monkey-patched version
        from koopmans.workflows import (KoopmansDFPTWorkflow,
                                        KoopmansDSCFWorkflow)
        if indices is None:
            indices = range(0, self.number_of_snapshots)

        for i in indices:
            subdirectory = f'snapshot_{i+1}'
            snapshot = self.snapshots[i]
            self.parameters.current_snapshot = i
            self.atoms.set_positions(snapshot.positions)

            self.print(
                f'Performing Koopmans calculation on snapshot {i+1} / {self.number_of_snapshots}', style='heading')

            workflow: Union[KoopmansDFPTWorkflow, KoopmansDSCFWorkflow]
            if self.parameters.method == 'dfpt':
                workflow = KoopmansDFPTWorkflow(**self.wf_kwargs)
                raise NotImplementedError("TODO Yannick: What is the final_calculator of dftp")
            elif self.parameters.method == 'dscf':
                workflow = KoopmansDSCFWorkflow(**self.wf_kwargs)
                final_calculator = calculators.KoopmansCPCalculator
            else:
                raise NotImplementedError("The trajectory workflow is currently only implemented with dfpt and dscf.")
            # reset the bands to the initial guesses (i.e. either from file or to 0.6 but not from the previous calculation)
            self.bands = workflow.bands
            if delete_final_dir:
                utils.system_call(f'rm -r {subdirectory}/final')

            self.run_subworkflow(workflow, subdirectory=subdirectory)

            # since we have deleted the final directory and therefore had to rerun, we must now make sure, that from_scratch is again set to false
            if delete_final_dir:
                self.parameters.from_scratch = False

            if save_dir is not None:
                final_calc = [c for c in workflow.calculations if isinstance(c, final_calculator)][-1]
                evs = final_calc.results['eigenvalues']
                alphas = self.bands.alphas
                np.savetxt(save_dir / f"evs_snapshot_{i+1}.txt", evs)
                np.savetxt(save_dir / f"alphas_snapshot_{i+1}.txt", alphas)

    @ classmethod
    def _fromjsondct(cls, bigdct: Dict[str, Any]):
        """
        Reads the atomic positions for each snapshot from the xyz-file specified by the user in the snapshots-file.
        """

        try:
            snapshots_file = bigdct['setup'].pop('snapshots')
        except:
            raise ValueError(
                f'To calculate a trajectory, please provide a xyz-file containing the atomic positions of the snapshots in the setup-block of the json-input file.')

        snapshots = io.read(snapshots_file, index=':')
        if isinstance(snapshots, Atoms):
            snapshots = [snapshots]
        bigdct['setup']['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])
        wf = super(TrajectoryWorkflow, cls)._fromjsondct(bigdct)
        wf.snapshots = snapshots
        wf.number_of_snapshots = len(snapshots)
        return wf

    # @classmethod
    # def fromdict(cls, dct: Dict):

    #     if dct['parameters']['use_ml'] and

    #     wf = cls(snapshots=dct.pop('snapshots'),
    #              parameters=dct.pop('parameters'),
    #              master_calc_params=dct.pop('master_calc_params'),
    #              pseudopotentials=dct.pop('_pseudopotentials'),
    #              gamma_only=dct.pop('_gamma_only'),
    #              kgrid=dct.pop('_kgrid'),
    #              kpath=dct.pop('_kpath'),
    #              projections=dct.pop('projections'),
    #              autogenerate_settings=False)

    #     for k, v in dct.items():
    #         setattr(wf, k, v)

    #     return wf

    #     try:
    #         snapshots = bigdct.pop('snapshots')
    #     except:
    #         raise ValueError(
    #             f'To calculate a trajectory, please provide a xyz-file containing the atomic positions of the snapshots in the setup-block of the json-input file.')

    #     # snapshots = io.read(snapshots_file, index=':')
    #     if isinstance(snapshots, Atoms):
    #         snapshots = [snapshots]
    #     bigdct['atomic_positions'] = utils.construct_atomic_positions_block(snapshots[0])
    #     import ipdb
    #     ipdb.set_trace()
    #     wf = super(TrajectoryWorkflow, cls).fromdict(bigdct)
    #     wf.snapshots = snapshots
    #     return wf

    # def todict(self):

    #     dct = dict(self.__dict__)

    #     items_to_pop = ['atoms']
    #     for item in items_to_pop:
    #         dct.pop(item)

    #     if self.parameters.use_ml and self.parameters.mode == 'convergence':
    #         items_to_pop = ['metrics', 'statistics']
    #         for item in items_to_pop:
    #             dct.pop(item)
    #     return dct
