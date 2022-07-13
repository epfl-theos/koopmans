'''
Create a "stumbling" workflow that deliberately crashes the code after every single calculation and attempts to restart (for testing purposes)
'''

from __future__ import annotations

from contextlib import contextmanager
import copy
from pathlib import Path
from typing import Generator, Optional, Union

from ase.calculators.calculator import CalculationFailed
from koopmans import workflows


class DeliberateCalculationFailed(CalculationFailed):
    '''
    An error speciflcally for when we deliberately crash a calculation
    '''


stumble_message = 'Deliberately crashing for testing purposes'


class StumblingWorkflow:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc_counter = 1

    @property
    def calc_counter(self) -> int:
        return self._calc_counter

    @calc_counter.setter
    def calc_counter(self, value: int):
        self._calc_counter = value

    def run_calculator_single(self, *args, **kwargs):
        if len(self.calculations) == self.calc_counter:
            self.print(stumble_message)
            raise DeliberateCalculationFailed(stumble_message)
        else:
            super().run_calculator_single(*args, **kwargs)

    def _run(self, *args, **kwargs):
        if self.parent is None:
            # Create a copy of the state of the workflow before we start running it
            dct_before_running = {k: copy.deepcopy(v) for k, v in self.__dict__.items() if k not in
                                  ['calc_counter']}
            self.print(f'Attempt {self.calc_counter}', style='heading')
            try:
                super()._run(*args, **kwargs)
            except DeliberateCalculationFailed:
                # Restore the workflow to the state it was in before we ran any calculations and attempt to restart
                # from where we left off (keeping the calc_counter, which we excluded above)

                calc_counter = self.calc_counter
                self.__dict__ = dct_before_running

                # Update the parameters related to stumbling
                self.parameters.from_scratch = False
                self.calc_counter = calc_counter + 1

                # Attempt to run again
                self._run(*args, **kwargs)
        else:
            # Prevent subworkflows from catching stumbles
            self.calc_counter = self.parent.calc_counter
            super()._run(*args, **kwargs)
            self.parent.calc_counter = self.calc_counter


class StumblingSinglepointWorkflow(StumblingWorkflow, workflows.SinglepointWorkflow):
    pass


class StumblingConvergenceWorkflow(StumblingWorkflow, workflows.ConvergenceWorkflow):
    pass


class StumblingWannierizeWorkflow(StumblingWorkflow, workflows.WannierizeWorkflow):
    pass


class StumblingFoldToSupercellWorkflow(StumblingWorkflow, workflows.FoldToSupercellWorkflow):
    pass


class StumblingKoopmansDSCFWorkflow(StumblingWorkflow, workflows.KoopmansDSCFWorkflow):
    pass


class StumblingDFTCPWorkflow(StumblingWorkflow, workflows.DFTCPWorkflow):
    pass


class StumblingDFTPWWorkflow(StumblingWorkflow, workflows.DFTPWWorkflow):
    pass


class StumblingDeltaSCFWorkflow(StumblingWorkflow, workflows.DeltaSCFWorkflow):
    pass


class StumblingKoopmansDFPTWorkflow(StumblingWorkflow, workflows.KoopmansDFPTWorkflow):
    pass


class StumblingUnfoldAndInterpolateWorkflow(StumblingWorkflow, workflows.UnfoldAndInterpolateWorkflow):
    pass
