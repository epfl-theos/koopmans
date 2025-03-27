"""Patches for testing the restart capabilities of `koopmans`."""

from __future__ import annotations

import copy

from ase.calculators.calculator import CalculationFailed

from koopmans import workflows


class DeliberateCalculationFailed(CalculationFailed):
    """An error speciflcally for when we deliberately crash a calculation."""


stumble_message = 'Deliberately crashing for testing purposes'


class StumblingWorkflow:
    """A workflow that deliberately crashes the code after every single calculation and attempts to restart."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.calc_counter = 1

    @property
    def calc_counter(self) -> int:
        """Count the number of calculations that have been run."""
        return self._calc_counter

    @calc_counter.setter
    def calc_counter(self, value: int):
        self._calc_counter = value

    def run_calculator_single(self, *args, **kwargs):
        """Run a single calculation, deliberately crashing at the first attempt."""
        if len(self.calculations) == self.calc_counter:
            self.print(stumble_message)
            raise DeliberateCalculationFailed(stumble_message)
        else:
            super().run_calculator_single(*args, **kwargs)

    def _run(self, *args, **kwargs):
        """Run the workflow."""
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


# TODO: monkeypatch only the methods, not the entire class; patch `Workflow` and not each subclass

class StumblingSinglepointWorkflow(StumblingWorkflow, workflows.SinglepointWorkflow):
    """A stumbling `SinglepointWorkflow`."""

    pass


class StumblingConvergenceWorkflow(StumblingWorkflow, workflows.ConvergenceWorkflow):
    """A stumbling `ConvergenceWorkflow`."""

    pass


class StumblingWannierizeWorkflow(StumblingWorkflow, workflows.WannierizeWorkflow):
    """A stumbling `WannierizeWorkflow`."""

    pass


class StumblingFoldToSupercellWorkflow(StumblingWorkflow, workflows.FoldToSupercellWorkflow):
    """A stumbling `FoldToSupercellWorkflow`."""

    pass


class StumblingKoopmansDSCFWorkflow(StumblingWorkflow, workflows.KoopmansDSCFWorkflow):
    """A stumbling `KoopmansDSCFWorkflow`."""

    pass


class StumblingDFTCPWorkflow(StumblingWorkflow, workflows.DFTCPWorkflow):
    """A stumbling `DFTCPWorkflow`."""

    pass


class StumblingDFTPhWorkflow(StumblingWorkflow, workflows.DFTPhWorkflow):
    """A stumbling `DFTPhWorkflow`."""

    pass


class StumblingDFTPWWorkflow(StumblingWorkflow, workflows.DFTPWWorkflow):
    """A stumbling `DFTPWWorkflow`."""

    pass


class StumblingDeltaSCFWorkflow(StumblingWorkflow, workflows.DeltaSCFWorkflow):
    """A stumbling `DeltaSCFWorkflow`."""

    pass


class StumblingKoopmansDFPTWorkflow(StumblingWorkflow, workflows.KoopmansDFPTWorkflow):
    """A stumbling `KoopmansDFPTWorkflow`."""

    pass


class StumblingUnfoldAndInterpolateWorkflow(StumblingWorkflow, workflows.UnfoldAndInterpolateWorkflow):
    """A stumbling `UnfoldAndInterpolateWorkflow`."""

    pass


class StumblingTrajectoryWorkflow(StumblingWorkflow, workflows.TrajectoryWorkflow):
    """A stumbling `TrajectoryWorkflow`."""

    pass
