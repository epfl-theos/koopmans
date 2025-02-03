import sys
from abc import ABC, abstractmethod
from typing import Generator, Optional

from upf_tools import UPFDict

from koopmans import utils
from koopmans.files import FilePointer
from koopmans.status import Status
from koopmans.step import Step


class Engine(ABC):

    def __init__(self, from_scratch: bool = True, npool: Optional[int] = None):
        self.from_scratch = from_scratch
        self.npool = npool

    def print(self, text: str, **kwargs):
        utils.indented_print(text, **kwargs, wrap=False)

    def _step_message(self, uid, symbol, suffix, **kwargs):
        self.print(f'- {symbol} `{uid}` {suffix}', **kwargs)

    def _step_completed_message(self, step):
        self._step_completed_message_by_uid(step.uid)

    def _step_failed_message(self, step):
        self._step_failed_message_by_uid(step.uid)

    def _step_running_message(self, step):
        self._step_running_message_by_uid(step.uid)

    def _step_skipped_message(self, step):
        self._step_skipped_message_by_uid(step.uid)

    def _step_completed_message_by_uid(self, uid):
        self._step_message(uid, '✅', 'completed  ')

    def _step_failed_message_by_uid(self, uid):
        self._step_message(uid, '❌', 'failed     ')

    def _step_running_message_by_uid(self, uid):
        if sys.stdout.isatty():
            self._step_message(uid, '🖥️ ', 'running...', end='\r', flush=True)

    def _step_skipped_message_by_uid(self, uid):
        self._step_message(uid, '⏭️ ', 'already complete  ')

    @abstractmethod
    def run(self, step: Step) -> None:
        ...

    @abstractmethod
    def get_status(self, step: Step) -> Status:
        ...

    @abstractmethod
    def set_status(self, step: Step, status: Status):
        ...

    @abstractmethod
    def update_statuses(self) -> None:
        ...

    @abstractmethod
    def load_results(self, step: Step) -> None:
        ...

    @abstractmethod
    def get_pseudopotential(self, library: str, element: str) -> UPFDict:
        ...

    @abstractmethod
    def read(self, file: FilePointer, binary: bool = False) -> str | bytes:
        ...

    @abstractmethod
    def write(self, content: str | bytes, file: FilePointer) -> None:
        ...

    @abstractmethod
    def glob(self, directory: FilePointer, pattern: str, recursive: bool = False) -> Generator[FilePointer, None, None]:
        ...

    @abstractmethod
    def available_pseudo_families(self) -> set[str]:
        ...

#
#        # Print the header
#        self.print(header())
#
#        # Print the bibliography
#        workflow.print_bib()
#
#        steps_generator = workflow.steps_generator()
#        while True:
#            steps_to_run: tuple[Step, ...]
#            try:
#                steps_to_run = next(steps_generator)
#            except StopIteration:
#                break
#
#            self.run_steps(steps_to_run)
#
#        # Print farewell message
#        self.print('\n**Workflow complete** 🎉')
#
#    # @abstractmethod
#    # def from_yaml(yaml_file):
#    #     ...
#
#    def run_steps(self, steps_list: tuple[Step, ...]) -> None:
#        '''
#        Wraps `self._run_steps` with pre- and post-run steps
#        '''
#
#        steps_to_run = []
#        for step in steps_list:
#            proceed = self._pre_run_step(step)
#            if proceed:
#                steps_to_run.append(step)
#
#        self._run_steps(tuple(steps_to_run))
#
#        for step in steps_to_run:
#            self._post_run_step(step)
#
#    def _pre_run_step(self, step: Step) -> bool:
#        """Perform operations that need to occur before a Step is run
#
#        :param step: The step to run
#        :return: Whether the step should be run
#        """
#
#        # Check that the step directory has been set
#        if step.directory is None:
#            raise ValueError(f'{step.__class__.__name__} directory must be set before running')
#
#        # Check that another step hasn't already been run in step.directory
#        if any([s.absolute_directory == step.absolute_directory for s in self.steps]):
#            raise ValueError(
#                f'A calculation has already been run in `{step.directory}`; this should not happen')
#
#        # If an output file already exists, check if the run completed successfully
#        if not self.from_scratch:
#            to_run = True
#            if isinstance(step, ImplementedCalc):
#                calc_file = step.directory / step.prefix
#
#                if calc_file.with_suffix(step.ext_out).is_file():
#                    loaded_step = self.load_old_calculator(step)
#
#                    if loaded_step.is_complete():
#                        to_run = False
#
#            elif isinstance(step, Process):
#                if step.is_complete():
#                    to_run = False
#                    step.load_outputs()
#            else:
#                raise ValueError(f'Unknown step type: {type(step)}')
#
#            if to_run:
#                # This and subsequent calculations should be performed from scratch
#                self.from_scratch = True
#            else:
#                self._step_skipped_message(step)
#                return False
#
#        # Remove the directory if it already exists
#        if step.directory.exists():
#            utils.remove(step.directory)
#
#        # Update postfix if relevant
#        if self.npool and isinstance(step, ImplementedCalc):
#            assert hasattr(step, 'command')
#            if isinstance(step.command, ParallelCommandWithPostfix):
#                step.command.postfix = f'-npool {self.npool}'
#
#        return True
#
#    @abstractmethod
#    def _run_steps(self, steps: tuple[Step, ...]) -> None:
#        """Run a list of steps, without doing anything else (other than printing messages)
#
#        Any pre- or post-processing of each step should be done in `_pre_run_step()` and
#        `_post_run_step()`.
#
#        :param steps: The steps to run
#        """
#        ...
#
#    def _post_run_step(self, step: Step) -> None:
#        pass
#
#    @abstractmethod
#    def load_old_calculator(self, calc: Calc) -> Calc:
#        ...
#
#
