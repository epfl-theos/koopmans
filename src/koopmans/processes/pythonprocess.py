from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Generic

from .process import InputModel, OutputModel, Process


class PythonProcess(Process[InputModel, OutputModel], ABC, Generic[InputModel, OutputModel]):
    """A class for running python code

    Roughly equivalent to CWL's `ExpressionTool`.
    """

    @abstractmethod
    def _run(self) -> None:
        """Run the python code for this process."""

    def run(self, **kwargs):
        """Run the process."""
        logger = logging.getLogger(__name__)
        logger.info(f'Running {self.directory}')
        self._pre_run()
        self._run(**kwargs)
        self._post_run()
