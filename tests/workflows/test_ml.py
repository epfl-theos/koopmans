"""
Test for to- and fromdict functions of the MLFittingWorkflow class
"""
import numpy as np

from koopmans import utils
from koopmans.calculators import KoopmansCPCalculator
from koopmans.io import read, write
from koopmans.workflows import MLFittingWorkflow, SinglepointWorkflow, Workflow


def test_mlfitting_write_then_read_empty_from_dict(tmpdir, water):
    with utils.chdir(tmpdir):
        dummy_calculator = KoopmansCPCalculator(**water)
        water['atoms'].calc = None
        mlfit_out = MLFittingWorkflow(calc_that_produced_orbital_densities=dummy_calculator, **water)
        mlfit_dict = mlfit_out.todict()
        mlfit_in = MLFittingWorkflow.fromdict(mlfit_dict)

        assert mlfit_out == mlfit_in


def test_mlfitting_write_then_read_non_empty_from_dict(tmpdir, water):
    with utils.chdir(tmpdir):
        dummy_calculator = KoopmansCPCalculator(**water)
        water['atoms'].calc = None
        mlfit_out = MLFittingWorkflow(calc_that_produced_orbital_densities=dummy_calculator, **water)
        mlfit_out.input_vectors_for_ml['test_data'] = np.array([[2.7, 5.9, 4.9], [4.2, 4.9, 6.9]])
        mlfit_dict = mlfit_out.todict()
        mlfit_in = MLFittingWorkflow.fromdict(mlfit_dict)

        assert mlfit_out == mlfit_in
