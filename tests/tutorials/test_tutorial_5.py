from pathlib import Path

import pytest
from deepdiff import DeepDiff

from koopmans.io import read
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json
from koopmans.utils import chdir

tutorial_dir = Path(__file__).parents[2] / 'tutorials' / 'tutorial_5'
benchmark_dir = Path(__file__).parents[1] / 'benchmarks'


@pytest.mark.tutorials
def test_run_trajectory_ml(tutorial_patch, tmpdir, pytestconfig):
    with chdir(tutorial_dir / 'tutorial_5a'):
        # run the tutorial
        wf = read('h2o_trajectory_ml.json')
        wf.run()

        # obtain all screening parameters
        alphas = wf.all_alphas

        # check that the screening parameters match the reference solution
        benchmark_file = benchmark_dir / 'test_run_trajectory_ml.json'

        if pytestconfig.getoption('generate_benchmark'):
            # Write the power spectrum tensor to file
            with open(benchmark_file, 'w') as fd:
                write_encoded_json(alphas, fd)
        else:
            # Compare with the power spectrum on file
            with open(benchmark_file, 'r') as fd:
                alphas_ref = read_encoded_json(fd)
            assert DeepDiff(alphas, alphas_ref, significant_digits=8,
                            number_format_notation='e', ignore_numeric_type_changes=True) == {}


@pytest.mark.tutorials
def test_run_convergence_ml(tutorial_patch, tmpdir, pytestconfig, sys2file):
    with chdir(tutorial_dir / 'tutorial_5b'):
        # run the tutorial
        wf = read('h2o_convergence_ml.json')
        wf.run()

        # obtain the result dictionary
        results = wf.result_dict

        # the eigenvalues can not be compared since they would require a final Koopmans calculation
        results['spin_0'].pop('evs')
        results['spin_1'].pop('evs')

        # check that all entries from the result dictionary match the reference solution
        benchmark_file = benchmark_dir / 'test_run_convergence_ml.json'

        if pytestconfig.getoption('generate_benchmark'):
            with open(benchmark_file, 'w') as fd:
                write_encoded_json(results, fd)
        else:
            with open(benchmark_file, 'r') as fd:
                results_ref = read_encoded_json(fd)
            assert DeepDiff(results, results_ref, significant_digits=8,
                            number_format_notation='e', ignore_numeric_type_changes=True) == {}
