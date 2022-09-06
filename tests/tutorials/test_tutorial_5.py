from pathlib import Path

import pytest
from deepdiff import DeepDiff

from koopmans import base_directory
from koopmans.io import read
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write
from koopmans.io import write_kwf as write_encoded_json
from koopmans.utils import chdir

tutorial_dir = base_directory / 'tutorials' / 'tutorial_5'


@pytest.mark.tutorials
def test_run_trajectory_ml(tutorial_patch, tmpdir, pytestconfig):
    with chdir(tutorial_dir / 'tutorial_5a'):
        wf = read('h2o_trajectory_ml.json')

        wf.run()
        alphas = wf.all_alphas

        benchmark_dir = base_directory / 'tests' / 'benchmarks'
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
        wf = read('h2o_convergence_ml.json')
        wf.run()

        results = wf.result_dict

        benchmark_dir = base_directory / 'tests' / 'benchmarks'
        benchmark_file = benchmark_dir / 'test_run_convergence_ml.json'

        if pytestconfig.getoption('generate_benchmark'):
            # Write the power spectrum tensor to file
            with open(benchmark_file, 'w') as fd:
                write_encoded_json(results, fd)
        else:
            # Compare with the power spectrum on file
            with open(benchmark_file, 'r') as fd:
                results_ref = read_encoded_json(fd)
            assert DeepDiff(results, results_ref, significant_digits=8,
                            number_format_notation='e', ignore_numeric_type_changes=True) == {}
