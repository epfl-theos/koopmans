import os
from pathlib import Path
from koopmans import base_directory


def benchmark_filename(calc):
    benchmark_dir = base_directory / 'tests' / 'benchmarks'
    if base_directory / 'tests' / 'tmp' in calc.directory.parents:
        parent = base_directory / 'tests' / 'tmp'
    else:
        parent = base_directory
    benchmark_name = calc.directory.relative_to(parent) / calc.prefix
    return benchmark_dir / (str(benchmark_name).replace(os.path.sep, '-') + '.json')
