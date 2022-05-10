import os
from pathlib import Path
from koopmans import __path__ as koopmans_base_directory

koopmans_base_directory = Path(koopmans_base_directory[0]).parent


def benchmark_filename(calc):
    benchmark_dir = Path(__file__).parent / 'benchmarks'
    benchmark_name = calc.directory.relative_to(koopmans_base_directory) / calc.prefix
    return benchmark_dir / (str(benchmark_name).replace(os.path.sep, '-') + '.json')
