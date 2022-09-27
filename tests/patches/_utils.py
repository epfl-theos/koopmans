import os
from pathlib import Path
from typing import Set, Tuple

from koopmans import calculators


def benchmark_filename(calc: calculators.CalculatorExt) -> Path:
    base_directory = Path(__path__).parents[2]
    benchmark_dir = base_directory / 'tests' / 'benchmarks'
    if base_directory / 'tests' / 'tmp' in calc.directory.parents:
        parent = base_directory / 'tests' / 'tmp'
    else:
        parent = base_directory
    benchmark_name = calc.directory.relative_to(parent) / calc.prefix
    return benchmark_dir / (str(benchmark_name).replace(os.path.sep, '-') + '.json')


def metadata_filename(calc: calculators.CalculatorExt) -> Path:
    benchmark_path = benchmark_filename(calc)
    return benchmark_path.with_name(benchmark_path.name.replace('.json', '_metadata.json'))


def find_subfiles_of_calc(calc: calculators.CalculatorExt) -> Set[Tuple[Path, float]]:
    files = find_subfiles_of_dir(calc.directory)
    if 'outdir' in calc.parameters.valid:
        files = files | find_subfiles_of_dir(calc.parameters.outdir)
    return files


def find_subfiles_of_dir(base_dir: Path) -> Set[Tuple[Path, float]]:
    files: Set[Tuple[Path, float]]
    if base_dir.exists():
        files = set([(x, x.stat().st_mtime) for x in base_dir.rglob('*') if x.is_file()])
    else:
        files = set([])
    return files
