import os
from pathlib import Path
from typing import Set, Tuple

from koopmans.files import FilePointer


def benchmark_filename(obj) -> Path:
    base_directory = Path(__file__).parents[3]
    benchmark_dir = base_directory / 'tests' / 'benchmarks'
    assert obj.directory is not None
    abs_obj_directory = obj.directory.resolve()
    if base_directory / 'tests' / 'tmp' in abs_obj_directory.parents:
        parent = base_directory / 'tests' / 'tmp'
    else:
        parent = base_directory
    name = getattr(obj, 'prefix', getattr(obj, 'name', None))
    assert isinstance(name, str)
    benchmark_name = abs_obj_directory.relative_to(parent) / name
    return benchmark_dir / benchmark_name.with_suffix('.pkl')


def metadata_filename(calc) -> Path:
    benchmark_path = benchmark_filename(calc)
    return benchmark_path.with_name(benchmark_path.name.replace('.pkl', '_metadata.json'))


def find_subfiles_of_calc(calc) -> Set[Tuple[Path, float]]:
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


def recursively_find_files(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from recursively_find_files(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from recursively_find_files(v)
    elif isinstance(obj, FilePointer):
        yield obj.aspath()
    elif isinstance(obj, Path):
        yield obj
    return
