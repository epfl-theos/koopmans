import json
import os
import shutil
from abc import ABC
from pathlib import Path
from typing import Union

import numpy as np
from ase import Atoms

from koopmans.files import FilePointer
from koopmans.io import read_pkl
from koopmans.utils import chdir, symlink, warn

from ._utils import (benchmark_filename, metadata_filename,
                     recursively_find_files)


def atoms_eq(self, other):
    # Patching the Atoms class to compare positions and cell with np.allclose rather than strict equality
    if not isinstance(other, Atoms):
        return False
    a = self.arrays
    b = other.arrays
    return (len(self) == len(other) and
            np.allclose(a['positions'], b['positions']) and
            (a['numbers'] == b['numbers']).all() and
            np.allclose(self.cell, other.cell) and
            (self.pbc == other.pbc).all())


def write_mock_file(filename: Union[Path, str], written_by: str):
    filename = Path(filename)
    with chdir(filename.parent):
        with open(filename.name, 'w') as fd:
            json.dump({'filename': filename.name,
                       'written_to': str(filename.parent),
                       'written_by': written_by}, fd)


def mock_calculator__calculate(self):
    # Write the input file
    self.write_input(self.atoms)

    with chdir(self.directory):
        # By moving into the directory where the calculation was run, we ensure when we read in the settings that
        # paths are interpreted relative to this particular working directory
        calc = read_pkl(benchmark_filename(self), base_directory=self.directory)

    # Compare the settings
    for key in set(list(self.parameters.keys()) + list(calc.parameters.keys())):
        if key not in self.parameters:
            raise ValueError(f'Error in {self.prefix}: {key} is in benchmark but not in test')
        elif key not in calc.parameters:
            raise ValueError(f'Error in {self.prefix}: {key} is in test but not in benchmark')
        else:
            val = self.parameters[key]
            ref_val = calc.parameters[key]

            if isinstance(val, np.ndarray):
                val = val.tolist()

            if isinstance(ref_val, np.ndarray):
                ref_val = ref_val.tolist()

            if val != ref_val:
                raise ValueError(f'Error in {self.prefix}: {key} differs ({val} != {ref_val})')

    # Compare the atoms and the cell
    # TODO
    assert all(calc.atoms.pbc == self.atoms.pbc)

    with open(metadata_filename(self), 'r') as fd:
        metadata = json.load(fd)

    # Check that the right files exist (currently does not check the file contents)
    for path_str in metadata['input_files']:
        path = self.directory / path_str
        if not path.is_file():
            raise FileNotFoundError(f'File {path} does not exist')

    # Copy over the results
    self.results = calc.results
    if hasattr(calc, 'kpts'):
        self.kpts = calc.kpts

    # Create dummy output files
    main_output_file = Path(f'{self.directory}/{self.prefix}{self.ext_out}')
    if main_output_file.is_file():
        # If this input file exists already, this means that in a real calculation it is skipped with from_scratch,
        # so we won't generate any input files
        pass
    else:
        for path_str in metadata['output_files']:
            path = Path(path_str)
            directory = self.directory / path.parent
            filename = path.name

            if path.suffix == '.upf':
                if (self.directory / path).exists():
                    continue

                # Find the actual upf file and symlink it (so that we can read its contents if need be)
                [upf_file] = [os.path.relpath(x, directory)
                              for x in (self.directory / 'pseudopotentials').glob(filename)]

                with chdir(directory):
                    symlink(upf_file, filename)
            else:
                with chdir(directory):
                    # Create this (nested) directory if it does not exist and...
                    with open(filename, 'w') as f:
                        # ... write the file
                        json.dump({'filename': filename,
                                   'written_to': os.path.relpath(path, self.directory),
                                   'written_by': f'{self.directory / self.prefix}'}, f)


def mock_calculator_is_complete(self):
    return (self.directory / f'{self.prefix}{self.ext_out}').is_file()


def mock_calculator_check_code_is_installed(self):
    pass


def mock_calculator_read_results(self) -> None:
    raise AssertionError('A MockCalc should not attempt to read results')


def patch_generate_dos(calc_class, monkeypatch):
    # Patch the generate_dos method to first copy the DOS files from the benchmark directory to the test directory
    original_generate_dos = calc_class.generate_dos

    def generate_dos(self):
        # Fetch all the DOS files from the benchmark directory
        for f in benchmark_filename(self).parent.rglob('*.pdos*'):
            shutil.copy(f, self.directory / f.name)
        original_generate_dos(self)

    monkeypatch.setattr(calc_class, 'generate_dos', generate_dos)


def mock_process_run(self):
    # Load the inputs from file
    bench_process = read_pkl(benchmark_filename(self), base_directory=self.base_directory)
    assert self.inputs == bench_process.inputs
    self.outputs = bench_process.outputs

    for f in recursively_find_files([o for _, o in self.outputs]):
        if f.name in ['power_spectrum.npy']:
            src = benchmark_filename(self).parent / f.name
            assert src.exists()
            shutil.copy(src, f.name)
        else:
            write_mock_file(f, self.name)


def monkeypatch_mock(monkeypatch):
    from koopmans.calculators import (EnvironCalculator, KoopmansCPCalculator,
                                      KoopmansHamCalculator,
                                      KoopmansScreenCalculator, PhCalculator,
                                      ProjwfcCalculator, PW2WannierCalculator,
                                      PWCalculator, Wann2KCCalculator,
                                      Wann2KCPCalculator, Wannier90Calculator)
    from koopmans.processes.bin2xml import Bin2XMLProcess
    from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                                ConvertFilesFromSpin2To1)
    from koopmans.processes.merge_evc import MergeEVCProcess
    from koopmans.processes.power_spectrum import (
        ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)
    from koopmans.processes.ui import UnfoldAndInterpolateProcess
    from koopmans.processes.wannier import ExtendProcess, MergeProcess

    # Replace calculators with mock versions that obtain results from the database
    for c in [KoopmansCPCalculator, Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PhCalculator, PWCalculator, KoopmansCPCalculator, EnvironCalculator, Wann2KCCalculator, KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator]:
        monkeypatch.setattr(c, '_calculate', mock_calculator__calculate)
        monkeypatch.setattr(c, 'is_complete', mock_calculator_is_complete)
        monkeypatch.setattr(c, 'check_code_is_installed', mock_calculator_check_code_is_installed)
        monkeypatch.setattr(c, 'read_results', mock_calculator_read_results)

    patch_generate_dos(ProjwfcCalculator, monkeypatch)

    # for c in [PWCalculator, KoopmansHamCalculator]:
    #     monkeypatch.setattr(c, 'generate_band_structure', mock_generate_band_structure)

    # Processes
    for p in [ExtractCoefficientsFromXMLProcess, ComputePowerSpectrumProcess, Bin2XMLProcess, ConvertFilesFromSpin1To2, ConvertFilesFromSpin2To1, ExtendProcess, MergeProcess, UnfoldAndInterpolateProcess, MergeEVCProcess]:
        monkeypatch.setattr(p, '_run', mock_process_run)

    # Patch the Atoms class
    monkeypatch.setattr(Atoms, '__eq__', atoms_eq)
