import json
import os
import shutil
from pathlib import Path
from typing import List

from koopmans.calculators import (EnvironCalculator, KoopmansCPCalculator,
                                  KoopmansHamCalculator,
                                  KoopmansScreenCalculator, PhCalculator,
                                  ProjwfcCalculator, PW2WannierCalculator,
                                  PWCalculator, Wann2KCCalculator,
                                  Wann2KCPCalculator, Wannier90Calculator)
from koopmans.io import write_pkl
from koopmans.processes.bin2xml import Bin2XMLProcess
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.processes.power_spectrum import (
    ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)
from koopmans.processes.ui import UnfoldAndInterpolateProcess
from koopmans.processes.wannier import ExtendProcess, MergeProcess

from ._utils import (benchmark_filename, find_subfiles_of_calc,
                     metadata_filename, recursively_find_files)

base_directory = Path(__file__).parents[3]


def patch_calculator(c, monkeypatch):

    unpatched_calculate = c._calculate

    def _calculate(self):
        # Before running the calculation, make a list of the files that exist
        files_before = find_subfiles_of_calc(self)

        # Run the calculation
        unpatched_calculate(self)

        # Store the calculator as a pickle file to use as a benchmark, Temporary wiping the parent attribute so that
        # the entire workflow doesn't get pickled
        self.results.pop('walltime', None)
        filename = benchmark_filename(self)
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True)

        self.parent, parent = None, self.parent
        self.linked_files, linked_files = [], self.linked_files
        write_pkl(self, filename, base_directory=base_directory)
        self.parent = parent
        self.linked_files = linked_files

        # After running the calculation, make a new list of the files, and then work out which files have been
        # modified by the calculation
        files_after = find_subfiles_of_calc(self)
        modified_files: List[str] = sorted([str(os.path.relpath(x[0], self.directory))
                                           for x in files_after - files_before])
        # Exclude the input file, since we don't want to create a dummy version of this for the mock calc (it is
        # already written in full)
        modified_files.remove(self.prefix + self.ext_in)

        # Exclude .wfc* files because these depend on the parallelism
        if isinstance(self, (PWCalculator, Wann2KCCalculator, PhCalculator)):
            assert isinstance(self.parameters.outdir, Path)
            modified_files = [f for f in modified_files if not Path(f).suffix.startswith('.wfc')]

        # For pDOS files, copy the files explicitly (because we need to parse them later)
        for pdos_file in [Path(f) for f in modified_files if Path(f).suffix.startswith('.pdos')]:
            shutil.copy(self.directory / pdos_file, benchmark_filename(self).parent / pdos_file)
        modified_files = [f for f in modified_files if not Path(f).suffix.startswith('.pdos')]

        # Write out the modified files to the "_metadata.pkl" file
        sanitised_files_before = sorted([str(os.path.relpath(x[0], self.directory)) for x in files_before])
        sanitised_files_before = [f for f in sanitised_files_before if not Path(f).suffix.startswith('.wfc')]
        fname = metadata_filename(self)
        with open(fname, 'w') as fd:
            json.dump({'input_files': sanitised_files_before, 'output_files': sorted(modified_files)}, fd)

    monkeypatch.setattr(c, '_calculate', _calculate)


def patch_process(p, monkeypatch):

    unpatched_run = p._run

    def _run(self):
        unpatched_run(self)

        # Write the benchmark to file
        filename = benchmark_filename(self)
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True)
        # Temporarily wipe the parent attribute so that the entire workflow doesn't get pickled
        self.parent, parent = None, self.parent
        write_pkl(self, filename, base_directory=base_directory)
        self.parent = parent

        # Copy over all files that are outputs of the process that need to be read
        for filepath in recursively_find_files([o for _, o in self.outputs]):
            if filepath.name in ['power_spectrum.npy']:
                shutil.copy(filepath, benchmark_filename(self).parent / filepath.name)

    # Patching the absolute_directory property
    unpatched_absolute_directory = p.absolute_directory

    def absolute_directory(self) -> Path:
        if self.parent is None:
            # Because we wipe parents when storing benchmarks (see above), this prevents us from being able to construct
            # an absolute directory to locate files. Usually, this would raise an error. For the purposes of the test suite,
            # instead simply use the base directory of the repo
            return Path().resolve().relative_to(base_directory)
        else:
            # Default behavior
            return unpatched_absolute_directory.__get__(self)

    monkeypatch.setattr(p, '_run', _run)
    monkeypatch.setattr(p, 'absolute_directory', property(absolute_directory))


def monkeypatch_bench(monkeypatch):
    for c in [EnvironCalculator, KoopmansCPCalculator, KoopmansHamCalculator, KoopmansScreenCalculator,
              PhCalculator, ProjwfcCalculator, PW2WannierCalculator, PWCalculator, Wann2KCCalculator,
              Wann2KCPCalculator, Wannier90Calculator]:
        patch_calculator(c, monkeypatch)

    for p in [Bin2XMLProcess, ComputePowerSpectrumProcess, ConvertFilesFromSpin1To2, ConvertFilesFromSpin2To1,
              ExtractCoefficientsFromXMLProcess, UnfoldAndInterpolateProcess, ExtendProcess, MergeProcess]:
        patch_process(p, monkeypatch)
