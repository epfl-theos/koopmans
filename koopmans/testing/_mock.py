from abc import ABC
import json
import numpy as np
import os
from pathlib import Path
from typing import List, Union, Optional
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator, Calc
from koopmans.workflows import WannierizeWorkflow, KoopmansDSCFWorkflow
from koopmans import utils, projections
from koopmans.io import read_kwf as read_encoded_json
from ._utils import benchmark_filename, metadata_filename


def write_mock_file(filename: Union[Path, str], written_by: str):
    filename = Path(filename)
    with utils.chdir(filename.parent):
        with open(filename.name, 'w') as fd:
            json.dump({'filename': filename.name,
                       'written_to': str(filename.parent),
                       'written_by': written_by}, fd)


class MockCalc(ABC):
    def calculate(self):
        # Write the input file
        self.write_input(self.atoms)

        with utils.chdir(self.directory):
            # By moving into the directory where the calculation was run, we ensure when we read in the settings that
            # paths are interpreted relative to this particular working directory
            with open(benchmark_filename(self), 'r') as fd:
                calc = read_encoded_json(fd)

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

        # Check that the right files exist
        # TODO

        # Copy over the results
        self.results = calc.results

        # Create dummy output files
        main_output_file = Path(f'{self.directory}/{self.prefix}{self.ext_out}')
        if main_output_file.is_file():
            # If this input file exists already, this means that in a real calculation it is skipped with from_scratch,
            # so we won't generate any input files
            pass
        else:
            with open(metadata_filename(self), 'r') as fd:
                metadata = json.load(fd)
            for path_str in metadata['output_files']:
                path = Path(path_str)
                directory = (self.directory / path.parent).resolve()
                filename = path.name

                with utils.chdir(directory):
                    # Create this (nested) directory if it does not exist and...
                    with open(filename, 'w') as f:
                        # ... write the file
                        json.dump({'filename': filename,
                                   'written_to': os.path.relpath(path, self.directory),
                                   'written_by': f'{self.directory / self.prefix}'}, f)

    def is_complete(self):
        return (self.directory / f'{self.prefix}{self.ext_out}').is_file()

    def check_code_is_installed(self):
        # Don't check if the code is installed
        return

    def check_convergence(self):
        # Monkeypatched version of check_convergence to avoid any convergence check
        return

    def read_results(self) -> None:
        raise AssertionError('A MockCalc should not attempt to read results')


class MockKoopmansCPCalculator(MockCalc, KoopmansCPCalculator):
    # Monkeypatched KoopmansCPCalculator class which never actually calls kcp.x
    pass


class MockEnvironCalculator(MockCalc, EnvironCalculator):
    # Monkeypatched EnvironCalculator class which never actually calls pw.x
    pass


class MockPWCalculator(MockCalc, PWCalculator):
    # Monkeypatched PWCalculator class which never actually calls pw.x
    def generate_band_structure(self):
        pass


class MockWannier90Calculator(MockCalc, Wannier90Calculator):
    # Monkeypatched Wannier90Calculator class which never actually calls wannier90.x
    pass


class MockPW2WannierCalculator(MockCalc, PW2WannierCalculator):
    # Monkeypatched PW2WannierCalculator class which never actually calls pw2wannier90.x
    pass


class MockWann2KCPCalculator(MockCalc, Wann2KCPCalculator):
    # Monkeypatched Wann2KCPCalculator class which never actually calls wann2kcp.x
    pass


class MockUnfoldAndInterpolateCalculator(MockCalc, UnfoldAndInterpolateCalculator):
    def write_results(self):
        # Do nothing when it goes to write out the results (because this relies on self.bg, among others)
        return


class MockWann2KCCalculator(MockCalc, Wann2KCCalculator):
    pass


class MockKoopmansScreenCalculator(MockCalc, KoopmansScreenCalculator):
    pass


class MockKoopmansHamCalculator(MockCalc, KoopmansHamCalculator):
    def generate_band_structure(self):
        pass


class MockProjwfcCalculator(MockCalc, ProjwfcCalculator):
    def generate_dos(self):
        return


class MockWorkflow:

    calculations: List[Calc]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parameters.from_scratch = True

    def load_old_calculator(self, calc: Calc) -> bool:
        # Load old calculators by looking through the workflow's list of previous calculations
        # (During the test suite, the only scenario where we want to reload an old calculation will arise
        # because we did it earlier in the same workflow)

        # Load the dummy output file
        calc_fname = calc.directory / f'{calc.prefix}{calc.ext_out}'
        with open(calc_fname, 'r') as fd:
            output_file_info = json.load(fd)

        # Find out the calculation that generated the output file
        written_by = Path(output_file_info['written_by'])
        directory, prefix = written_by.parent, written_by.name
        matches = [c for c in self.calculations if c.directory == directory and c.prefix == prefix]

        # Copy that calculation into the record of all calculations
        if len(matches) == 1:
            # Fetch the results from the match
            calc.results = matches[0].results
            self.calculations.append(calc)
        elif len(matches) == 0:
            raise ValueError(f'Could not find a calculator matching {calc.directory}/{calc.prefix}')
        else:
            raise ValueError(f'Found multiple calculators for {calc.directory}/{calc.prefix}')

        return calc.is_complete()


class MockKoopmansDSCFWorkflow(MockWorkflow, KoopmansDSCFWorkflow):
    pass


class MockWannierizeWorkflow(MockWorkflow, WannierizeWorkflow):

    @staticmethod
    def _merge_wannier_files(dirs_in: List[Path], dir_out: Path, fname: str):
        for dir_in in dirs_in:
            fname_in = dir_in.resolve() / fname
            assert fname_in.exists()
        write_mock_file(dir_out / fname, 'workflow')

    @staticmethod
    def merge_wannier_hr_files(dirs_in: List[Path], dir_out: Path, prefix: str):
        MockWannierizeWorkflow._merge_wannier_files(dirs_in, dir_out, prefix + '_hr.dat')

    @staticmethod
    def merge_wannier_u_files(dirs_in: List[Path], dir_out: Path, prefix: str):
        MockWannierizeWorkflow._merge_wannier_files(dirs_in, dir_out, prefix + '_u.mat')

    @staticmethod
    def merge_wannier_centres_files(dirs_in: List[Path], dir_out: Path, prefix: str):
        MockWannierizeWorkflow._merge_wannier_files(dirs_in, dir_out, prefix + '_centres.xyz')

    def extend_wannier_u_dis_file(self, block: List[projections.ProjectionBlock], prefix: str = 'wann'):
        raise NotImplementedError()
