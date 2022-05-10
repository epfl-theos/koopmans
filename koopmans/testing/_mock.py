import os
from typing import List
from abc import ABC, abstractproperty
import json
from koopmans.calculators import KoopmansCPCalculator
from koopmans import utils
from pathlib import Path
from koopmans.io import read_kwf as read_encoded_json
from ._utils import benchmark_filename


class MockCalc(ABC):
    def _ase_calculate(self):
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
            for path in self.output_files + [main_output_file]:
                directory = path.parent
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

    @abstractproperty
    def input_files(self) -> List[Path]:
        ...

    @abstractproperty
    def output_files(self) -> List[Path]:
        ...


class MockKoopmansCPCalculator(MockCalc, KoopmansCPCalculator):
    @property
    def __files(self) -> List[Path]:
        files = [fname for ispin in range(1, self.parameters.nspin + 1) for fname in [f'evc0{ispin}.dat',
                                                                                      f'evc{ispin}.dat',
                                                                                      f'evcm{ispin}.dat',
                                                                                      f'hamiltonian{ispin}.xml',
                                                                                      f'eigenval{ispin}.xml',
                                                                                      f'lambda0{ispin}.dat',
                                                                                      f'lambdam{ispin}.dat']]
        for ispin in range(self.parameters.nspin):
            if self.has_empty_states(ispin):
                files += [fname for fname in [f'evc0_empty{ispin + 1}.dat', f'evc_empty{ispin + 1}.dat']]
        return [Path(f) for f in files]

    @property
    def output_files(self) -> List[Path]:
        files = self.__files
        if self.parameters.print_wfc_anion:
            for ispin in range(self.parameters.nspin):
                files.append(Path(f'evcfixed_empty{ispin + 1}.dat'))
        return [self.parameters.outdir
                / Path(f'{self.parameters.prefix}_{self.parameters.ndw}.save/K00001/{fname}') for fname in files]

    @property
    def input_files(self) -> List[Path]:
        files = self.__files
        if self.parameters.restart_from_wannier_pwscf:
            for ispin in range(self.parameters.nspin):
                files.append(Path(f'evc_occupied{ispin + 1}.dat'))
        return [self.parameters.outdir
                / Path(f'{self.parameters.prefix}_{self.parameters.ndr}.save/K00001/{fname}') for fname in files]
