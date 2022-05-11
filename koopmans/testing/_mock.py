from abc import ABC, abstractproperty
import json
import numpy as np
import os
from pathlib import Path
from typing import List
from ase.dft.kpoints import BandPath
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator
from koopmans import utils
from koopmans.io import read_kwf as read_encoded_json
from ._utils import benchmark_filename


class MockCalc(ABC):
    def calculate(self):
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
    # Monkeypatched KoopmansCPCalculator class which never actually calls kcp.x

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


class MockEnvironCalculator(MockCalc, EnvironCalculator):
    # Monkeypatched EnvironCalculator class which never actually calls pw.x

    @property
    def output_files(self) -> List[Path]:
        files = []
        if 'kpts' in self.parameters:
            assert self.parameters.kpts == [
                1, 1, 1], 'Have not implemented dummy environ calculations for kpts != [1, 1, 1]'

            i_kpoints = range(1)
            files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/wfc{i}.dat' for i in i_kpoints]
            if self.parameters.nosym:
                files += [f'{self.parameters.outdir}/wfc{i}.dat' for i in i_kpoints]
        files += [f'{self.parameters.outdir}/{self.parameters.prefix}.xml']
        return [Path(f) for f in files]

    @property
    def input_files(self) -> List[Path]:
        # Not yet implemented
        return []


class MockPWCalculator(MockCalc, PWCalculator):
    # Monkeypatched PWCalculator class which never actually calls pw.x

    @property
    def output_files(self) -> List[Path]:
        files = []
        if 'kpts' in self.parameters:
            if isinstance(self.parameters.kpts, BandPath):
                n_kpoints = len(self.parameters.kpts.kpts)
            elif self.parameters.nosym:
                n_kpoints = np.prod(self.parameters['kpts']) + 1
            else:
                n_kpoints = len(BandPath(self.parameters['kpts'], self.atoms.cell).kpts)
            i_kpoints = range(1, n_kpoints + 1)

            files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/wfc{i}.dat' for i in i_kpoints]
        files += [f'{self.parameters.outdir}/{self.parameters.prefix}.xml']
        files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/{f}'
                  for f in ['data-file-schema.xml', 'charge-density.dat']]
        return [Path(f) for f in files]

    @property
    def input_files(self) -> List[Path]:
        files = []
        if self.parameters.calculation == 'nscf':
            files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/{f}'
                      for f in ['data-file-schema.xml', 'charge-density.dat']]
        return [Path(f) for f in files]

    def generate_band_structure(self):
        pass


class MockWannier90Calculator(MockCalc, Wannier90Calculator):
    # Monkeypatched Wannier90Calculator class which never actually calls wannier90.x

    @property
    def output_files(self) -> List[Path]:
        if '-pp' in self.command.flags:
            files = [self.directory / f'{self.prefix}.nnkp']
        else:
            suffixes = ['.chk', '_wsvec.dat', '_hr.dat']
            if self.parameters.get('write_u_matrices', False):
                suffixes += ['_u.mat', '_u_dis.mat']
            if self.parameters.get('write_xyz', False):
                suffixes += ['_centres.xyz']
            files = [self.directory / f'{self.prefix}{suffix}' for suffix in suffixes]
        return [f for f in files]

    @property
    def input_files(self) -> List[Path]:
        if '-pp' in self.command.flags:
            files = []
        else:
            files = [f'{self.directory}/{self.prefix}{suffix}' for suffix in ['.eig', '.mmn', '.amn']]
        return [Path(f) for f in files]


class MockPW2WannierCalculator(MockCalc, PW2WannierCalculator):
    # Monkeypatched PW2WannierCalculator class which never actually calls pw2wannier90.x

    @property
    def output_files(self) -> List[Path]:
        files = [self.directory / f'{self.parameters.seedname}{suffix}' for suffix in ['.eig', '.mmn', '.amn']]
        return [Path(f) for f in files]

    @property
    def input_files(self) -> List[Path]:
        i_kpoints = range(1, np.prod(self.parameters.kpts) + 1)
        files = [f'{self.parameters.outdir}/{self.parameters.prefix}.save/wfc{i}.dat' for i in i_kpoints]
        files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/{f}'
                  for f in ['data-file-schema.xml', 'charge-density.dat']]
        files.append(f'{self.directory}/{self.parameters.seedname}.nnkp')
        return [Path(f) for f in files]


class MockWann2KCPCalculator(MockCalc, Wann2KCPCalculator):
    # Monkeypatched Wann2KCPCalculator class which never actually calls wann2kcp.x

    @property
    def output_files(self) -> List[Path]:
        if self.parameters.wan_mode == 'wannier2kcp':
            files = [self.directory / fname for fname in ['evcw1.dat', 'evcw2.dat']]
        else:
            files = [self.directory / fname for i in [1, 2]
                     for fname in [f'evc_occupied{i}.dat', f'evc0_empty{i}.dat']]
        return [Path(f) for f in files]

    @property
    def input_files(self) -> List[Path]:
        i_kpoints = range(1, np.prod(self.parameters.kpts) + 1)
        files = [f'{self.parameters.outdir}/{self.parameters.prefix}.save/wfc{i}.dat' for i in i_kpoints]
        files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/{f}'
                  for f in ['data-file-schema.xml', 'charge-density.dat']]
        if self.parameters.wan_mode == 'wannier2kcp':
            files.append(f'{self.directory}/{self.parameters.seedname}.nnkp')
            files.append(f'{self.directory}/{self.parameters.seedname}.chk')
        return [Path(f) for f in files]


class MockUnfoldAndInterpolateCalculator(MockCalc, UnfoldAndInterpolateCalculator):
    def write_results(self):
        # Do nothing when it goes to write out the results (because this relies on self.bg, among others)
        return

    # We don't ever construct self.bg during a mock calculation, but we refer to it later, so give it a dummy value
    @property
    def bg(self):
        return None

    @bg.setter
    def bg(self, value):
        return

    @property
    def output_files(self):
        return []

    @property
    def input_files(self):
        return []


class MockWann2KCCalculator(MockCalc, Wann2KCCalculator):
    @property
    def input_files(self) -> List[Path]:
        return []

    @property
    def output_files(self) -> List[Path]:
        return []


class MockKoopmansScreenCalculator(MockCalc, KoopmansScreenCalculator):
    @property
    def input_files(self) -> List[Path]:
        return []

    @property
    def output_files(self) -> List[Path]:
        return []


class MockKoopmansHamCalculator(MockCalc, KoopmansHamCalculator):
    @property
    def input_files(self) -> List[Path]:
        return []

    @property
    def output_files(self) -> List[Path]:
        return []

    def generate_band_structure(self):
        pass


class MockProjwfcCalculator(MockCalc, ProjwfcCalculator):
    @property
    def input_files(self) -> List[Path]:
        return []

    @property
    def output_files(self) -> List[Path]:
        return []

    def generate_dos(self):
        return
