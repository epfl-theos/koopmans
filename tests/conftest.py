'''
Functionality for testing koopmans code with pytest

Written by Edward Linscott, Dec 2020
'''

import os
import json
import numpy as np
import copy
from typing import List, Union, Dict, Tuple
from pathlib import Path
import pytest
from ase.calculators.calculator import CalculationFailed
from ase.dft.kpoints import BandPath
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, PWCalculator, KoopmansCPCalculator, \
    EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, KoopmansScreenCalculator, \
    KoopmansHamCalculator
from koopmans.io import read
from koopmans.io import read_kwf as read_encoded_json
from koopmans.io import write_kwf as write_encoded_json
from koopmans import utils


# A hard and a soft tolerance for checking floats
default_tolerances = {'alpha': (2e-3, 2e-5),
                      'homo_energy': (2e-3, 2e-5),
                      'lumo_energy': (2e-3, 2e-5),
                      'self-hartree': (2e-3, 2e-5),
                      'array': (2e-2, 2e-4),
                      'default': (1e-4, 1e-7)}


# Set this to true to construct a benchmark file
construct_exceptions = False


class WorkflowTest:

    def __init__(self, json_input: Union[str, Path], capsys: pytest.CaptureFixture[str], mock: bool = False,
                 tolerances: Dict[str, Tuple[float, float]] = default_tolerances):

        if isinstance(json_input, str):
            json_input = Path(json_input)

        self.directory = json_input.parent.resolve()
        self.json = json_input.name
        self.tests_directory = self.directory.parent
        self.subdirectory = self.directory.name
        self.mock = mock
        self.capsys = capsys
        self.tolerances = tolerances

        # Fetch the benchmarks
        with open('tests/benchmarks.json', 'r') as f:
            benchmarks = read_encoded_json(f)

        self.benchmark = {k: v for k, v in benchmarks.items() if self.subdirectory in k}

    def run(self):
        # Move into the test directory
        with utils.chdir(self.directory):

            # Load the input file
            workflow = read(self.json)

            # Give it access to the benchmark
            workflow.benchmark = self.benchmark

            # Make sure QC printing is turned on
            workflow.parameters.print_qc = True

            try:
                # Run the workflow
                workflow.run()
            finally:
                # Write the output, even if the workflow gave an error
                out, err = self.capsys.readouterr()
                stdout = self.json.replace('.json', '.stdout')
                with open(stdout, 'w') as f:
                    f.write(out)
                    f.write(err)

            if not self.mock:
                # Print out the contents of the QE i/o files
                utils.system_call(f'cat_workflow_io.sh >> {stdout}')

                # Check numerical results against the benchmarks
                self.check_qc_results(workflow)

    def check_qc_result(self, result, ref_result, result_name, tols):
        # Compare the calculated result to the reference result
        if isinstance(ref_result, float):
            diff = result - ref_result
            message = f'{result_name} = {result:.5f} differs from benchmark {ref_result:.5f} by ' \
                      f'{diff:.2e}'
            if abs(diff) > tols[0]:
                return {'kind': 'error', 'message': message}
            elif abs(diff) > tols[1]:
                return {'kind': 'warning', 'message': message}
        elif isinstance(ref_result, np.ndarray):
            # For arrays, perform a mixed error test. If Delta = |x - x_ref| then in the limit of large Delta,
            # then this reduces to testing relative error, whereas in the limit of small Delta it reduces to
            # testing absolute error. We use 0.1*max(ref_result) as a reference scale factor.

            # Sanitise the input
            result = result.flatten()
            ref_result = ref_result.flatten()

            # Calculate the absolute difference
            abs_diffs = np.abs(result - ref_result)

            # Calculate the mixed difference
            scale_factor = 0.1 * np.max(np.abs(ref_result))
            if scale_factor == 0.0:
                scale_factor = 1.0
            mixed_diffs = abs_diffs / (scale_factor + np.abs(ref_result))

            # Locate the datapoint with the larges mixed difference
            i_max = np.argmax(mixed_diffs)
            mixed_diff = mixed_diffs[i_max]
            abs_diff = abs_diffs[i_max]

            # Generate an error message if necessary
            message = f'{result_name}[{i_max}] = {result[i_max]:.5f} differs from benchmark ' \
                      f'{ref_result[i_max]:.5f} by {abs_diff:.2e}'
            if mixed_diff > tols[0]:
                return {'kind': 'error', 'message': message}
            elif mixed_diff > tols[1]:
                return {'kind': 'warning', 'message': message}
        else:
            if result != ref_result:
                message = f'{result_name} = {result} differs from benchmark {ref_result}'
                return {'kind': 'error', 'message': message}

        return None

    def check_qc_results(self, workflow):

        # Loop through the stored QC results of all the calculations
        log = {}
        for calc in workflow.calculations:

            # Skip checking n-1 calculations since these are known to be unreliable
            if 'n-1' in calc.prefix:
                continue

            # Skip calculations in the base postproc directory because this corresponds to the calculation where we
            # simply merge the two previous calculations together
            if calc.directory.name == 'postproc':
                continue

            # Prepare the log dict entry for this calc
            calc_relpath = os.path.relpath(f'{calc.directory}/{calc.prefix}', self.tests_directory)
            log[calc_relpath] = []

            # If there has been no monkey-patching, calc will not have a benchmark attribute yet, so add that
            assert calc_relpath in workflow.benchmark, f'Could not find an entry for {calc_relpath} in benchmarks.json'
            calc.benchmark = workflow.benchmark[calc_relpath]

            # Loop through results that require checking
            for result_name, result in calc.qc_results.items():

                # Fetch the corresponding benchmark
                if result_name in ['band structure', 'dos']:
                    if result_name == 'band structure':
                        result = result.energies[0].flatten()
                        ref_result = np.array(calc.benchmark['results'][result_name].energies[0].flatten())
                    else:
                        result = result.get_dos()
                        ref_result = np.array(calc.benchmark['results'][result_name].get_dos())
                    assert result.shape == ref_result.shape, 'Array shape mismatch between result and benchmark'
                    tols = self.tolerances['array']
                else:
                    assert result_name in calc.benchmark['results']
                    ref_result = calc.benchmark['results'][result_name]
                    tols = self.tolerances.get(result_name, self.tolerances['default'])

                message = self.check_qc_result(result, ref_result, result_name, tols)
                if message is not None:
                    log[calc_relpath].append(message)

            # Also check the alphas at this stage
            if 'alphas' in calc.benchmark or hasattr(calc, 'alphas'):
                assert 'alphas' in calc.benchmark, f'{calc.directory/calc.prefix} benchmark is missing alphas'
                assert hasattr(calc, 'alphas'), f'{calc.directory/calc.prefix} is missing alphas'

                # Sanitise the alphas (flatten the arrays and convert to numpy arrays)
                if isinstance(calc.alphas[0], list):
                    alphas = np.array([x for row in calc.alphas for x in row])
                    ref_alphas = np.array([x for row in calc.benchmark['alphas'] for x in row])
                else:
                    alphas = np.array(calc.alphas)
                    ref_alphas = np.array(calc.benchmark['alphas'])
                message = self.check_qc_result(alphas, ref_alphas, 'alphas', self.tolerances['alpha'])
                if message is not None:
                    log[calc_relpath].append(message)

        errors = []
        warnings = []
        with open(self.json.replace('.json', '.stdout'), 'a') as f:
            for calc_relpath, messages in log.items():
                if len(messages) > 0:
                    # Construct a list of warnings
                    warning_list = [f'  {m["message"]}' for m in messages if m['kind'] == 'warning']
                    if len(warning_list) > 0:
                        warnings.append(calc_relpath + '\n' + '\n'.join(warning_list))

                    # Construct a list of errors
                    error_list = [f'  {m["message"]}' for m in messages if m['kind'] == 'error']
                    if len(error_list) > 0:
                        errors.append(calc_relpath + '\n' + '\n'.join(error_list))

                    # Write out error messages to the .stdout file
                    f.write(f'{calc.directory}/{calc_relpath}\n' + '\n'.join([f'  {m["kind"]}: {m["message"]}' for m
                                                                              in messages]) + '\n\n')

        # Warn for warnings:
        if len(warnings) > 0:
            message = 'Minor disagreements with benchmark detected\n' + '\n'.join(warnings)
            if len(warnings) == 1:
                message = message.replace('disagreements', 'disagreement')
            utils.warn(message)

        # Raise errors for errors:
        if len(errors) > 0:
            message = 'Major disagreements with benchmark detected\n' + '\n'.join(errors)
            if len(errors) == 1:
                message = message.replace('disagreements', 'disagreement')
            raise ValueError(message)


@pytest.fixture
def stumble(monkeypatch):
    '''
    Deliberately crash the code after every single calculation and attempt to restart
    '''

    class DeliberateCalculationFailed(CalculationFailed):
        '''
        An error speciflcally for when we deliberately crash a calculation
        '''

    stumble_message = 'Deliberately crashing for testing purposes'

    class StumblingWorkflow(object):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.calc_counter = 1
            self.catch_stumbles = True

        def run_calculator_single(self, *args, **kwargs):
            if len(self.calculations) == self.calc_counter:
                self.print(stumble_message)
                raise DeliberateCalculationFailed(stumble_message)
            else:
                super().run_calculator_single(*args, **kwargs)

        def run(self, *args, **kwargs):
            if self.catch_stumbles:
                # Create a copy of the state of the workflow before we start running it
                dct_before_running = {k: copy.deepcopy(v) for k, v in self.__dict__.items() if k not in
                                      ['benchmark', 'calc_counter']}
                self.print(f'ATTEMPT {self.calc_counter}', style='heading')
                try:
                    super().run(*args, **kwargs)
                except DeliberateCalculationFailed:
                    # Restore the workflow to the state it was in before we ran any calculations and attempt to restart
                    # from where we left off
                    self.__dict__.update(**dct_before_running)
                    self.parameters.from_scratch = False
                    self.calc_counter += 1
                    self.print_stumble = True
                    self.run(*args, **kwargs)
            else:
                super().run(*args, **kwargs)

        def run_subworkflow(self, workflow, *args, **kwargs):
            # Prevent subworkflows from catching stumbles
            workflow.catch_stumbles = False
            workflow.calc_counter = self.calc_counter
            super().run_subworkflow(workflow, *args, **kwargs)

    from koopmans.workflows import WannierizeWorkflow

    class StumblingWannierizeWorkflow(StumblingWorkflow, WannierizeWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', StumblingWannierizeWorkflow)

    from koopmans.workflows import FoldToSupercellWorkflow

    class StumblingFoldToSupercellWorkflow(StumblingWorkflow, FoldToSupercellWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.FoldToSupercellWorkflow', StumblingFoldToSupercellWorkflow)

    from koopmans.workflows import KoopmansDSCFWorkflow

    class StumblingKoopmansDSCFWorkflow(StumblingWorkflow, KoopmansDSCFWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', StumblingKoopmansDSCFWorkflow)

    from koopmans.workflows import DFTCPWorkflow

    class StumblingDFTCPWorkflow(StumblingWorkflow, DFTCPWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DFTCPWorkflow', StumblingDFTCPWorkflow)

    from koopmans.workflows import DFTPWWorkflow

    class StumblingDFTPWWorkflow(StumblingWorkflow, DFTPWWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DFTPWWorkflow', StumblingDFTPWWorkflow)

    from koopmans.workflows import DeltaSCFWorkflow

    class StumblingDeltaSCFWorkflow(StumblingWorkflow, DeltaSCFWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DeltaSCFWorkflow', StumblingDeltaSCFWorkflow)

    from koopmans.workflows import KoopmansDFPTWorkflow

    class StumblingKoopmansDFPTWorkflow(StumblingWorkflow, KoopmansDFPTWorkflow):
        def plot_bandstructure(self):
            # We don't want the mock test to generate files (nor does it have access to the band_structure result)
            return

    monkeypatch.setattr('koopmans.workflows.KoopmansDFPTWorkflow', StumblingKoopmansDFPTWorkflow)


@pytest.fixture
def mock_quantum_espresso(monkeypatch, pytestconfig):
    '''
    Replace all calls to pw.x, kcp.x, etc with pretend calculations that check the input is correct and then fetch
    the results from a lookup table rather than running any actual calculations.
    '''

    def tests_directory() -> Path:
        directory = pytestconfig.rootpath
        if directory.parts[-1] != 'tests':
            directory /= 'tests'
        return directory.resolve()

    def relative_directory(path) -> Path:
        return path.relative_to(tests_directory())

    def generic_mock_calculate(calc, from_scratch=True):

        # Write the input file
        calc.write_input(calc.atoms)

        # Load the input file
        calc.read_input()

        # Check that we are using the correct settings
        input_file_name = f'{relative_directory(calc.directory)}/{calc.prefix}{calc.ext_in}'
        for key in set(list(calc.benchmark['parameters'].keys()) + list(calc.parameters.keys())):
            # Don't check starting magnetization because ASE doesn't parse it
            if key.startswith('starting_magnetization'):
                continue

            if key == 'ibrav':
                continue

            assert key in calc.benchmark['parameters'].keys(), f'{key} in {input_file_name} not found in benchmark'
            ref_val = calc.benchmark['parameters'][key]

            assert key in calc.parameters.keys(), f'{key} missing from {input_file_name}'
            val = calc.parameters[key]

            if key in calc.parameters.are_paths:
                # Compare the path relative to the location of the input file (mirroring behaviour of
                # construct_benchmark.py)
                val = Path(os.path.relpath(val, calc.directory))

            if isinstance(val, np.ndarray):
                val = val.tolist()

            if isinstance(ref_val, np.ndarray):
                ref_val = ref_val.tolist()

            if isinstance(val, BandPath):
                assert val.path == ref_val.path, f'{key}.path = {val.path} (test) != {ref_val.path} (benchmark)'
                assert len(val.kpts) == len(ref_val.kpts), \
                    f'Mismatch between len({key}) = {len(val.kpts)} (test) != {len(ref_val.kpts)} (benchmark)'
            else:
                assert val == ref_val, f'{key} = {val} (test) != {ref_val} (benchmark)'

        # Copy over the results
        calc.results = calc.benchmark['results']

        main_output_file = Path(f'{calc.directory}/{calc.prefix}{calc.ext_out}')
        # Create dummy output files
        if main_output_file.is_file():
            # If this input file exists already, this means that in a real calculation it is skipped with from_scratch,
            # so we won't generate any input files
            pass
        else:
            for path in calc.output_files + [main_output_file]:
                directory = path.parent
                filename = path.name

                with utils.chdir(directory):
                    # Create this (nested) directory if it does not exist and...
                    with open(filename, 'w') as f:
                        # ... write the file
                        json.dump({'filename': filename, 'written_to': str(relative_directory(path)),
                                   'written_by': f'{input_file_name}'}, f)

    class MockCalc(object):
        def _ase_calculate(self):
            # Monkeypatching _ase_calculate to replace it with a pretend calculation
            generic_mock_calculate(self)

        def is_complete(self):
            return (self.directory / f'{self.prefix}{self.ext_out}').is_file()

        def check_code_is_installed(self):
            # Don't check if the code is installed
            return

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
            if self.parameters.empty_states_nbnd > 0:
                files += [fname for ispin in range(1, self.parameters.nspin + 1)
                          for fname in [f'evc0_empty{ispin}.dat', f'evc_empty{ispin}.dat']]
            return [Path(f) for f in files]

        @property
        def output_files(self) -> List[Path]:
            files = self.__files
            if self.parameters.print_wfc_anion:
                files.append(Path('evcfixed_empty.dat'))
            return [self.parameters.outdir
                    / Path(f'{self.parameters.prefix}_{self.parameters.ndw}.save/K00001/{fname}') for fname in files]

        @property
        def input_files(self) -> List[Path]:
            files = self.__files
            if self.parameters.restart_from_wannier_pwscf:
                files.append(Path('evc_occupied.dat'))
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
            if self.parameters.wan_mode == 'wannier2odd':
                if self.parameters.split_evc_file:
                    files = [
                        f'{self.directory}/{fname}' for fname in ['evcw1.dat', 'evcw2.dat']]
                else:
                    files = [f'{self.directory}/{fname}' for fname in [
                        'evcw.dat', 'charge-density-x.dat']]
            else:
                files = [
                    f'{self.directory}/{self.parameters.seedname}{suffix}' for suffix in ['.eig', '.mmn', '.amn']]
            return [Path(f) for f in files]

        @property
        def input_files(self) -> List[Path]:
            i_kpoints = range(1, np.prod(self.parameters.kpts) + 1)
            files = [f'{self.parameters.outdir}/{self.parameters.prefix}.save/wfc{i}.dat' for i in i_kpoints]
            files += [f'{self.parameters.outdir}/{self.parameters.prefix}.save/{f}'
                      for f in ['data-file-schema.xml', 'charge-density.dat']]
            files.append(f'{self.directory}/{self.parameters.seedname}.nnkp')
            if self.parameters.wan_mode == 'wannier2odd':
                files.append(f'{self.directory}/{self.parameters.seedname}.chk')
            return [Path(f) for f in files]

    class MockUnfoldAndInterpolateCalculator(MockCalc, UnfoldAndInterpolateCalculator):
        # For the UI calculator, _calculate() plays the role of _ase_calculate()
        def _calculate(self):
            generic_mock_calculate(self)

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

    monkeypatch.setattr('koopmans.calculators.KoopmansCPCalculator', MockKoopmansCPCalculator)
    monkeypatch.setattr('koopmans.calculators.PWCalculator', MockPWCalculator)
    monkeypatch.setattr('koopmans.calculators.EnvironCalculator', MockEnvironCalculator)
    monkeypatch.setattr('koopmans.calculators.Wannier90Calculator', MockWannier90Calculator)
    monkeypatch.setattr('koopmans.calculators.PW2WannierCalculator', MockPW2WannierCalculator)
    monkeypatch.setattr('koopmans.calculators.UnfoldAndInterpolateCalculator', MockUnfoldAndInterpolateCalculator)
    monkeypatch.setattr('koopmans.calculators.Wann2KCCalculator', MockWann2KCCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansScreenCalculator', MockKoopmansScreenCalculator)
    monkeypatch.setattr('koopmans.calculators.KoopmansHamCalculator', MockKoopmansHamCalculator)

    # MockWorkflow class only intended to be used to generate other MockWorkflow subclasses with multiple inheritance
    class MockWorkflow(object):
        def run_calculator_single(self, qe_calc):
            # Monkeypatching workflows to provide individual calculations with their benchmarks and getting them
            # to check if the required input files exist and make sense

            # Check we have a benchmark entry for this calculation, and connect it to the calculator
            qe_calc_seed = str(relative_directory(qe_calc.directory / qe_calc.prefix))
            assert qe_calc_seed in self.benchmark, \
                f'Could not find an entry for {qe_calc_seed} in tests/benchmarks.json'
            qe_calc.benchmark = self.benchmark[qe_calc_seed]

            # Check required input files exist and come from where we expect
            if construct_exceptions:
                input_file_exceptions = {}
            else:
                input_file_exceptions = qe_calc.benchmark.get('input files', {})

            # If this calculator is a pw2wannier object, it need to know how many kpoints there are (usually done via
            # the contents of .nnkp)
            if isinstance(qe_calc, PW2WannierCalculator):
                qe_calc.parameters.kpts = self.kgrid

            # We only need to check input files for calculations...
            # a) not starting from scratch, and
            # b) not being skipped (i.e. self.from_scratch is True)
            if qe_calc.parameters.get('restart_mode', 'restart') != 'from_scratch' \
                    and self.parameters.from_scratch:
                for input_file in qe_calc.input_files:

                    # Check the input file exists
                    assert os.path.isfile(input_file), f'Could not find the required input file {input_file}'

                    # Load the contents of the dummy input file (it is a json file)
                    with open(input_file, 'r') as fd:
                        input_file_info = json.load(fd)
                    del input_file_info['filename']

                    input_file = relative_directory(input_file)

                    # Check this file has come from where we expect
                    if str(input_file) in input_file_exceptions:
                        if not construct_exceptions:
                            # These files are overwritten during the workflow, so we must check these against a lookup
                            # table
                            for k, v in input_file_info.items():
                                assert v == input_file_exceptions[str(input_file)][k]

                    else:
                        if not construct_exceptions:
                            # These files are not overwritten during the workflow
                            # Check it was written by the most recent calculation that wrote to this location
                            match_found = False
                            for c in self.calculations[::-1]:
                                if input_file in [relative_directory(f) for f in c.output_files]:
                                    c_input_file = relative_directory(Path.cwd() / c.directory / (c.prefix + c.ext_in))
                                    c_output_file = c.directory / (c.prefix + c.ext_out)
                                    c_output_file = c_output_file.resolve()
                                    assert c_output_file.is_file()
                                    with open(c_output_file, 'r') as fd:
                                        c_output_file_info = json.load(fd)

                                    # Check that this file wrote its own output file (if it didn't it was skipped
                                    # and won't have produced any output files, so it is not a valid match)
                                    if c_output_file_info['written_by'] == str(c_input_file):
                                        match_found = True

                                        # Check that this calculator did indeed write this file
                                        assert str(c_input_file) == input_file_info['written_by'], \
                                            f'{input_file} has been mistakenly overwritten'
                                        break

                            assert match_found, f'Could not find {input_file_info["written_by"]} that generated' \
                                f' {input_file}'

                            # Check this file was written to its current location
                            assert str(input_file) == input_file_info['written_to'], \
                                f'{input_file} has been overwritten by {input_file_info["written_to"]}'
                        else:
                            # Populate benchmarks.json with exceptions
                            if not relative_directory(Path(input_file).resolve()) == input_file_info['written_to']:
                                # Add entry to exceptions dict
                                fname = tests_directory() / 'benchmarks.json'
                                assert fname.is_file()
                                with open(fname, 'r') as fd:
                                    exceptions = read_encoded_json(fd)

                                assert qe_calc_seed in exceptions, f'Benchmark for {qe_calc_seed} is missing'
                                if 'input files' not in exceptions[qe_calc_seed]:
                                    exceptions[qe_calc_seed]['input files'] = {}
                                exceptions[qe_calc_seed]['input files'][str(input_file)] = input_file_info

                                with open(fname, 'w') as fd:
                                    write_encoded_json(exceptions, fd)
            super().run_calculator_single(qe_calc)

        def load_old_calculator(self, qe_calc):
            # Load old calculators by looking through the workflow's list of previous calculations
            # (During the test suite, the only scenario where we want to reload an old calculation will arise
            # because we did it earlier in the same workflow)

            # Load the dummy output file
            calc_fname = qe_calc.directory / f'{qe_calc.prefix}{qe_calc.ext_out}'
            with open(calc_fname, 'r') as fd:
                output_file_info = json.load(fd)

            # Find out the calculation that generated the output file
            directory, prefix = output_file_info['written_by'].rsplit('.', 1)[0].rsplit('/', 1)
            directory = tests_directory() / directory
            matches = [c for c in self.calculations if c.directory == directory and c.prefix == prefix]

            # Copy that calculation into the record of all calculations
            if len(matches) == 1:
                # Fetch the results from the match
                qe_calc.results = matches[0].results
                self.calculations.append(qe_calc)
            elif len(matches) == 0:
                raise ValueError(f'Could not find a calculator matching {qe_calc.directory}/{qe_calc.prefix}')
            else:
                raise ValueError(f'Found multiple calculators for {qe_calc.directory}/{qe_calc.prefix}')

            return qe_calc.is_complete()

    from koopmans.workflows import WannierizeWorkflow

    class MockWannierizeWorkflow(MockWorkflow, WannierizeWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.WannierizeWorkflow', MockWannierizeWorkflow)

    from koopmans.workflows import FoldToSupercellWorkflow

    class MockFoldToSupercellWorkflow(MockWorkflow, FoldToSupercellWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.FoldToSupercellWorkflow', MockFoldToSupercellWorkflow)

    from koopmans.workflows import KoopmansDSCFWorkflow

    class MockKoopmansDSCFWorkflow(MockWorkflow, KoopmansDSCFWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.KoopmansDSCFWorkflow', MockKoopmansDSCFWorkflow)

    from koopmans.workflows import DFTCPWorkflow

    class MockDFTCPWorkflow(MockWorkflow, DFTCPWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DFTCPWorkflow', MockDFTCPWorkflow)

    from koopmans.workflows import DFTPWWorkflow

    class MockDFTPWWorkflow(MockWorkflow, DFTPWWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DFTPWWorkflow', MockDFTPWWorkflow)

    from koopmans.workflows import DeltaSCFWorkflow

    class MockDeltaSCFWorkflow(MockWorkflow, DeltaSCFWorkflow):
        pass

    monkeypatch.setattr('koopmans.workflows.DeltaSCFWorkflow', MockDeltaSCFWorkflow)

    from koopmans.workflows import KoopmansDFPTWorkflow

    class MockKoopmansDFPTWorkflow(MockWorkflow, KoopmansDFPTWorkflow):
        def plot_bandstructure(self):
            # We don't want the mock test to generate files (nor does it have access to the band_structure result)
            return

    monkeypatch.setattr('koopmans.workflows.KoopmansDFPTWorkflow', MockKoopmansDFPTWorkflow)

    # Monkeypatch find_executable, since we don't want to actually try and find the program
    def mock_find_executable(program):
        return './' + program

    monkeypatch.setattr('koopmans.utils.find_executable', mock_find_executable)
