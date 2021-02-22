'''
Functionality for testing koopmans code with pytest

Written by Edward Linscott, Dec 2020
'''

import os
import json
import numpy as np
import pytest
from koopmans.calculators.wannier90 import W90_calc
from koopmans.calculators.pw2wannier import PW2Wannier_calc
from koopmans.calculators.pw import PW_calc
from koopmans.calculators.kcp import KCP_calc
from koopmans.calculators.environ import Environ_calc
from koopmans.calculators.ui import UI_calc
from koopmans.io import read_json
from koopmans import utils
from ase import io as ase_io
from ase.dft.kpoints import BandPath


# A hard and a soft tolerance for checking floats
default_tolerances = {'alpha': (2e-3, 2e-5),
                      'homo_energy': (2e-3, 2e-5),
                      'lumo_energy': (2e-3, 2e-5),
                      'self-hartree': (2e-3, 2e-5),
                      'default': (1e-4, 1e-7)}


# Set this to true to construct a benchmark file
construct_exceptions = False


class WorkflowTest:

    def __init__(self, json_input, capsys, mock=False, tolerances=default_tolerances):

        self.directory, self.json = os.path.abspath(json_input).rsplit('/', 1)
        self.tests_directory, self.subdirectory = self.directory.rsplit('/', 1)
        self.mock = mock
        self.capsys = capsys
        self.tolerances = tolerances

        # Fetch the benchmarks
        with open('tests/benchmarks.json', 'r') as f:
            benchmarks = json.load(f)

        self.benchmark = {k: v for k, v in benchmarks.items() if self.subdirectory in k}

    def run(self):
        # Move into the test directory
        with utils.chdir(self.directory):

            # Load the input file
            workflow = read_json(self.json)

            # Give it access to the benchmark
            workflow.benchmark = self.benchmark

            # Make sure QC printing is turned on
            workflow.print_qc = True

            # Run the workflow
            workflow.run()

            # Write the output
            out, _ = self.capsys.readouterr()
            stdout = self.json.replace('.json', '.stdout')
            with open(stdout, 'w') as f:
                f.write(out)

            if not self.mock:
                # Print out the contents of the QE i/o files
                utils.system_call(f'cat_workflow_io.sh >> {stdout}')

                # Check numerical results against the benchmarks
                self.check_qc_results(workflow)

    def check_qc_results(self, workflow):

        # Loop through the stored QC results of all the calculations
        log = {}
        for calc in workflow.all_calcs:

            # Skip checking n-1 calculations since these are known to be unreliable
            if 'n-1' in calc.name:
                continue

            # Prepare the log dict entry for this calc
            calc_relpath = os.path.relpath(f'{calc.directory}/{calc.name}', self.tests_directory)
            log[calc_relpath] = []

            # If there has been no monkey-patching, calc will not have a benchmark attribute yet, so add that
            assert calc_relpath in workflow.benchmark, f'Could not find an entry for {calc_relpath} in benchmarks.json'
            calc.benchmark = workflow.benchmark[calc_relpath]

            # Loop through results that require checking
            for result_name, result in calc.qc_results.items():

                # Fetch the corresponding benchmark
                if result_name.startswith('alpha'):
                    # alphas are not stored in results but in a separate field
                    assert 'alphas' in calc.benchmark
                    i_orb = int(result_name.split('(')[-1][:-1])
                    ref_result = calc.benchmark['alphas'][i_orb - 1]
                    tols = self.tolerances['alpha']
                elif result_name.startswith('orbital_self_Hartree'):
                    assert 'orbital_data' in calc.benchmark['results']
                    [i_orb, i_spin] = [int(s.split('=')[1]) for s in result_name.split('(')[-1][:-1].split(',')]
                    ref_result = calc.benchmark['results']['orbital_data']['self-Hartree'][i_spin - 1][i_orb - 1]
                    tols = self.tolerances['self-hartree']
                else:
                    assert result_name in calc.benchmark['results']
                    ref_result = calc.benchmark['results'][result_name]
                    tols = self.tolerances.get(result_name, self.tolerances['default'])

                # Compare the calculated result to the reference result
                if isinstance(ref_result, float):
                    message = f'{result_name} = {result:.5f} differs from benchmark {ref_result:.5f} by {result - ref_result:.2e}'
                    if abs(result - ref_result) > tols[0]:
                        log[calc_relpath].append({'kind': 'error', 'message': message})
                    elif abs(result - ref_result) > tols[1]:
                        log[calc_relpath].append({'kind': 'warning', 'message': message})
                else:
                    if result != ref_result:
                        message = f'{result_name} = {result} differs from benchmark {ref_result}'
                        log[calc_relpath].append({'kind': 'error', 'message': message})

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
                    f.write(f'{calc.directory}/{calc_relpath}\n' + '\n'.join([f'  {m["kind"]}: {m["message"]}' for m in messages]) + '\n\n')

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
def mock_quantum_espresso(monkeypatch, pytestconfig):
    '''
    Replace all calls to pw.x, kcp.x, etc with pretend calculations that check the input is correct and then fetch
    the results from a lookup table rather than running any actual calculations.
    '''

    # Prevent run_calculator_single from attempting to read output files
    import koopmans.workflows.generic as wf_generic
    wf_generic.skip_loading_outputs = True

    def tests_directory():
        directory = pytestconfig.rootpath
        if directory.parts[-1] != 'tests':
            directory /= 'tests'
        return os.path.abspath(directory)

    def relative_directory(path):
        return os.path.relpath(path, tests_directory())

    def generic_mock_calculate(calc, from_scratch=True):

        # Write the input file
        calc.write_input_file()

        # Load the input file
        calc.load_input_file()

        # Check that we are using the correct settings
        input_file_name = f'{relative_directory(calc.directory)}/{calc.name}{calc.ext_in}'
        for key in set(list(calc.benchmark['settings'].keys()) + list(calc.settings.keys())):
            # Don't check starting magnetization because ASE doesn't parse it
            if key.startswith('starting_magnetization'):
                continue

            assert key in calc.benchmark['settings'].keys(), f'{key} in {input_file_name} not found in benchmark'
            ref_val = calc.benchmark['settings'][key]

            assert key in calc.settings.keys(), f'{key} missing from {input_file_name}'
            val = calc.settings[key]

            if key in calc._settings_that_are_paths:
                # Compare the path relative to the location of the input file (mirroring behaviour of
                # construct_benchmark.py)
                val = os.path.relpath(val, calc.directory)
            assert val == ref_val, f'{key} = {val} (test) != {ref_val} (benchmark)'

        # Copy over the results
        calc.results = calc.benchmark['results']

        main_output_file = f'{calc.directory}/{calc.name}{calc.ext_out}'
        # Create dummy output files
        if os.path.isfile(main_output_file):
            # If this input file exists already, this means that in a real calculation it is skipped with from_scratch,
            # so we won't generate any input files
            pass
        else:
            for path in calc.output_files + [main_output_file]:
                directory, filename = path.rsplit('/', 1)
                with utils.chdir(directory):
                    # Create this (nested) directory if it does not exist and...
                    with open(filename, 'w') as f:
                        # ... write the file
                        json.dump({'filename': filename, 'written_to': relative_directory(path),
                                   'written_by': f'{input_file_name}'}, f)

    class mock_calc(object):
        def _ase_calculate(self):
            # Monkeypatching _ase_calculate to replace it with a pretend calculation
            generic_mock_calculate(self)

        def is_complete(self):
            return os.path.isfile(f'{self.directory}/{self.name}{self.ext_out}')

        def check_code_is_installed(self):
            return True

    class mock_KCP_calc(mock_calc, KCP_calc):
        # Monkeypatched KCP_calc class which never actually calls kcp.x

        @property
        def __files(self):
            files = [fname for ispin in range(1, self.nspin + 1) for fname in [f'evc0{ispin}.dat', f'evc{ispin}.dat',
                                                                               f'evcm{ispin}.dat',
                                                                               f'hamiltonian{ispin}.xml',
                                                                               f'eigenval{ispin}.xml',
                                                                               f'lambda0{ispin}.dat',
                                                                               f'lambdam{ispin}.dat']]
            if self.empty_states_nbnd > 0:
                files += [fname for ispin in range(1, self.nspin + 1) for fname in [f'evc0_empty{ispin}.dat',
                                                                                    f'evc_empty{ispin}.dat']]
            return files

        @property
        def output_files(self):
            files = self.__files
            if self.print_wfc_anion:
                files.append('evcfixed_empty.dat')
            return [f'{self.outdir}/{self.prefix}_{self.ndw}.save/K00001/{fname}' for fname in files]

        @property
        def input_files(self):
            files = self.__files
            if self.restart_from_wannier_pwscf:
                files.append('evc_occupied.dat')
            return [f'{self.outdir}/{self.prefix}_{self.ndr}.save/K00001/{fname}' for fname in files]

    class mock_Environ_calc(mock_calc, Environ_calc):
        # Monkeypatched Environ_calc class which never actually calls pw.x

        @property
        def output_files(self):
            files = []
            if 'kpts' in self.calc.parameters:
                assert self.calc.parameters['kpts'] == [
                    1, 1, 1], 'Have not implemented dummy environ calculations for kpts != [1, 1, 1]'

                i_kpoints = range(1)
                files += [f'{self.outdir}/{self.prefix}.save/wfc{i}.dat' for i in i_kpoints]
                if self.nosym:
                    files += [f'{self.outdir}/wfc{i}.dat' for i in i_kpoints]
            files += [f'{self.outdir}/{self.prefix}.xml']
            return files

        @property
        def input_files(self):
            # Not yet implemented
            return []

    class mock_PW_calc(mock_calc, PW_calc):
        # Monkeypatched PW_calc class which never actually calls pw.x

        @property
        def output_files(self):
            files = []
            if 'kpts' in self.calc.parameters:
                if self.nosym:
                    n_kpoints = np.prod(self.calc.parameters['kpts']) + 1
                else:
                    n_kpoints = len(BandPath(self.calc.parameters['kpts'], self.calc.atoms.cell).kpts)
                i_kpoints = range(1, n_kpoints + 1)

                files += [f'{self.outdir}/{self.prefix}.save/wfc{i}.dat' for i in i_kpoints]
            files += [f'{self.outdir}/{self.prefix}.xml']
            files += [f'{self.outdir}/{self.prefix}.save/{f}' for f in ['data-file-schema.xml', 'charge-density.dat']]
            return files

        @property
        def input_files(self):
            files = []
            if self.calculation == 'nscf':
                files += [f'{self.outdir}/{self.prefix}.save/{f}' for f in ['data-file-schema.xml', 'charge-density.dat']]
            return files

    class mock_W90_calc(mock_calc, W90_calc):
        # Monkeypatched W90_calc class which never actually calls wannier90.x

        @property
        def output_files(self):
            if '-pp' in self.preprocessing_flags:
                files = [f'{self.directory}/{self.name}.nnkp']
            else:
                files = [f'{self.directory}/{self.name}{suffix}' for suffix in [
                    '.chk', '_wsvec.dat', '_hr.dat']]
            return files

        @property
        def input_files(self):
            if '-pp' in self.preprocessing_flags:
                files = []
            else:
                files = [f'{self.directory}/{self.name}{suffix}' for suffix in ['.eig', '.mmn', '.amn']]
            return files

    class mock_PW2Wannier_calc(mock_calc, PW2Wannier_calc):
        # Monkeypatched PW_calc class which never actually calls pw.x

        @property
        def output_files(self):
            if self.wan_mode == 'wannier2odd':
                if self.split_evc_file:
                    files = [
                        f'{self.directory}/{fname}' for fname in ['evcw1.dat', 'evcw2.dat']]
                else:
                    files = [f'{self.directory}/{fname}' for fname in [
                        'evcw.dat', 'charge-density-x.dat']]
            else:
                files = [
                    f'{self.directory}/{self.seedname}{suffix}' for suffix in ['.eig', '.mmn', '.amn']]
            return files

        @property
        def input_files(self):
            i_kpoints = range(1, np.prod(self.calc.parameters['kpts']) + 1)
            files = [f'{self.outdir}/{self.prefix}.save/wfc{i}.dat' for i in i_kpoints]
            files += [f'{self.outdir}/{self.prefix}.save/{f}' for f in ['data-file-schema.xml', 'charge-density.dat']]
            files.append(f'{self.directory}/{self.seedname}.nnkp')
            if self.wan_mode == 'wannier2odd':
                files.append(f'{self.directory}/{self.seedname}.chk')
            return files

    class mock_UI_calc(mock_calc, UI_calc):
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

    monkeypatch.setattr('koopmans.calculators.kcp.KCP_calc', mock_KCP_calc)
    monkeypatch.setattr('koopmans.calculators.pw.PW_calc', mock_PW_calc)
    monkeypatch.setattr('koopmans.calculators.environ.Environ_calc', mock_Environ_calc)
    monkeypatch.setattr('koopmans.calculators.wannier90.W90_calc', mock_W90_calc)
    monkeypatch.setattr('koopmans.calculators.pw2wannier.PW2Wannier_calc', mock_PW2Wannier_calc)
    monkeypatch.setattr('koopmans.calculators.ui.UI_calc', mock_UI_calc)

    # Monkeypatching workflows to provide individual calculations with their benchmarks and getting them
    # to check if the required input files exist and make sense

    def generic_mock_calculator_single(workflow, qe_calc):
        # Check we have a benchmark entry for this calculation, and connect it to the calculator
        qe_calc_seed = relative_directory(qe_calc.directory + '/' + qe_calc.name)
        assert qe_calc_seed in workflow.benchmark, f'Could not find an entry for {qe_calc_seed} in tests/benchmarks.json'
        qe_calc.benchmark = workflow.benchmark[qe_calc_seed]

        # Check required input files exist and come from where we expect
        if construct_exceptions:
            input_file_exceptions = {}
        else:
            input_file_exceptions = qe_calc.benchmark.get('input files', {})

        # If this calculator is a pw2wannier object, it need to know how many kpoints there are (usually done via the contents of .nnkp)
        if isinstance(qe_calc, PW2Wannier_calc):
            recent_pw_calc = [c for c in workflow.all_calcs if isinstance(c, PW_calc)][-1]
            qe_calc.calc.parameters['kpts'] = recent_pw_calc.calc.parameters['kpts']

        # We only need to check input files for calculations...
        # a) not starting from scratch, and
        # b) not being skipped (i.e. workflow.from_scratch is True)
        if getattr(qe_calc, 'restart_mode', 'restart') != 'from_scratch' and workflow.from_scratch:
            for input_file in qe_calc.input_files:

                # Check the input file exists
                assert os.path.isfile(input_file), f'Could not find the required input file {input_file}'

                # Load the contents of the dummy input file (it is a json file)
                with open(input_file, 'r') as fd:
                    input_file_info = json.load(fd)
                del input_file_info['filename']

                input_file = relative_directory(input_file)

                # Check this file has come from where we expect
                if input_file in input_file_exceptions:
                    if not construct_exceptions:
                        # These files are overwritten during the workflow, so we must check these against a lookup
                        # table
                        for k, v in input_file_info.items():
                            assert v == input_file_exceptions[input_file][k]

                else:
                    if not construct_exceptions:
                        # These files are not overwritten during the workflow
                        # Check it was written by the most recent calculation that wrote to this location
                        match_found = False
                        for c in workflow.all_calcs[::-1]:
                            if input_file in [relative_directory(f) for f in c.output_files]:
                                # Check that this file wrote its own output file (if it didn't it was skipped
                                # and won't have produced any output files, so it is not a valid match)
                                c_input_file = relative_directory(os.getcwd() + '/' + c.directory + '/' + c.name + c.ext_in)
                                c_output_file = os.path.abspath(c.directory + '/' + c.name + c.ext_out)
                                assert os.path.isfile(c_output_file)
                                with open(c_output_file, 'r') as fd:
                                    c_output_file_info = json.load(fd)

                                if c_output_file_info['written_by'] == c_input_file:
                                    match_found = True
                                    assert relative_directory(c_input_file) == input_file_info['written_by'], \
                                        f'{input_file} has been mistakenly overwritten'
                                    break

                        assert match_found, f'Could not find {input_file_info["written_by"]} that generated' \
                            f' {input_file}'

                        # Check this file was written to its current location
                        assert relative_directory(input_file) == input_file_info['written_to'], \
                            f'{input_file} has been overwritten by {input_file_info["written_to"]}'
                    else:
                        # Populate benchmarks.json with exceptions
                        if not relative_directory(input_file) == input_file_info['written_to']:
                            # Add entry to exceptions dict
                            fname = tests_directory() + '/benchmarks.json'
                            assert os.path.isfile(fname)
                            with open(fname, 'r') as fd:
                                exceptions = json.load(fd)

                            assert qe_calc_seed in exceptions, f'Benchmark for {qe_calc_seed} is missing'
                            if 'input files' not in exceptions[qe_calc_seed]:
                                exceptions[qe_calc_seed]['input files'] = {}
                            exceptions[qe_calc_seed]['input files'][input_file] = input_file_info
                            with open(fname, 'w') as fd:
                                json.dump(exceptions, fd, indent=2)

    from koopmans.workflows.wf_with_w90 import WannierizeWorkflow

    class MockWannierizeWorkflow(WannierizeWorkflow):

        def run_calculator_single(self, qe_calc):
            generic_mock_calculator_single(self, qe_calc)
            super().run_calculator_single(qe_calc)

    monkeypatch.setattr('koopmans.workflows.wf_with_w90.WannierizeWorkflow', MockWannierizeWorkflow)

    from koopmans.workflows.kc_with_cp import KoopmansWorkflow

    class MockKoopmansWorkflow(KoopmansWorkflow):

        def run_calculator_single(self, qe_calc):
            generic_mock_calculator_single(self, qe_calc)
            super().run_calculator_single(qe_calc)

    monkeypatch.setattr('koopmans.workflows.kc_with_cp.KoopmansWorkflow', MockKoopmansWorkflow)

    from koopmans.workflows.pbe_with_cp import PBEWorkflow

    class MockPBEWorkflow(PBEWorkflow):

        def run_calculator_single(self, qe_calc):
            generic_mock_calculator_single(self, qe_calc)
            super().run_calculator_single(qe_calc)

    monkeypatch.setattr('koopmans.workflows.pbe_with_cp.PBEWorkflow', MockPBEWorkflow)

    from koopmans.workflows.pbe_dscf_with_pw import DeltaSCFWorkflow

    class MockDeltaSCFWorkflow(DeltaSCFWorkflow):

        def run_calculator_single(self, qe_calc):
            generic_mock_calculator_single(self, qe_calc)
            super().run_calculator_single(qe_calc)

    monkeypatch.setattr('koopmans.workflows.pbe_dscf_with_pw.DeltaSCFWorkflow', MockDeltaSCFWorkflow)

    # Monkeypatch find_executable, since we don't want to actually try and find the program
    def mock_find_executable(program):
        return program

    monkeypatch.setattr('koopmans.utils.find_executable', mock_find_executable)
