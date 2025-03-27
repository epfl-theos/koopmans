import itertools
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import numpy as np
from ase_koopmans.calculators.calculator import CalculationFailed
from ase_koopmans.dft.dos import DOS
from ase_koopmans.spectrum.band_structure import BandStructure
from ase_koopmans.spectrum.doscollection import GridDOSCollection

from koopmans import utils
from koopmans.calculators import (Calc, EnvironCalculator,
                                  KoopmansCPCalculator, KoopmansHamCalculator,
                                  KoopmansScreenCalculator, PhCalculator,
                                  ProjwfcCalculator, PW2WannierCalculator,
                                  PWCalculator, Wann2KCCalculator,
                                  Wann2KCPCalculator, Wannier90Calculator)
from koopmans.files import File
from koopmans.io import read_pkl
from koopmans.processes.bin2xml import Bin2XMLProcess
from koopmans.processes.koopmans_cp import (ConvertFilesFromSpin1To2,
                                            ConvertFilesFromSpin2To1)
from koopmans.processes.power_spectrum import (
    ComputePowerSpectrumProcess, ExtractCoefficientsFromXMLProcess)
from koopmans.processes.ui import UnfoldAndInterpolateProcess
from koopmans.processes.wannier import ExtendProcess, MergeProcess

from ._utils import benchmark_filename, metadata_filename

# A hard and a soft tolerance for checking floats
tolerances = {'alphas': (2e-3, 2e-5),
              'eigenenergies': (2e-3, 2e-5),
              'centersandspreads': (2e-2, 2e-4),
              'array': (2e-2, 2e-4),
              'default': (2e-4, 2e-6)}


def compare(result: Any, ref_result: Any, result_name: str) -> Optional[Dict[str, str]]:
    """Compare the result of a calculation to a reference result."""
    # Sanitize the input and fetch the corresponding tolerances
    if isinstance(result, BandStructure):
        result = result.energies
        ref_result = ref_result.energies
        tols = tolerances['array']
    elif isinstance(result, GridDOSCollection):
        result = result._weights
        ref_result = ref_result._weights
        tols = tolerances['array']
    elif isinstance(result, DOS):
        result = result.get_dos()
        ref_result = ref_result.get_dos()
        tols = tolerances['array']
    elif result_name == 'dielectric tensor':
        tols = tolerances['array']
    elif result_name in ['homo_energy', 'lumo_energy', 'eigenvalues', 'ki_eigenvalues_on_grid']:
        tols = tolerances.get('eigenenergies', tolerances['eigenenergies'])
    else:
        tols = tolerances.get(result_name, tolerances['default'])

    if isinstance(ref_result, float):
        diff = result - ref_result
        message = f'{result_name} = {result:.5f} differs from benchmark {ref_result:.5f} by {diff:.2e}'
        if abs(diff) > tols[0]:
            return {'kind': 'error', 'message': message}
        elif abs(diff) > tols[1]:
            return {'kind': 'warning', 'message': message}

    elif isinstance(ref_result, np.ndarray):
        # For arrays, perform a mixed error test. If Delta = |x - x_ref| then in the limit of large Delta,
        # then this reduces to testing relative error, whereas in the limit of small Delta it reduces to
        # testing absolute error. We use 0.1*max(ref_result) as a reference scale factor.

        # Sanitize the input
        if result.shape != ref_result.shape:
            return {'kind': 'error', 'message': f'Array shape mismatch between result and benchmark for {result_name}'}
        result = result.flatten()
        ref_result = ref_result.flatten()

        # Calculate the absolute difference
        abs_diffs = np.abs(result - ref_result)

        # Calculate the mixed difference
        scale_factor = 0.1 * np.max(np.abs(ref_result))
        if scale_factor < 1e-10:
            # The array appears to be empty
            scale_factor = 1.0
        mixed_diffs = abs_diffs / (scale_factor + np.abs(ref_result))

        # Locate the datapoint with the largest mixed difference
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

    elif isinstance(ref_result, list):
        # For lists, the list is first flattened and the elements are compared one by one
        for a, b in zip(utils.flatten(ref_result), utils.flatten(result)):
            diff = b - a
            message = f'{result_name} = {b:.5f} differs from benchmark {a:.5f} by {diff:.2e}'
            if abs(diff) > tols[0]:
                return {'kind': 'error', 'message': message}
            elif abs(diff) > tols[1]:
                return {'kind': 'warning', 'message': message}

    else:
        if result != ref_result:
            message = f'{result_name} = {result} differs from benchmark {ref_result}'
            return {'kind': 'error', 'message': message}

    return None


def generate_messages(self, benchmark: Calc, results_for_qc: List[str]) -> List[Dict[str, str]]:
    """Generate a list of messages that describe the differences between the calculation and the benchmark."""
    messages: List[Dict[str, str]] = []

    # Loop through results that require checking
    for result_name, ref_result in benchmark.results.items():
        # Only inspect results listed in self.results_for_qc
        if result_name not in results_for_qc:
            continue
        assert result_name in self.results, f'Error in {self.name}: {result_name} is missing'
        result = self.results[result_name]

        # Check the result against the benchmark
        message = compare(result, ref_result, result_name)

        if message is not None:
            messages.append(message)

    return messages


def print_messages(self, messages: List[Dict[str, str]]) -> None:
    """Print messages that describe the differences between the results of a calculation and the benchmark."""
    # Warn for warnings:
    warnings = [f'  {m["message"]}' for m in messages if m['kind'] == 'warning']
    if len(warnings) > 0:
        message = f'Minor disagreements with benchmark detected for {self.name}\n' + '\n'.join(warnings)
        if len(warnings) == 1:
            message = message.replace('disagreements', 'disagreement')
        utils.warn(message)

    # Raise errors for errors:
    errors = [f'  {m["message"]}' for m in messages if m['kind'] == 'error']
    if len(errors) > 0:
        message = f'Major disagreements with benchmark detected for {self.name}\n' + '\n'.join(errors)
        if len(errors) == 1:
            message = message.replace('disagreements', 'disagreement')

        raise CalculationFailed(message)


def load_benchmark(self) -> Calc:
    """Load the benchmark for a calculation."""
    with utils.chdir(self.directory):  # type: ignore[attr-defined]
        # By moving into the directory where the calculation was run, we ensure when we read in the settings that
        # paths are interpreted relative to this particular working directory
        benchmark = read_pkl(benchmark_filename(self))
    return benchmark


def centers_spreads_allclose(center0, spread0, center1, spread1, tol):
    """Check if two sets of centers and spreads are close to each other."""
    if abs(spread0 - spread1) > tol:
        return False
    for x in itertools.product([0, 1, -1], repeat=3):
        if np.allclose(center0, center1 + x, rtol=0, atol=tol):
            return True
    return False


def wannier_generate_messages(self, benchmark: Calc) -> List[Dict[str, str]]:
    """Generate helpful messages for the Wannier90 calculator.

    For a Wannier90 calculator, we don't want to check the contents of self.results key-by-key. Instead, we want
    to check if the same wannier functions have been found, accounting for permutations and for periodic images

    """
    messages: List[Dict[str, str]] = []

    for key in ['centers', 'spreads']:
        assert len(self.results[key]) == len(benchmark.results[key]), f'The list of {key} is the wrong length'

    # If an empty list, return immediately
    if len(benchmark.results['centers']) == 0:
        return []

    # Translate wannier centers to the unit cell
    centers = self.atoms.cell.scaled_positions(np.array(self.results['centers'])) % 1
    ref_centers = self.atoms.cell.scaled_positions(np.array(benchmark.results['centers'])) % 1

    for i, (ref_center, ref_spread) in enumerate(zip(ref_centers, benchmark.results['spreads'])):
        match = False
        rough_match = False
        for j, (center, spread) in enumerate(zip(centers, self.results['spreads'])):
            if centers_spreads_allclose(center, spread, ref_center, ref_spread, tolerances['centersandspreads'][0]):
                match = True
                break
            elif centers_spreads_allclose(center, spread, ref_center, ref_spread,
                                          tolerances['centersandspreads'][1]):
                rough_match = True
                match_index = j
                match_spread = ref_spread
                match_center_str = ', '.join([f'{x:.5f}' for x in benchmark.results['centers'][match_index]])

        ref_center_str = ', '.join([f'{x:.5f}' for x in benchmark.results['centers'][i]])

        if match:
            pass
        elif rough_match:
            message = f'Wannier function #{i+1}, with center = {ref_center_str} and spread = {ref_spread:.5f} ' \
                      f'does not precisely match the benchmark Wannier function #{j+1}, with center = ' \
                      f'{match_center_str} and spread = {match_spread:.5f}'
            messages.append({'kind': 'warning', 'message': message})
        else:
            message = f'Wannier function #{i+1}, with center = {ref_center_str} and spread = {ref_spread:.5f} ' \
                'not found'
            messages.append({'kind': 'error', 'message': message})

    return messages


def compare_output(output, bench_output):
    """Recursively compare the contents of any Files or Paths within the output against the benchmark."""
    if isinstance(output, (File, Path)):
        # Compare the file contents
        binary_formats = ['.npy', '.dat', '.xml']
        if isinstance(output, File):
            binary = output.name.suffix in binary_formats
            if binary:
                bench_output_contents = bench_output.read_bytes()
                output_contents = output.read_bytes()
            else:
                bench_output_contents = bench_output.read_text()
                output_contents = output.read_text()
        else:
            mode = 'rb' if output.suffix in binary_formats else 'r'
            with open(output, mode) as fd:
                output_contents = fd.read()
            with open(bench_output, mode) as fd:
                bench_output_contents = fd.read()

        if isinstance(output_contents, np.ndarray):
            assert np.all(output_contents == bench_output_contents), \
                f'Contents of {output} differs from the benchmark'
        else:
            assert output_contents == bench_output_contents, \
                f'Contents of {output} differs from the benchmark'
    elif isinstance(output, dict):
        # Recurse over the dictionary
        for k, v in output.items():
            compare_output(v, bench_output[k])
    elif isinstance(output, list):
        # Recurse over the list
        for v1, v2 in zip(output, bench_output, strict=True):
            compare_output(v1, v2)
    else:
        pass


def patch_calculator(c, monkeypatch, results_for_qc, generate_messages_function):
    """Patch a calculator."""
    unpatched_pre_calculate = c._pre_calculate
    unpatched_post_calculate = c._post_calculate

    def _pre_calculate(self):
        """Check the settings are the same before running the calculation."""
        # Perform the pre_calculate first, as sometimes this function modifies the input parameters
        unpatched_pre_calculate(self)

        # Load the benchmark
        benchmark = load_benchmark(self)

        # Compare the settings
        unique_keys: Set[str] = set(list(self.parameters.keys()) + list(benchmark.parameters.keys()))
        for key in unique_keys:
            if key not in self.parameters:
                raise ValueError(f'Error in {self.prefix}: {key} is in benchmark but not in test')
            elif key not in benchmark.parameters:
                raise ValueError(f'Error in {self.prefix}: {key} is in test but not in benchmark')
            else:
                val = self.parameters[key]
                ref_val = benchmark.parameters[key]

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

    def _post_calculate(self, generate_messages_function=generate_messages):
        """Patch the post_calculate method to check the results against the benchmark."""
        # Perform the post_calculate first, as sometimes this function adds extra entries to self.results
        unpatched_post_calculate(self)

        # Check the results
        if self.skip_qc:
            # For some calculations (e.g. dummy calculations) we don't care about the results and actively don't
            # want to compare them against a benchmark. For these calculations, we set self.skip_qc to False inside
            # the corresponding workflow
            pass
        else:
            benchmark = load_benchmark(self)
            messages = generate_messages_function(self, benchmark, results_for_qc)
            print_messages(self, messages)

        # Check the expected files were produced
        with open(metadata_filename(self), 'r') as fd:
            metadata: Dict[str, str] = json.load(fd)

        for output_file in metadata['output_files']:
            assert (self.directory
                    / output_file).exists(), f'Error in {self.name}: {output_file} was not generated'

    monkeypatch.setattr(c, '_pre_calculate', _pre_calculate)
    monkeypatch.setattr(c, '_post_calculate', _post_calculate)


def patch_process(p, monkeypatch):
    """Patch a process."""
    unpatched_run = p._run

    def _run(self):
        """Run the process and additionally check the results against the benchmark."""
        # Load the benchmark
        bench_process = read_pkl(benchmark_filename(self), base_directory=self.base_directory)

        # Compare the inputs
        assert bench_process.inputs == self.inputs

        unpatched_run(self)

        # Compare the outputs
        assert bench_process.outputs == self.outputs

        # If any outputs are files, compare the file contents
        for k, output in self.outputs:
            bench_output = getattr(bench_process.outputs, k)

            # Recurse over the output and compare against the benchmark (it might be a dictionary or a list)
            compare_output(output, bench_output)

    monkeypatch.setattr(p, '_run', _run)


def monkeypatch_check(monkeypatch):
    """Monkeypatch all the calculators and processes to check their results against the benchmarks."""
    # Calculators
    results_for_qc = {KoopmansCPCalculator: ['energy', 'homo_energy', 'lumo_energy'],
                      PhCalculator: ['dielectric tensor'],
                      PWCalculator: ['energy', 'eigenvalues', 'band structure'],
                      EnvironCalculator: ['energy', 'electrostatic embedding'],
                      Wannier90Calculator: ['centers', 'spreads'],
                      PW2WannierCalculator: [],
                      Wann2KCPCalculator: [],
                      Wann2KCCalculator: [],
                      KoopmansScreenCalculator: ['alphas'],
                      KoopmansHamCalculator: ['ki_eigenvalues_on_grid', 'band structure'],
                      ProjwfcCalculator: ['dos']}
    for c, r_for_qc in results_for_qc.items():
        if c == Wannier90Calculator:
            gm_func = wannier_generate_messages
        else:
            gm_func = generate_messages
        patch_calculator(c, monkeypatch, r_for_qc, gm_func)

    # Processes
    for p in [ExtractCoefficientsFromXMLProcess, ComputePowerSpectrumProcess, Bin2XMLProcess, ConvertFilesFromSpin1To2,
              ConvertFilesFromSpin2To1, ExtendProcess, MergeProcess, UnfoldAndInterpolateProcess]:
        patch_process(p, monkeypatch)
