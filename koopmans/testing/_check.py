from abc import ABC, abstractmethod
import json
import numpy as np
from pathlib import Path
from typing import List, Optional, Any, Dict, Set
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.dft.dos import DOS
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator, Calc
from koopmans import utils, base_directory
from koopmans.io import read_kwf as read_encoded_json
from ._utils import benchmark_filename, metadata_filename

# A hard and a soft tolerance for checking floats
tolerances = {'alpha': (2e-3, 2e-5),
              'homo_energy': (2e-3, 2e-5),
              'lumo_energy': (2e-3, 2e-5),
              'self-hartree': (2e-3, 2e-5),
              'array': (2e-2, 2e-4),
              'default': (1e-4, 1e-7)}


def compare(result: Any, ref_result: Any, result_name: str) -> Optional[Dict[str, str]]:
    # Compare the calculated result to the reference result

    # Sanitise the input and fetch the corresponding tolerances
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

        # Sanitise the input
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


class CheckCalc:
    results_for_qc: List[str]
    prefix: str
    results: Dict[Any, Any]

    @property
    def _calcname(self) -> Path:
        calcname = (self.directory / self.prefix).relative_to(base_directory  # type: ignore[attr-defined]
                                                              / 'tests' / 'tmp')
        return calcname.relative_to(calcname.parts[0])

    def _check_results(self, benchmark: Calc):
        messages = self._generate_messages(benchmark)
        self._print_messages(messages)

    def _generate_messages(self, benchmark: Calc) -> List[Dict[str, str]]:
        messages: List[Dict[str, str]] = []

        # Loop through results that require checking
        for result_name, ref_result in benchmark.results.items():
            # Only inspect results listed in self.results_for_qc
            if result_name not in self.results_for_qc:
                continue
            assert result_name in self.results, f'Error in {self._calcname}: {result_name} is missing'
            result = self.results[result_name]

            # Check the result against the benchmark
            message = compare(result, ref_result, result_name)

            if message is not None:
                messages.append(message)
        return messages

    def _print_messages(self, messages: List[Dict[str, str]]) -> None:
        # Warn for warnings:
        warnings = [f'  {m["message"]}' for m in messages if m['kind'] == 'warning']
        if len(warnings) > 0:
            message = f'Minor disagreements with benchmark detected for {self._calcname}\n' + '\n'.join(warnings)
            if len(warnings) == 1:
                message = message.replace('disagreements', 'disagreement')
            utils.warn(message)

        # Raise errors for errors:
        errors = [f'  {m["message"]}' for m in messages if m['kind'] == 'error']
        if len(errors) > 0:
            message = f'Major disagreements with benchmark detected for {self._calcname}\n' + '\n'.join(errors)
            if len(errors) == 1:
                message = message.replace('disagreements', 'disagreement')
            raise ValueError(message)

    def calculate(self):
        # Before running the calculation, check the settings are the same

        with utils.chdir(self.directory):
            # By moving into the directory where the calculation was run, we ensure when we read in the settings that
            # paths are interpreted relative to this particular working directory
            with open(benchmark_filename(self), 'r') as fd:
                benchmark = read_encoded_json(fd)

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

        # Run the calculation
        super().calculate()

        # Check the results
        if self.skip_qc:
            # For some calculations (e.g. dummy calculations) we don't care about the results and actively don't
            # want to compare them against a benchmark. For these calculations, we set self.skip_qc to False inside
            # the corresponding workflow
            pass
        else:
            self._check_results(benchmark)

        # Check the expected files were produced
        with open(metadata_filename(self), 'r') as fd:
            metadata: Dict[str, str] = json.load(fd)

        for output_file in metadata['output_files']:
            assert (self.directory
                    / output_file).exists(), f'Error in {self._calcname}: {output_file} was not generated'


class CheckKoopmansCPCalculator(CheckCalc, KoopmansCPCalculator):
    results_for_qc = ['energy', 'homo_energy', 'lumo_energy']
    pass


class CheckPWCalculator(CheckCalc, PWCalculator):
    results_for_qc = ['energy']
    pass


class CheckEnvironCalculator(CheckCalc, EnvironCalculator):
    results_for_qc = ['energy', 'electrostatic embedding']
    pass


class CheckWannier90Calculator(CheckCalc, Wannier90Calculator):
    results_for_qc = ['centers', 'spreads']

    def _generate_messages(self, benchmark: Calc) -> List[Dict[str, str]]:
        # For a Wannier90 calculator, we don't want to check the contents of self.results key-by-key. Instead, we want
        # to check if the same wannier functions have been found, accounting for permutations and for periodic images

        messages: List[Dict[str, str]] = []

        for key in ['centers', 'spreads']:
            assert len(self.results[key]) == len(benchmark.results[key]), f'The list of {key} is the wrong length'

        # If an empty list, return immediately
        if len(benchmark.results['centers']) == 0:
            return []

        # Translate wannier centres to the unit cell
        centers = self.atoms.cell.scaled_positions(np.array(self.results['centers'])) % 1
        ref_centers = self.atoms.cell.scaled_positions(np.array(benchmark.results['centers'])) % 1

        for i, (center, spread) in enumerate(zip(centers, self.results['spreads'])):
            result = center + [spread]
            match = False
            rough_match = False
            for j, (ref_center, ref_spread) in enumerate(zip(ref_centers, benchmark.results['spreads'])):
                ref_result = ref_center + [ref_spread]
                if np.allclose(result, ref_result, tolerances['default'][0]):
                    match = True
                    break
                elif np.allclose(result, ref_result, tolerances['default'][1]):
                    rough_match = True
                    match_index = j
                    match_spread = ref_spread
                    match_center_str = ', '.join([f'{x:.5f}' for x in benchmark.results['centers'][match_index]])

            center_str = ', '.join([f'{x:.5f}' for x in self.results['centers'][i]])

            if match:
                pass
            elif rough_match:
                message = f'Wannier function #{i+1}, with center = {center_str} and spread = {spread:.5f} does not precisely match' \
                          f'the benchmark Wannier function #{j+1}, with center = {match_center_str} and spread = {match_spread:.5f}'
                messages.append({'kind': 'warning', 'message': message})
            else:
                message = f'Wannier function #{i+1}, with center = {center_str} and spread = {spread:.5f} not found'
                messages.append({'kind': 'error', 'message': message})

        return messages


class CheckPW2WannierCalculator(CheckCalc, PW2WannierCalculator):
    results_for_qc: List[str] = []
    pass


class CheckWann2KCPCalculator(CheckCalc, Wann2KCPCalculator):
    results_for_qc: List[str] = []
    pass


class CheckUnfoldAndInterpolateCalculator(CheckCalc, UnfoldAndInterpolateCalculator):
    results_for_qc = ['band structure', 'dos']
    pass


class CheckWann2KCCalculator(CheckCalc, Wann2KCCalculator):
    results_for_qc: List[str] = []
    pass


class CheckKoopmansScreenCalculator(CheckCalc, KoopmansScreenCalculator):
    results_for_qc = ['alphas']
    pass


class CheckKoopmansHamCalculator(CheckCalc, KoopmansHamCalculator):
    results_for_qc = ['ki_eigenvalues_on_grid', 'band structure']
    pass


class CheckProjwfcCalculator(CheckCalc, ProjwfcCalculator):
    results_for_qc = ['dos']
    pass
