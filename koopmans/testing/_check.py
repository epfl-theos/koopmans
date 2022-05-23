from abc import ABC, abstractproperty
import json
import numpy as np
import os
from pathlib import Path
from typing import List, Union
from ase.dft.kpoints import BandPath
from ase.spectrum.band_structure import BandStructure
from ase.spectrum.doscollection import GridDOSCollection
from ase.dft.dos import DOS
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator
from koopmans.workflows import WannierizeWorkflow
from koopmans import utils, projections, base_directory
from koopmans.io import read_kwf as read_encoded_json
from ._utils import benchmark_filename

# A hard and a soft tolerance for checking floats
tolerances = {'alpha': (2e-3, 2e-5),
              'homo_energy': (2e-3, 2e-5),
              'lumo_energy': (2e-3, 2e-5),
              'self-hartree': (2e-3, 2e-5),
              'array': (2e-2, 2e-4),
              'default': (1e-4, 1e-7)}


def compare(result, ref_result, result_name):
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
        assert result.shape == ref_result.shape, f'Array shape mismatch between result and benchmark for {result_name}'
        result = result.flatten()
        ref_result = ref_result.flatten()

        # Calculate the absolute difference
        abs_diffs = np.abs(result - ref_result)

        # Calculate the mixed difference
        scale_factor = 0.1 * np.max(np.abs(ref_result))
        if scale_factor == 0.0:
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


class CheckCalc(ABC):
    @abstractproperty
    def results_for_qc() -> List[str]:
        ...

    def calculate(self):
        # Before running the calculation, check the settings are the same

        with utils.chdir(self.directory):
            # By moving into the directory where the calculation was run, we ensure when we read in the settings that
            # paths are interpreted relative to this particular working directory
            with open(benchmark_filename(self), 'r') as fd:
                benchmark = read_encoded_json(fd)

        # Compare the settings
        for key in set(list(self.parameters.keys()) + list(benchmark.parameters.keys())):
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
            messages = []

            # Loop through results that require checking
            for result_name, ref_result in benchmark.results.items():
                # Only inspect results listed in self.results_for_qc
                if result_name not in self.results_for_qc:
                    continue
                assert result_name in self.results, f'Error in {calc_path}: {result_name} is missing'
                result = self.results[result_name]

                # Check the result against the benchmark
                message = compare(result, ref_result, result_name)

                if message is not None:
                    messages.append(message)

            # Construct the name of this calculation
            calcname = (self.directory / self.prefix).relative_to(base_directory / 'tests' / 'tmp')
            calcname = calcname.relative_to(calcname.parts[0])

            # Warn for warnings:
            warnings = [f'  {m["message"]}' for m in messages if m['kind'] == 'warning']
            if len(warnings) > 0:
                message = f'Minor disagreements with benchmark detected for {calcname}\n' + '\n'.join(warnings)
                if len(warnings) == 1:
                    message = message.replace('disagreements', 'disagreement')
                utils.warn(message)

            # Raise errors for errors:
            errors = [f'  {m["message"]}' for m in messages if m['kind'] == 'error']
            if len(errors) > 0:
                message = f'Major disagreements with benchmark detected for {calcname}\n' + '\n'.join(errors)
                if len(errors) == 1:
                    message = message.replace('disagreements', 'disagreement')
                raise ValueError(message)


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
    pass


class CheckPW2WannierCalculator(CheckCalc, PW2WannierCalculator):
    results_for_qc = []
    pass


class CheckWann2KCPCalculator(CheckCalc, Wann2KCPCalculator):
    results_for_qc = []
    pass


class CheckUnfoldAndInterpolateCalculator(CheckCalc, UnfoldAndInterpolateCalculator):
    results_for_qc = ['band structure', 'dos']
    pass


class CheckWann2KCCalculator(CheckCalc, Wann2KCCalculator):
    results_for_qc = []
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
