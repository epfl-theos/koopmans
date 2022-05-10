from koopmans.calculators import KoopmansCPCalculator
from koopmans.io import write_kwf as write_encoded_json
from ._utils import benchmark_filename


class BenchmarkGenCalc():
    def _ase_calculate(self):
        super()._ase_calculate()

        # Make sure we store all paths as relative paths
        self.parameters.use_relative_paths = True

        # Save the calculator as an encoded json in tests/benchmarks/
        with open(benchmark_filename(self), 'w') as fd:
            write_encoded_json(self, fd)


class BenchGenKoopmansCPCalculator(BenchmarkGenCalc, KoopmansCPCalculator):
    pass
