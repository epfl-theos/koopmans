from koopmans.calculators import KoopmansCPCalculator
from koopmans.io import write_kwf as write_encoded_json
from ._utils import benchmark_filename
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator


class BenchmarkGenCalc():
    def calculate(self):
        super().calculate()

        # Make sure we store all paths as relative paths
        tmp, self.parameters.use_relative_paths = self.parameters.use_relative_paths, True

        # Don't store the walltime
        self.results.pop('walltime', None)

        # Save the calculator as an encoded json in the benchmarks directory
        with open(benchmark_filename(self), 'w') as fd:
            write_encoded_json(self, fd)

        # Restore the behaviour of use_relative_paths
        self.parameters.use_relative_paths = tmp


class BenchGenWannier90Calculator(BenchmarkGenCalc, Wannier90Calculator):
    pass


class BenchGenPW2WannierCalculator(BenchmarkGenCalc, PW2WannierCalculator):
    pass


class BenchGenWann2KCPCalculator(BenchmarkGenCalc, Wann2KCPCalculator):
    pass


class BenchGenPWCalculator(BenchmarkGenCalc, PWCalculator):
    pass


class BenchGenKoopmansCPCalculator(BenchmarkGenCalc, KoopmansCPCalculator):
    pass


class BenchGenEnvironCalculator(BenchmarkGenCalc, EnvironCalculator):
    pass


class BenchGenUnfoldAndInterpolateCalculator(BenchmarkGenCalc, UnfoldAndInterpolateCalculator):
    pass


class BenchGenWann2KCCalculator(BenchmarkGenCalc, Wann2KCCalculator):
    pass


class BenchGenKoopmansScreenCalculator(BenchmarkGenCalc, KoopmansScreenCalculator):
    pass


class BenchGenKoopmansHamCalculator(BenchmarkGenCalc, KoopmansHamCalculator):
    pass


class BenchGenProjwfcCalculator(BenchmarkGenCalc, ProjwfcCalculator):
    pass
