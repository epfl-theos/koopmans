import os
import json
from typing import List
from koopmans.calculators import KoopmansCPCalculator
from koopmans.io import write_kwf as write_encoded_json
from ._utils import benchmark_filename, metadata_filename, find_subfiles_of_calc
from koopmans.calculators import Wannier90Calculator, PW2WannierCalculator, Wann2KCPCalculator, PWCalculator, \
    KoopmansCPCalculator, EnvironCalculator, UnfoldAndInterpolateCalculator, Wann2KCCalculator, \
    KoopmansScreenCalculator, KoopmansHamCalculator, ProjwfcCalculator


class BenchmarkGenCalc():
    def calculate(self):
        # Before running the calculation, make a list of the files that exist
        files_before = find_subfiles_of_calc(self)

        super().calculate()

        # After running the calculation, make a new list of the files, and then work out which files have been
        # modified by the calculation
        files_after = find_subfiles_of_calc(self)
        modified_files: List[str] = sorted([str(os.path.relpath(x[0], self.directory))
                                           for x in files_after - files_before])
        # Exclude the input file, since we don't want to create a dummy version of this for the mock calc (it is
        # already written in full)
        modified_files.remove(self.prefix + self.ext_in)

        # Exclude .wfc* files because these depend on the parallelism
        if isinstance(self, PWCalculator):
            wfc_prefix = os.path.relpath(self.parameters.outdir, self.directory) + f'/{self.parameters.prefix}.wfc'
            modified_files = [f for f in modified_files if not f.startswith(wfc_prefix)]

        # Write out the calculator itself to file
        # Make sure we store all paths as relative paths
        tmp, self.parameters.use_relative_paths = self.parameters.use_relative_paths, True

        # Don't store the walltime
        self.results.pop('walltime', None)

        # Save the calculator as an encoded json in the benchmarks directory
        with open(benchmark_filename(self), 'w') as fd:
            write_encoded_json(self, fd)

        # Restore the behaviour of use_relative_paths
        self.parameters.use_relative_paths = tmp

        # Write out the modified files to the "_metadata.json" file
        fname = metadata_filename(self)
        with open(fname, 'w') as fd:
            json.dump({'output_files': modified_files}, fd)


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
