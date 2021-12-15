"""

projwfc.x calculator module for koopmans

"""

from ._utils import CalculatorExt, CalculatorABC


class ProjwfcCalculator(CalculatorExt, Projwfc, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with projwfc.x
    ext_in = '.pri'
    ext_out = '.pro'

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = ProjwfcSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Projwfc.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.results_for_qc = ['energy']
        if not isinstance(self.command, Command):
            self.command = ParallelCommandWithPostfix(os.environ.get(
                'ASE_PROJWFC_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = []
