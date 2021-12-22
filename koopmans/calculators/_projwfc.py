"""

projwfc.x calculator module for koopmans

"""

import os
from ase import Atoms
from koopmans.commands import Command, ParallelCommand
from koopmans.settings import ProjwfcSettingsDict
from ase.calculators.espresso import Projwfc
from ase.spectrum.dosdata import GridDOSData
from ase.spectrum.doscollection import GridDOSCollection
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory


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

        self.results_for_qc = ['DOS']
        if not isinstance(self.command, Command):
            self.command = ParallelCommand(os.environ.get(
                'ASE_PROJWFC_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = []

    def calculate(self):
        super().calculate()

        # Marija add code here to
        # 1) work out the names of the pDOS files to read
        #
        # 2) read in the pDOS file contents
        # dos = self.read_pdos(filename)
        #
        # 3) add pDOS to self.results
        # self.results['dos'] = GridDOSCollection(...)

    def read_pdos(filename: str) -> GridDOSData:
        # Marija: implement in this function how to extract from a DOS filename the contents of that file
        with open(filename, 'r') as fd:
            flines = fd.readlines()

        # Logic here to turn flines into a DOS object
        # dos = GridDOSData(...)

        return dos

    def is_complete(self):
        # Change me to check calc.results['job done']!
        return True

    def is_converged(self):
        return True
