"""

pw2wannier calculator module for koopmans

Written by Edward Linscott Sep 2020

"""

import os
from ase import Atoms
from ase.calculators.espresso import PW2Wannier
from koopmans.settings import PW2WannierSettingsDict
from koopmans.commands import ParallelCommand
from ._utils import ExtendedCalculator, qe_bin_directory


class PW2WannierCalculator(ExtendedCalculator, PW2Wannier):

    # Adding all wannier90 keywords as decorated properties of the Wannier90Calculator class.
    # This means one can set and get wannier90 keywords as self.<keyword> but
    # internally they are stored as self._settings['keyword'] rather than
    # self.<keyword>

    def __init__(self, atoms: Atoms, *args, **kwargs):
        self.parameters = PW2WannierSettingsDict()

        # Define the appropriate file extensions
        self.ext_in = '.p2wi'
        self.ext_out = '.p2wo'

        # Initialise first using the ASE parent and then ExtendedCalculator
        PW2Wannier.__init__(self, atoms=atoms)
        ExtendedCalculator.__init__(self, *args, **kwargs)

        self.command = ParallelCommand(os.environ.get('ASE_PW2WANNIER_COMMAND',
                                                      str(qe_bin_directory) + os.path.sep + self.command))

    def is_converged(self):
        return True

    def is_complete(self):
        return self.results['job done']
