"""

projwfc.x calculator module for koopmans

"""

import os
import numpy as np
import re
from ase import Atoms
from koopmans.commands import Command, ParallelCommand
from koopmans.settings import ProjwfcSettingsDict
from ase.calculators.espresso import Projwfc
from ase.spectrum.dosdata import GridDOSData
from ase.spectrum.doscollection import GridDOSCollection
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory
from glob import glob


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

        self.results_for_qc = ['dos']
        if not isinstance(self.command, Command):
            self.command = ParallelCommand(os.environ.get(
                'ASE_PROJWFC_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        self.results_for_qc = []

    def calculate(self):
        super().calculate()
        self.generate_dos()

    def generate_dos(self):
        dos_list = []
        for atom in self.atoms:

            for filename, orbital in zip(sorted(glob(self.parameters.filpdos + f'.pdos_atm#{atom.index+1}*')), self.expected_orbitals[atom.symbol]):
                dos_list += self.read_pdos(filename, orbital)

        #  add pDOS to self.results
        self.results['dos'] = GridDOSCollection(dos_list)

    def read_pdos(self, filename: str, expected_subshell: str) -> GridDOSData:
        # Marija: implement in this function how to extract from a DOS filename the contents of that file
        with open(filename, 'r') as fd:
            flines = fd.readlines()
        [_, index, symbol, _, _, subshell, _] = re.split("#|\(|\)", filename)
        if subshell != expected_subshell[1]:
            raise ValueError(
                f"Unexpected pdos file {filename}, a pdos file corresponding to {expected_subshell} was expected")
        dos_list = []
        data = np.array([l.split() for l in flines[1:]], dtype=float).transpose()
        energy = data[0]
        orbital_order = {"s": ["s"], "p": ["pz", "px", "py"], "d": ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]}
        orbitals = [(o, spin) for o in orbital_order[subshell] for spin in ["up", "down"]]
        for weight, (label, spin) in zip(data[-len(orbitals):], orbitals):
            dos = GridDOSData(energy, weight, info={"symbol": symbol, "index": int(index), "n": expected_subshell[0], "l": subshell, "m": label,
                                                    "spin": spin})
            dos_list.append(dos)
        return dos_list

    def is_complete(self):
        return self.results['job done']

    def is_converged(self):
        return True
