"""

projwfc.x calculator module for koopmans

"""

import os
import numpy as np
import re
from glob import glob
from pathlib import Path
from typing import List, Dict, Optional
from ase import Atoms
from ase.calculators.espresso import Projwfc
from ase.spectrum.dosdata import GridDOSData
from ase.spectrum.doscollection import GridDOSCollection
from koopmans import pseudopotentials
from koopmans.commands import Command, ParallelCommand
from koopmans.settings import ProjwfcSettingsDict
from ._utils import CalculatorExt, CalculatorABC, qe_bin_directory


class ProjwfcCalculator(CalculatorExt, Projwfc, CalculatorABC):
    # Subclass of CalculatorExt for performing calculations with projwfc.x
    ext_in = '.pri'
    ext_out = '.pro'

    def __init__(self, atoms: Atoms, pseudopotentials: Dict[str, str], pseudo_dir: Path,
                 spin_polarised: bool, *args, **kwargs):
        # Define the valid settings
        self.parameters = ProjwfcSettingsDict()

        # Initialise first using the ASE parent and then CalculatorExt
        Projwfc.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        self.results_for_qc = ['dos']
        if not isinstance(self.command, Command):
            self.command = ParallelCommand(os.environ.get(
                'ASE_PROJWFC_COMMAND', str(qe_bin_directory) + os.path.sep + self.command))

        # We need pseudopotentials and pseudo dir in order to work out the number of valence electrons for each
        # element, and therefore what pDOS files to expect
        self.pseudopotentials = pseudopotentials
        self.pseudo_dir = pseudo_dir

        # We need spin-polarised to know if the pDOS files will contain columns for each spin channel
        self.spin_polarised = spin_polarised

        self.results_for_qc = []

    def calculate(self):
        super().calculate()
        self.generate_dos()

    @property
    def _expected_orbitals(self):
        """
        Generates a list of orbitals (e.g. 1s, 2s, 2p, ...) expected for each element in the system, based on the
        corresponding pseudopotential.
        """
        expected_orbitals = {}
        z_core_to_first_orbital = {0: '1s', 2: '2s', 4: '2p', 10: '3s', 12: '3p', 18: '4s', 20: '3d', 30: '4p'}
        for atom in self.atoms:
            if atom.symbol in expected_orbitals:
                continue
            pseudo_file = self.pseudopotentials[atom.symbol]
            z_core = atom.number - pseudopotentials.valence_from_pseudo(self.pseudo_dir, pseudo_file)
            first_orbital = z_core_to_first_orbital[z_core]
            all_orbitals = list(z_core_to_first_orbital.values())
            expected_orbitals[atom.symbol] = all_orbitals[all_orbitals.index(first_orbital):]
        return expected_orbitals

    def generate_dos(self) -> GridDOSCollection:
        """
        Parse all of the pdos files and add these densities of state to self.results
        """
        dos_list = []
        for atom in self.atoms:
            filenames = sorted(glob(self.parameters.filpdos + f'.pdos_atm#{atom.index+1}*'))
            orbitals = self._expected_orbitals[atom.symbol]
            for filename, orbital in zip(filenames, orbitals):
                dos_list += self.read_pdos(filename, orbital)

        #  add pDOS to self.results
        self.results['dos'] = GridDOSCollection(dos_list)

    def read_pdos(self, filename: str, expected_subshell: str) -> List[GridDOSData]:
        """
        Function for reading a pDOS file
        """

        # Read the file contents
        with open(filename, 'r') as fd:
            flines = fd.readlines()

        # Parse important information from the filename
        [_, index, symbol, _, _, subshell, _] = re.split(r"#|\(|\)", filename)

        # Compare against the expected subshell
        if subshell != expected_subshell[1]:
            raise ValueError(
                f"Unexpected pdos file {filename}, a pdos file corresponding to {expected_subshell} was expected")

        # Work out what orbitals will be contained within the pDOS file
        orbital_order = {"s": ["s"], "p": ["pz", "px", "py"], "d": ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]}

        spins: List[Optional[str]]
        if self.spin_polarised:
            spins = ["up", "down"]
        else:
            spins = [None]
        orbitals = [(o, spin) for o in orbital_order[subshell] for spin in spins]

        # Parse the rest of the file
        data = np.array([l.split() for l in flines[1:]], dtype=float).transpose()
        energy = data[0]

        # Looping over each pDOS in the file...
        dos_list = []
        for weight, (label, spin) in zip(data[-len(orbitals):], orbitals):
            # Assemble all of the information about this particular pDOS
            info = {"symbol": symbol, "index": int(index), "n": expected_subshell[0], "l": subshell, "m": label}
            if self.spin_polarised:
                info['spin'] = spin

            # Create and store the pDOS as a GridDOSData object
            dos = GridDOSData(energy, weight, info=info)
            dos_list.append(dos)

        # Return the list of GridDOSData objects
        return dos_list

    def is_complete(self):
        return self.results['job done']

    def is_converged(self):
        return True
