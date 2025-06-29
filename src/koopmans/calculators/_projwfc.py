"""projwfc.x calculator module for koopmans."""

import copy
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from ase_koopmans import Atoms
from ase_koopmans.calculators.espresso import Projwfc
from ase_koopmans.spectrum.doscollection import GridDOSCollection
from ase_koopmans.spectrum.dosdata import GridDOSData
from upf_tools import UPFDict

from koopmans import pseudopotentials
from koopmans.files import File
from koopmans.settings import ProjwfcSettingsDict

from ._calculator import CalculatorABC, CalculatorExt


class ProjwfcCalculator(CalculatorExt, Projwfc, CalculatorABC):
    """Subclass of CalculatorExt for performing calculations with projwfc.x."""

    ext_in = '.pri'
    ext_out = '.pro'
    code = "projwfc"

    def __init__(self, atoms: Atoms, *args, **kwargs):
        # Define the valid settings
        self.parameters = ProjwfcSettingsDict()
        self.parent_process = None

        # Initialize first using the ASE parent and then CalculatorExt
        Projwfc.__init__(self, atoms=atoms)
        CalculatorExt.__init__(self, *args, **kwargs)

        # We need pseudopotentials and pseudo dir in order to work out the number of valence electrons for each
        # element, and therefore what pDOS files to expect. We also need spin-polarized to know if the pDOS files will
        # contain columns for each spin channel

        # These must be provided post-initialization (because we want to be allowed to initialize the calculator)
        # without providing these arguments
        self.pseudopotentials: Optional[Dict[str, UPFDict]] = None
        self.pseudo_dir: Optional[Path] = None
        self.spin_polarized: Optional[bool] = None

    def _pre_calculate(self):
        super()._pre_calculate()
        for attr in ['pseudopotentials', 'pseudo_dir', 'spin_polarized']:
            if not hasattr(self, attr):
                raise ValueError(f'Please set `{self.__class__.__name__}.{attr}` before calling '
                                 f'`{self.__class__.__name__}calculate()`')

    def _post_calculate(self):
        super()._post_calculate()
        self.generate_dos()

    @property
    def _expected_orbitals(self):
        """Generate a list of orbitals (e.g. 1s, 2s, 2p, ...) expected for each element in the system.

        This information can be found in the corresponding pseudopotential.
        """
        return pseudopotentials.expected_subshells(self.atoms, self.pseudopotentials)

    def generate_dos(self) -> GridDOSCollection:
        """Parse all of the pdos files and add these densities of state to self.results."""
        dos_list = []
        for atom in self.atoms:
            parent_directory = File(self, Path())
            filenames = sorted(parent_directory.glob(pattern=self.parameters.filpdos + f'.pdos_atm#{atom.index + 1}(*'))
            # The filename does not encode the principal quantum number n. In order to recover this number, we compare
            # the reported angular momentum quantum number l against the list of expected orbitals, and infer n
            # assuming only that the file corresponding to nl will come before (n+1)l
            label = atom.symbol + str(atom.tag) if atom.tag > 0 else atom.symbol
            orbitals = copy.copy(self._expected_orbitals[label])
            for filename in filenames:
                # Infer l from the filename
                subshell = filename.name.name[-2]
                # Find the orbital with matching l with the smallest n
                orbital = [o for o in orbitals if o[-1] == subshell][0]
                orbitals.remove(orbital)

                # Load the DOSs
                dos_list += self.read_pdos(filename, orbital)

        #  add pDOS to self.results
        self.results['dos'] = GridDOSCollection(dos_list)

    def read_pdos(self, filename: File, expected_subshell: str) -> List[GridDOSData]:
        """Read a pDOS file."""
        # Read the file contents
        assert self.engine is not None
        content = self.engine.read_file(filename)
        assert isinstance(content, str)
        flines = content.split('\n')

        # Parse important information from the filename
        [_, index, symbol, _, _, subshell, _] = re.split(r"#|\(|\)", filename.name.name)

        # Compare against the expected subshell
        if subshell != expected_subshell[1]:
            raise ValueError(
                f"Unexpected pdos file `{filename.name.name}`, a pdos file corresponding to {expected_subshell} "
                "was expected")

        # Work out what orbitals will be contained within the pDOS file
        orbital_order = {"s": ["s"], "p": ["pz", "px", "py"], "d": ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]}

        spins: List[Optional[str]]
        if self.spin_polarized:
            spins = ["up", "down"]
        else:
            spins = [None]
        orbitals = [(o, spin) for o in orbital_order[subshell] for spin in spins]

        # Parse the rest of the file
        data = np.array([l.split() for l in flines[1:-1]], dtype=float).transpose()
        energy = data[0]

        # Looping over each pDOS in the file...
        dos_list = []
        for weight, (label, spin) in zip(data[-len(orbitals):], orbitals):
            # Assemble all of the information about this particular pDOS
            info = {"symbol": symbol, "index": int(index), "n": expected_subshell[0], "l": subshell, "m": label}
            if self.spin_polarized:
                info['spin'] = spin

            # Create and store the pDOS as a GridDOSData object
            dos = GridDOSData(energy, weight, info=info)
            dos_list.append(dos)

        # Return the list of GridDOSData objects
        return dos_list

    def is_complete(self):
        """Return True if the calculation is complete."""
        return self.results['job done']

    def is_converged(self):
        """Return True; a projwfc calculation is always converged."""
        return True
