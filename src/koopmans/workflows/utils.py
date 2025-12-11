"""Useful functions when constructing and manipulating Workflows."""

from collections import OrderedDict

from ase_koopmans import Atoms
from upf_tools.pydantic import UPFModel

from koopmans.engines.engine import Engine


def get_pseudopotentials(atoms: Atoms, engine: Engine, pseudo_library: str):
    """Fetch the pseudopotentials for the given system, using the specified engine and pseudopotential library."""
    # Work out the unique atom labels. This needs to be done like this rather than with set() to ensure
    # the order is consistent with how the atoms appear in the input file
    unique_labels = []
    for a in atoms:
        label_tuple: tuple[str, int] = (a.symbol, a.tag)
        if label_tuple not in unique_labels:
            unique_labels.append(label_tuple)

    # Fetch the pseudopotentials
    upfs: OrderedDict[str, UPFModel] = OrderedDict()
    for element, tag in unique_labels:
        pseudo = engine.get_pseudopotential(
            library=pseudo_library, element=element)

        pseudo_type = pseudo.content['header']['pseudo_type']
        if pseudo_type != 'NC':
            raise ValueError('Koopmans functionals only currently support norm-conserving pseudopotentials; '
                             f'{pseudo.filename} is {pseudo_type}')

        symbol = element + str(tag) if tag > 0 else element
        upfs[symbol] = pseudo

    return upfs
