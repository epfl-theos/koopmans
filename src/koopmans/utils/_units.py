"""Units for `koopmans.utils`."""


from ase_koopmans.units import create_units

# Quantum Espresso -- and koopmans -- uses CODATA 2006 internally
units = create_units('2006')
