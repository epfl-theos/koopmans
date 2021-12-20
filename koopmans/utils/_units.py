'''

Units koopmans.utils

Written by Edward Linscott May 2020

'''


from ase.units import create_units

# Quantum Espresso -- and koopmans -- uses CODATA 2006 internally
units = create_units('2006')
