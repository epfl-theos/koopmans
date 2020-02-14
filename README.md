# `python KI`
For performing KI calculations with cp.x

## Directories
`bin` scripts for running calculations and performing useful tasks  
`examples` example calculations  
`koopmans_cp` source code  
`tests` test suite  

## Installation
Make sure you have the following python modules installed:
 * ``argparse``
 * ``pandas``
 * ``pickle``
 * ``warnings``
 * a modified version of ASE (``https://github.com/elinscott/ase_koopmans``; contact to be given access)

Add the following to your ~/.bashrc:

``export PYTHONPATH=/path/to/ase/:$PYTHONPATH``  
``export ASE_ESPRESSO_CP_COMMAND="srun cp.x -in PREFIX.cpi > PREFIX.cpo"`` (or whatever command is appropriate for your machine. Do not change anything after '``cp.x``')  
``export PYTHONPATH=/path/to/python_KI/:$PYTHONPATH``  
``export PATH=/path/to/python_KI/bin:$PATH``  
``export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"`` (optional)

## Running
To see a list of options, run ``./koopmans.py --help``

## Contact
Written by Edward Linscott, Dec 2019

For help and feedback email edward.linscott@gmail.com
