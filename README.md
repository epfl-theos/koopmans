# `python KI`
For performing KI and KIPZ calculations with ``cp.x``

## Directories
`ase_koopmans` a fork of ASE that manages `cp.x` reading and writing (for installation details see below)  
`bin` scripts for running calculations and performing useful tasks  
`examples` example calculations  
`koopmans_cp` source code  
`pseudos` pseudopotentials  
`tests` test suite  

## Installation
Make sure you have the following python modules installed:
 * ``argparse``
 * ``pandas``
 * ``pickle``
 * ``warnings``

You will also need to install the git repositories `testcode` and `ase_koopmans`. These have been added as submodules to this repository, so to install these run

```
git submodule init
git submodule update
```

Add the following to your ~/.bashrc:

```
export PATH=/path/to/python_KI/bin:$PATH 
export PYTHONPATH=/path/to/python_KI/:$PYTHONPATH  
export PYTHONPATH=/path/to/python_KI/ase_koopmans/:$PYTHONPATH  
export ASE_ESPRESSO_CP_COMMAND="srun cp.x -in PREFIX.cpi > PREFIX.cpo"
export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"
```

N.B. 
 - for ``ASE_ESPRESSO_CP_COMMAND``, replace ``srun`` with whatever command is appropriate for your machine. Do not change anything after '``cp.x``')
 - ``ESPRESSO_PSEUDO`` is an *optional* way of directing ASE to a central pseudopotentials directory

## Running
To see a list of options, run ``koopmans.py --help``

## Contact
Written by Edward Linscott, May 2020

For help and feedback email edward.linscott@gmail.com
