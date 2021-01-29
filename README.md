# `python KI`
For performing KI and KIPZ calculations with ``quantum espresso``

## Directories
`ase_koopmans` a fork of ASE that manages reading and writing of ``quantum espresso`` (for installation details see below)  
`bin` scripts for running calculations and performing useful tasks  
`koopmans` source code  
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
git submodule update --remote --merge
```

Add the following to your ~/.bashrc:

```
export PATH=/path/to/python_KI/bin:$PATH 
export PYTHONPATH=/path/to/python_KI/:$PYTHONPATH  
export PYTHONPATH=/path/to/python_KI/ase_koopmans/:$PYTHONPATH  
```

For each code you want to use (e.g. ``kcp.x``, ``pw.x``, etc) you can define the command to run this code with e.g.
```
export ASE_ESPRESSO_KCP_COMMAND="srun kcp.x -in PREFIX.cpi > PREFIX.cpo"
export ASE_ESPRESSO_COMMAND="srun pw.x -in PREFIX.pwi > PREFIX.pwo"
export ASE_WANNIER90_COMMAND="wannier90.x PREFIX.win"
export ASE_PW2WANNIER_COMMAND="srun pw2wannier90.x -in PREFIX.p2wi > PREFIX.p2wo"
```
Replace ``srun`` with whatever command is appropriate for your machine. Do not change anything after the binary names)

You can also *optionally* direct ASE to a central pseudopotentials directory by adding
```
export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"
```

## Running
To see a list of options, run ``run_koopmans.py --help``

## Contact
Written by Edward Linscott, May 2020

For help and feedback email edward.linscott@gmail.com
