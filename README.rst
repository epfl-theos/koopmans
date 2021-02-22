=========
python KI
=========

| |GH Actions| |Coverage Status| |MIT License| |Documentation Status|

For performing KI and KIPZ calculations with ``quantum espresso``

Directories
-----------
| ``ase_koopmans/`` a fork of ASE that manages reading and writing of ``quantum espresso`` (for installation details see below)
| ``scripts/`` scripts for performing useful tasks  
| ``koopmans/`` source code  
| ``pseudos/`` pseudopotentials  
| ``requirements/`` python dependencies
| ``tests/`` test suite  

Installation
------------

Having checked out the git repository, there are a few final steps to set up ``python KI``. First, make sure you have installed the submodule ``ase_koopmans``. To install this, run

.. code-block:: bash

   git submodule init
   git submodule update --remote --merge

Then, create and activate a virtual environment, for example:

.. code-block:: bash

   virtualenv ~/venvs/koopmans
   source ~/venvs/koopmans/bin/activate

And then finally install ``python_KI`` using ``pip``:

.. code-block:: bash

   python3 -m pip install --upgrade pip
   python3 -m pip install -e .

For each code you want to use (e.g. ``kcp.x``, ``pw.x``, etc) you can define the command to run this code with

.. code-block:: bash

   export ASE_ESPRESSO_KCP_COMMAND="srun kcp.x -in PREFIX.cpi > PREFIX.cpo"
   export ASE_ESPRESSO_COMMAND="srun pw.x -in PREFIX.pwi > PREFIX.pwo"
   export ASE_WANNIER90_COMMAND="wannier90.x PREFIX.win"
   export ASE_PW2WANNIER_COMMAND="srun pw2wannier90.x -in PREFIX.p2wi > PREFIX.p2wo"

Replace ``srun`` with whatever command is appropriate for your machine.

You can also *optionally* direct ASE to a central pseudopotentials directory by adding

.. code-block:: bash

   export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"

Running
-------
To see a list of options, run ``koopmans --help``

Contact
-------
Written and maintained by Edward Linscott and Riccardo De Gennaro (2020-)

For help and feedback email edward.linscott@gmail.com

.. |GH Actions| image:: https://img.shields.io/github/workflow/status/elinscott/python_KI/Run%20tests/master?label=master&logo=github
   :target: https://github.com/elinscott/python_KI/actions?query=branch%3Amaster
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/elinscott/python_KI/master?logo=codecov
   :target: https://codecov.io/gh/elinscott/python_KI
.. |MIT License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/elinscott/python_KI/blob/master/LICENSE
.. |Documentation Status| image:: https://readthedocs.org/projects/koopmans-docs/badge/?version=latest
   :target: https://koopmans-docs.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

