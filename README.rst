=========
python KI
=========

| |GH Actions| |Coverage Status| |MIT License| |Documentation Status|

For performing KI and KIPZ calculations with ``Quantum ESPRESSO``

Directories
-----------
| ``ase_koopmans/`` a fork of ASE that manages reading and writing of ``Quantum ESPRESSO`` (for installation details see below)
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

In order to run the code in parallel, define the environment variables ``PARA_PREFIX`` and ``PARA_POSTFIX``. These are defined in the same way as in Quantum ESPRESSO, e.g.

.. code-block:: bash

   export PARA_PREFIX="srun"
   export PARA_POSTFIX="-npool 4"

You can also *optionally* direct the code to use a central pseudopotentials directory by defining

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

