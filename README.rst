========
koopmans
========

| |GH Actions| |Coverage Status| |MIT License| |Documentation Status|

For performing Koopmans spectral functional calculations with ``Quantum ESPRESSO``

Directories
-----------
This repository contains...

| ``ase/`` a fork of ASE that manages reading and writing of ``Quantum ESPRESSO``
| ``bin/`` executables (N.B. this directory does not need to be added to ``$PATH``)  
| ``quantum_espresso/`` modified versions of ``Quantum ESPRESSO`` that contain implementations of the Koopmans functionals 
| ``pseudos/`` pseudopotentials
| ``requirements/`` python dependencies
| ``tests/`` test suite  
| ``workflows/`` source code of the workflow manager

Installation
------------

You are encouraged (but it is not necessary) to first create and activate a virtual environment as follows:

.. code-block:: bash

   sudo apt-get install python3-pip
   pip3 install virtualenv
   virtualenv ~/venvs/koopmans
   source ~/venvs/koopmans/bin/activate

For a quick installation one can then simply run ``make install``.

For a more careful installation, first make sure you have downloaded the various ``git`` submodules. To install this, run ``make submodules``, or equivalently

.. code-block:: bash

   git submodule init
   git submodule update

Then you need to compile the copies of ``Quantum ESPRESSO``. To do this, run

.. code-block:: bash

   make espresso MPIF90=<mpif90>

where ``<mpif90>`` should be replaced by the name of your chosen MPI Fortran90 compiler e.g. ``MPIF90=mpiifort``. The code should automatically detect and link the requisite libraries. (If this fails you will need to manually compile the two copies of ``Quantum ESPRESSO`` contained in the ``quantum_espresso/``.)

Finally, install the python workflow manager, either via ``make workflow``, or

.. code-block:: bash

   python3 -m pip install --upgrade pip
   python3 -m pip install -e . -e ase/

Running
-------
Calculations are run with the command

.. code-block:: bash

   koopmans <seed>.json

where <seed>.json is the ``koopmans`` input file. For a description of the contents of this file, refer to the documentation (`available online <https://koopmans-docs.readthedocs.io>`_). The keywords of ``koopmans`` keywords can be readily listed by running

.. code-block:: bash
   
   koopmans --help

In order to run the code in parallel, define the environment variables ``PARA_PREFIX`` and ``PARA_POSTFIX``. These are defined in the same way as in ``Quantum ESPRESSO``, e.g.

.. code-block:: bash

   export PARA_PREFIX="srun"
   export PARA_POSTFIX="-npool 4"

You can also *optionally* direct the code to use a central pseudopotentials directory by defining

.. code-block:: bash

   export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"

Contact
-------
Written and maintained by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna (2020-)

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

