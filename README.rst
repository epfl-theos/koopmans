========
koopmans
========

| |GH Actions| |Coverage Status| |GPL License| |Documentation Status|

For performing Koopmans spectral functional calculations with ``Quantum ESPRESSO``

Directories
-----------
This repository contains...

| ``bin/`` executables (N.B. this directory does not need to be added to ``$PATH``)  
| ``docs/`` documentation (see https://koopmans-functionals.org/)  
| ``src/`` source code
| ``quantum_espresso/`` modified versions of ``Quantum ESPRESSO`` that contain implementations of the Koopmans functionals 
| ``pseudos/`` pseudopotentials
| ``requirements/`` python dependencies
| ``tests/`` test suite  
| ``tutorials/`` tutorials  

Installation
------------

Quick installation
^^^^^^^^^^^^^^^^^^
For a quick installation one can simply run ``make; sudo make install``

Detailed installation
^^^^^^^^^^^^^^^^^^^^^

Setting up a virtual environment
""""""""""""""""""""""""""""""""

You are encouraged (but it is not necessary) to first create and activate a virtual environment as follows:

.. code-block:: bash

   sudo apt-get install python3-pip
   pip3 install virtualenv
   virtualenv ~/venvs/koopmans
   source ~/venvs/koopmans/bin/activate

Note that ``koopmans`` requires python v3.7 or later. If your computer's version of python3 corresponds to an earlier version, install python v3.7 or later, and then direct ``virtualenv`` to create the virtual environment using that specific installation of python via

.. code-block:: bash

   virtualenv ~/venvs/koopmans -p /usr/bin/python3.x

Fetching the submodules
"""""""""""""""""""""""

Now, ensure you have downloaded the various ``git`` submodules. To do so, run ``make submodules``, or equivalently

.. code-block:: bash

   git submodule init
   git submodule update

Compiling Quantum ESPRESSO
""""""""""""""""""""""""""

Then you need to compile the copies of ``Quantum ESPRESSO``. To do this, run

.. code-block:: bash

   make espresso MPIF90=<mpif90>

where ``<mpif90>`` should be replaced by the name of your chosen MPI Fortran90 compiler e.g. ``MPIF90=mpiifort``. The code should automatically detect and link the requisite libraries. (If this fails you may need to manually compile the two versions of ``Quantum ESPRESSO`` contained in the ``quantum_espresso/`` directory.)

Adding Quantum ESPRESSO to your path
""""""""""""""""""""""""""""""""""""

To add all of the Quantum ESPRESSO binaries to your path, run

.. code-block:: bash

   sudo make install

By default this will copy the Quantum ESPRESSO binaries to ``/usr/local/bin``. This requires sudo privileges. If you do not have sudo privileges, you can either (a) install the codes in a different location by running ``make install PREFIX=/path/to/bin/`` (substitute ``/path/to/bin/`` with any directory of your choosing that is on your path) or (b) append ``bin/`` from the current directory to your path.

Installing the workflow manager
"""""""""""""""""""""""""""""""

Finally, install the python workflow manager, either via ``make workflow``, or

.. code-block:: bash

   python3 -m pip install --upgrade pip
   python3 -m pip install -e .

Running
-------
Calculations are run with the command

.. code-block:: bash

   koopmans <seed>.json

where <seed>.json is the ``koopmans`` input file. For more details, refer to the `online documentation <https://koopmans-docs.readthedocs.io>`_.

Parallelism
^^^^^^^^^^^

In order to run the code in parallel, define the environment variables ``PARA_PREFIX`` and ``PARA_POSTFIX``. These are defined in the same way as in ``Quantum ESPRESSO``, e.g.

.. code-block:: bash

   export PARA_PREFIX="srun"
   export PARA_POSTFIX="-npool 4"

Pseudopotentials
^^^^^^^^^^^^^^^^

Currently, Koopmans functionals only works with norm-conserving pseudopotentials. We suggest you use optimized norm-conserving Vanderbilt pseudopotentials, such as

- the `SG15 library <http://www.quantum-simulation.org/potentials/sg15_oncv/index.htm>`_
- the `Pseudo Dojo library <http://www.pseudo-dojo.org/index.html>`_

For convenience, ``koopmans`` already ships with both of these pseudopotential libraries and you can simply select the one you want to use using the ``pseudo_library`` keyword.

If you prefer to use your own pseudopotentials, add them to ``src/koopmans/pseudopotentials/<my_pseudos>/<functional>``, where ``<my_pseudos>`` is a name of your choosing and ``<functional>`` is the functional used to generate your pseudopotentials. You can then direct ``koopmans`` to use these pseudopotentials by setting the keywords ``pseudo_library`` and ``base_functional`` to ``<my_pseudos>`` and ``<functional>`` respectively.

Alternatively, you can direct the code to always use your personal pseudopotentials directory by defining the variable

.. code-block:: bash

   export ESPRESSO_PSEUDO="/path/to/pseudopotential/folder/"

Contact
-------
Written and maintained by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna (2020-)

For help and feedback email edward.linscott@gmail.com

.. |GH Actions| image:: https://img.shields.io/github/workflow/status/epfl-theos/koopmans/Run%20tests/master?label=master&logo=github
   :target: https://github.com/epfl-theos/koopmans/actions?query=branch%3Amaster
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/epfl-theos/koopmans/master?logo=codecov
   :target: https://codecov.io/gh/epfl-theos/koopmans
.. |GPL License| image:: https://img.shields.io/badge/license-GPL-blue
   :target: https://github.com/epfl-theos/koopmans/blob/master/LICENSE
.. |Documentation Status| image:: https://readthedocs.org/projects/koopmans/badge/?version=latest
   :target: https://koopmans-functionals.org/en/latest/?badge=latest
   :alt: Documentation Status

