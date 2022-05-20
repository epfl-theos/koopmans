Test suite for koopmans
=======================

Quick start
-----------

To run the tests, run

.. code-block:: bash

    pytest tests/

from the base ``koopmans`` directory (i.e. the parent directory of this one)

Directories
-----------
This directory contains...

| ``benchmarks/`` the benchmarks against which the tests are compared
| ``data/`` Quantum ESPRESSO files required by the tests
| ``tmp/`` the temporary directory where the results of tests can be found
| ... as well as various directories containing the python scripts for the tests themselves

Generating a new benchmark
--------------------------

To generate new benchmarks, run with the additional flag ``--generate_benchmark``