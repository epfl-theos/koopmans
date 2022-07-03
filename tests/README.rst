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
| ... as well as various directories containing the python scripts for the tests themselves. These mirror the file structure of ``koopmans`` itself

Generating a new benchmark
--------------------------

To generate new benchmarks, run with the additional flag ``--generate_benchmark``.

Creating new tests
------------------

To create a new test...

- identify which functionality you want to test, and in which module it is defined
- add a function called ``test_<name>`` to the relevant ``test_<module>.py`` file (the directory structure of ``tests/`` mirrors that of the source code)
- run ``pytest tests/ --generate_benchmark -k test_<name>`` to generate a benchmark
- check that ``pytest tests/``, ``pytest tests/ --mock``, and ``pytest tests/ --stumble`` pass

