========
koopmans
========

| |GH Actions| |Coverage Status| |GPL License| |Documentation Status|

For performing Koopmans spectral functional calculations with ``Quantum ESPRESSO``.

Installation
------------

This python package ``koopmans`` only contains a subset of the code that you will require to run Koopmans functional calculations. To install everything, download the full source code from `github <https://github.com/epfl-theos/koopmans>` and follow the installation instructions.

Running
-------
Calculations are run with the command

.. code-block:: bash

   koopmans <seed>.json

where <seed>.json is the ``koopmans`` input file. For more details, refer to the `online documentation <https://koopmans-functionals.org>`_.

Contact
-------
Written and maintained by Edward Linscott, Riccardo De Gennaro, and Nicola Colonna (2020-)

For help and feedback email edward.linscott@gmail.com

.. |GH Actions| image:: https://img.shields.io/github/actions/workflow/status/epfl-theos/koopmans/tests.yml?master&label=master&logo=github
   :target: https://github.com/epfl-theos/koopmans/actions?query=branch%3Amaster
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/epfl-theos/koopmans/master?logo=codecov
   :target: https://codecov.io/gh/epfl-theos/koopmans
.. |GPL License| image:: https://img.shields.io/badge/license-GPL-blue
   :target: https://github.com/epfl-theos/koopmans/blob/master/LICENSE
.. |Documentation Status| image:: https://readthedocs.org/projects/koopmans/badge/?version=latest
   :target: https://koopmans-functionals.org/en/latest/?badge=latest
   :alt: Documentation Status

