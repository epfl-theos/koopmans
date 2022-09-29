How to run
==========
To run a calculation from the command line, all that is required is 

.. code-block:: bash

  $ koopmans <seed>.json

where ``<seed>.json`` is a ``koopmans`` input file. The format of the input file is documented :ref:`here <input_file:The input file>`.

Running via python
^^^^^^^^^^^^^^^^^^
It is possible to run ``koopmans`` workflows from within python, bypassing the need for an input file entirely. To do this, all you need to do is create a ``SinglepointWorkflow`` object

.. code-block:: python

  wf = SinglepointWorkflow(...)

and then simply call

.. code-block:: python

  wf.run()

For details of how to initialize a workflow object, see the :ref:`workflow class documentation <modules:The workflow module>`. After a calculation has finished, you can access the individual calculations e.g.

.. code-block:: python

  final_calc = wf.calculations[-1]

and fetch their results e.g.

.. code-block:: python

  total_energy = final_calc.results['energy']

.. include:: ../README.rst
    :start-line: 102
    :end-line: 129

