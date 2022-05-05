Running
=======
To run a calculation, all that is required is 

.. code-block:: bash

  $ koopmans <seed>.json

The input file
^^^^^^^^^^^^^^
The input file is a JSON-formatted file that contains the calculation parameters, divided into the :ref:`workflow <The workflow block>`, :ref:`setup <The setup block>`, :ref:`kcp <The kcp block>`, :ref:`pw <The pw block>`, :ref:`w90 <The w90 block>`, :ref:`pw2wannier <The pw2wannier block>`, and :ref:`ui <The ui block>` blocks.

.. toctree::
  :hidden:

  running/workflow
  running/setup
  running/kcp
  running/pw
  running/w90
  running/pw2wannier
  running/ui
  running/plot

Environment variables
^^^^^^^^^^^^^^^^^^^^^

.. include:: ../README.rst
    :start-line: 68
    :end-line: 80
