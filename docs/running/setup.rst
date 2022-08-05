The setup block
^^^^^^^^^^^^^^^

The ``setup`` block contains variables common to all of the quantum espresso calculations e.g.

.. literalinclude:: ../../tutorials/tutorial_4/h2o_conv.json
  :lines: 12-36

It uses the same keywords and cards as ``pw.x`` input files, with namelists and cards are provided as subdictionairies. The one exception to this rule is the ``k_points`` block, which has a slightly customised layout:

.. literalinclude:: ../../tutorials/tutorial_2/si.json
  :lines: 49-61

``kpath`` specifies the path for band structures plotted during post-processing.

