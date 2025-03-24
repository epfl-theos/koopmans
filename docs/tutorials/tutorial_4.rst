Tutorial 4: running convergence tests
=====================================
While ``koopmans`` is a package primarily oriented towards performing Koopmans functional calculations, it does have a couple of other useful functionalities. Among these functionalities is the ability to perform arbitrary covnergence tests.


A simple convergence test
-------------------------

In this tutorial, we will make use of its ``convergence`` task to determine how large a cell size and energy cutoff is required to converge the PBE energy of the highest occupied molecular orbital (HOMO) of a water molecule. This workflow was chosen for its simplicity; it is possible to run convergence tests on any workflow implemented in the ``koopmans`` package.

The input file
''''''''''''''
In order to run this calculation, in the ``workflow`` block we need to set the ``converge`` parameter to true:

.. literalinclude:: ../../tutorials/tutorial_4/h2o_conv.json
  :lines: 1-9
  :emphasize-lines: 6

and then provide the convergence information in the ``convergence`` block:

.. literalinclude:: ../../tutorials/tutorial_4/h2o_conv.json
  :lines: 10-17

These settings state that we are going to converge the HOMO energy to within 0.01 eV, with respect to *both* the energy cutoff ``ecutwfc`` and the cell size. The full input file can be found :download:`here <../../tutorials/tutorial_4/h2o_conv.json>`.

The output file
'''''''''''''''
When you run the calculation, you should see something like this after the header:

----

.. include:: ../_static/tutorials/tutorial_4/md_excerpts/h2o_conv.md
  :parser: myst_parser.sphinx_

----

Here, the code is attempting to use progressively larger energy cutoffs and cell sizes. It will ultimately arrive at a converged solution, with a ``ecutwfc`` of 50.0 Ha and a cell slightly larger than that provided in the ``.json`` input file.

Plotting
''''''''

The individual ``Quantum ESPRESSO`` calculations reside in nested subdirectories. If you plot the HOMO energies from each of these, you should get something like this:

.. figure:: ../../tutorials/tutorial_4/convergence.png
  :width: 640
  :alt: Plot of HOMO energy with respect to ecutwfc and cell size
  :align: center
  
  Plot of the HOMO energy of water with respect to the energy cutoff and cell size (generated using :download:`this script <../../tutorials/tutorial_4/plot.py>`)

We can see that indeed the calculation with ``ecutwfc = 50.0`` and ``celldm(1) = 13.3`` is the one where the energy of the HOMO goes within (and stays within) 0.01 eV of the most accurate calculation.

A custom convergence test
-------------------------

In the previous example, we performed a convergence test with respect to ``ecutwfc`` and ``celldm1``. A full list of supported convergence variables can be found :ref:`here <input_file/convergence_keywords:Valid keywords>`. You will see that only a couple of variables are implemented by default. But don't let that limit you! It is possible to perform a convergence test on arbitrary keywords using factories.

First, try taking the input file from the first part of the tutorial, and switch the pseudopotential library to ``pseudo_dojo_standard``. What do you notice?

Hopefully, the first thing you will see is that there are now some warnings about the small box parameters ``nrb`` like this,

.. code-block:: text

  🚨 Small box parameters `nrb` not provided in input: these will be automatically set to safe default values. These 
  values can probably be decreased, but this would require convergence tests. Estimated real mesh dimension `(nr1, 
  nr2, nr3) = 36 36 45`. Small box mesh dimension `(nr1b, nr2b, nr3b) = 20 18 20`.

These parameters are associated with the way ``Quantum ESPRESSO`` handles non-local core corrections in pseudopotentials, and corrections are present in the new set of pseudopotentials but absent in the SG15 pseudopotentials.

So, let's perform a convergence test! The ``nrb`` parameters are not directly implemented as a convergence variable in ``koopmans``, but we can use a factory to perform a convergence test on them, by making use of the ``ConvergenceVariable`` and ``ConvergenceWorkflowfactory`` classes.

.. literalinclude:: ../../tutorials/tutorial_4/custom_convergence.py

Running this script will perform a convergence test with respect to ``nrb 1-3``.

.. warning:: 
  This tutorial performs convergence tests in a slightly incorrect way, To see this, add the keyword ``length = 10`` to the ``ConvergenceVariable`` in the above script. How is the behavior different? Which behavior is correct? Why?