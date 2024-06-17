.. _tutorial_6:

Tutorial 6: the band structure of bulk CrI3
=========================================================================

In this tutorial we will calculate the band structure of the bulk ferromagnetic semiconductor CrI3 using the supercell formulation of Koopmans. The input file for this tutorial can be downloaded :download:`here <../../tutorials/tutorial_6/cri3.json>`.


Calculating the Koopmans band structure
---------------------------------------

The input file
^^^^^^^^^^^^^^

First, let us inspect the input file:

.. literalinclude:: ../../tutorials/tutorial_6/cri3.json
  :lines: 1-10
  :lineno-start: 1
  :emphasize-lines: 4,6,10

Here we tell the code to calculate the KI bandstructure using the DSCF supercell cell approach. We will not actually calculate the screening parameters in this tutorial (because this calculation takes a bit of time) so we have set ``calculate_alpha`` to ``False`` and we have provided some reasonable screening parameters in the ``alpha_guess`` field. This value is given by the inverse of the (average of the) macroscopic dielectric function computed at DFT level. 

In its low temperature phase, bulk CrI3 is a ferromagnetic semiconductor with 3 unpaired d electrons on each of the two Cr atoms in the primitive cell. 
We provide this information to Koopmans by setting ``spin_polarized`` to ``true`` and by specifing the expected total magnetization in the calculator parameters. We also provide some initial values for the starting_magnetization on the Cr and I sites: 

.. literalinclude:: ../../tutorials/tutorial_6/cri3.json
  :lines: 45-53
  :lineno-start: 45
  :emphasize-lines: 3, 7,8

For a magnetic systems two sets of Wannier projections need to be provided, one for each spin channel. This is specified in the wannier parameter by splitting the ``w90`` block into an ``up`` and ``down`` sub-blocks:

.. literalinclude:: ../../tutorials/tutorial_6/cri3.json
  :lines: 55-75
  :lineno-start: 55
  :emphasize-lines: 2, 11

As already seen for ZnO in :ref:`Tutorial 3 <tutorial_3>`, we will use of the block-Wannierization functionality to wannierize each block of bands separately. 
The pojections provided above have been determined looking at the pDOS of the material as explained in :ref:`Tutorial 3 <projections_blocks_explanation>`. 
The first 4 up and down projections span the occupied manifold. Note that the 4th one differs between the two spin channels reflecting 
the different number of electrons in the up and down channels. The 5th projection in the up channel, and the 5th and 6th projections in the down channel span
the low lying part of the empty manifold. 

The rest of the file contains the atomic coordinates and k-point configuration.


Running the calculation
^^^^^^^^^^^^^^^^^^^^^^^

Running ``koopmans cri3.json`` should produce an output with several sections printing information on the different steps: after the header, the Wannierization is performed for the two spin channel, one at the time:

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 23-25
  :lineno-start: 23
  :language: text

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 38-40
  :lineno-start: 38
  :language: text


The Wannier functions are then written in the kcp format and the DFT inizialization is run. If we had instructed the code to calculate the alpha parameters, this would be followed by an extra block where these are calculated. But since we have already provided these, the workflow progresses immediately to the final step

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 90-92
  :lineno-start: 90
  :language: text

where the KI Hamiltonian is constructed and diagonalized in the supercell. To get the final band structure plot on the path specified in the input file, a last step is needed to unfold the SC bands into the Brillouin zone of the primitive cell: 

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 95-100
  :lineno-start: 95
  :language: text


Plotting the results
^^^^^^^^^^^^^^^^^^^^

To plot the final KI ban structure, we will load all of the information from the ``cri3.kwf`` file, as we already did for ZnO in :ref:`Tutorial 3 <tutorial_3>`. 
This is done in the ``plot_bands.py`` script

.. figure:: ../../tutorials/tutorial_6/cri3_bandstructures.png
  :width: 600
  :align: center
  :alt: the band structure plot for CrI3

  The band structure plot for CrI3

.. note:: 

   Altough qualitatively correct, this calculation and the final band structure is far from convergence. The energy cut-off and the supercell size needs to be increased for a quantitative calculation.
