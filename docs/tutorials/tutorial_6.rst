.. _tutorial_6:

Tutorial 3: the band structure of bulk CrI3
=========================================================================

In this tutorial we will calculate the band structure of the bulk ferromagnetic semiconductor CrI3 using supercell formulation of Koopmans. The input file for this tutorial can be downloaded :download:`here <../../tutorials/tutorial_6/cri3.json>`.

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
CrI3 in its bulk low temperature phase is a Ferromagnetic semiconductor with 3 umpaired d electrons on each of the two Cr atoms in the primitive cell. We provide this information to Koopmans by setting ``spin_polarized`` to ``true`` and by specify the expected total magnetization in the calculator parameters. We also provide some initial values for the starting_magnetization on the Cr and I sites: 

.. literalinclude:: ../../tutorials/tutorial_6/cri3.json
  :lines: 45-53
  :lineno-start: 45
  :emphasize-lines: 47, 51, 52


The rest of the file contains the atomic coordinates, k-point configuration, and the Wannier projectors (one set for each spin channel), which we will discuss later.

Running the calculation
^^^^^^^^^^^^^^^^^^^^^^^

Running ``koopmans cri3.json`` should produce an output with several sections: after the header the Wannierization is performed, the Wannier functions are written in the kcp format and the DFT inizialization.
If we had instructed the code to calculate the alpha parameters, this would be followed by an extra block where these are calculated. But since we have already provided these, the workflow progresses immediately to the final step

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 94-100
  :lineno-start: 94
  :language: text

where the KI Hamiltonian is constructed and diagonalized in the supercell. To get the final band structure plot on the path specified in the input file, a last step is needed to unfold the SC bands into the Brillouin zone of the primitive cell: 

.. literalinclude:: ../../tutorials/tutorial_6/cri3.out
  :lines: 103-108
  :lineno-start: 103
  :language: text


Plotting the results
^^^^^^^^^^^^^^^^^^^^

To plot the final KI ban structure, we will load all of the information from the ``cri3.kwf`` file, as we already did for ZnO in tutorial_3. 
This is done in the ``plot_bands.py`` script

.. figure:: ../../tutorials/tutorial_6/cri3_bandstructures.png
  :width: 600
  :align: center
  :alt: the band structure plot for CrI3

  The band structure plot for CrI3

