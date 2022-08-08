Tutorial 3: KI on ZnO with DFPT
===============================

In this tutorial we will calculate the band structure of bulk zinc oxide using the k-space formulation of Koopmans. The input file for this tutorial can be downloaded :download:`here <../../tutorials/tutorial_3/zno.json>`.

Choosing Wannier projectors
---------------------------

The first step to any calculation on a periodic system is to obtain a good set of Wannier functions. These depend strongly on our choice of the projections, which (for the moment) we must specify manually.

To determine a good choice for the Wannier projections, we will first calculate a projected density of states (pDOS). The input file

.. literalinclude:: ../../tutorials/tutorial_3/zno.json
  :lines: 1-4
  :emphasize-lines: 3

is set up to perform this calculation already. Running ``koopmans zno.json`` will run the DFT bandstructure workflow. It will produce a directory called ``dft_bands`` that contains various files, including a ``png`` of the bandstructure and pDOS. If you look at this file, you will see that the DFT band structure of ZnO consists of several sets of bands, each well-separated in energy space. As the pDOS shows us, the filled bands correspond to zinc 3s, zinc 3p, oxygen 2s, and then zinc 3d hybridized with oxygen 2p. Meanwhile, the lowest empty bands correspond to Zn 4s bands.

A sensible choice for the occupied projectors is therefore

.. code:: json

    "w90": {
       "occ": {
       "projections_blocks": [
           [{"site": "Zn", "ang_mtm": "l=0"}],
           [{"site": "Zn", "ang_mtm": "l=1"}],
           [{"site": "O", "ang_mtm": "l=0"}],
           [{"site": "Zn", "ang_mtm": "l=2"},
            {"site": "O", "ang_mtm": "l=1"}]
       ]  
  
Here we will use of the ``projections_blocks`` functionality to wannierize each block of bands separately.

For the empty bands we want to obtain two bands corresponding to the Zn 4s orbitals. These must be disentangled from the rest of the empty bands, which is achieved via the following ``Wannier90`` keywords.

``dis_win_max``
  defines the upper bound of the disentanglement energy window. This window should entirely encompass the lowest two bands corresponding to our Zn 4s projectors. Consequently, it will inevitably include some weights from higher bands

``dis_froz_max``
  defines the upper bound of the frozen energy window. This window should be as large as possible while excluding any bands that do not correspond to our Zn 4s projectors

To determine good values for these keywords, we clearly need a more zoomed-in band structure than the default. We can obtain this via the ``*.fig.pkl`` files that ``koopmans`` generates. Here is a short code snippet that replots the band structure over a narrower energy range

.. literalinclude:: ../../tutorials/tutorial_3/replot_dft_bandstructure.py

Based on this figure, choose values for these two keywords and add them to your input file, alongside the definition of your projectors, as follows:

.. code-block:: json
   :emphasize-lines: 2-3

       "emp": {
          "dis_froz_max": "?",
          "dis_win_max": "?",
          "projections": [
              {"site": "Zn", "ang_mtm": "l=0"}
          ]  

.. note::

  ``dis_froz_max`` and ``dis_win_max`` should **not** be provided relative to the valence band edge. Meanwhile the band structure plots have set the valence band edge to zero. Make sure to account for this by shifting the values of ``dis_froz_max`` and ``dis_win_max`` by 9.3 eV (the valence band edge energy; you can get this value yourself via ``grep 'highest occupied level' dft_bands/scf.pwo``)

To test your wannierization, you can now switch to the ``wannierize`` task and once again run ``koopmans zno.json``. This will generate a ``wannier`` directory as well as a band structure plot, this time with the interpolated band structure plotted on top of the explicitly evaluated band structure. Ideally, the interpolated band structure should lie on top of the explicit one. Play around with the values of ``dis_froz_max`` and ``dis_win_max`` until you are happy with the resulting band structure.

.. hint::

  Instead of using the ``*.fig.pkl`` file to obtain a zoomed-in band structure, use the ``plotting`` block to manually set the y-limits:

  .. code:: json
    
    "plotting": {
      "Emin": -5.0
      "Emax": 15
    }

Koopmans calculation on ZnO
---------------------------

.. note::

  We suggest you proceed using values similar to ``"dis_froz_max": 14.5`` and ``"dis_win_max": 17.0``

Finally, we are ready to perform the Koopmans calculation itself. To do so, switch the task to ``singlepoint``. The output should look something like this:

