.. _tutorial_2:

Tutorial 2: the band structure of bulk silicon (calculated via a supercell)
===========================================================================
In this tutorial, we will calculate the KI bandstructure of bulk silicon. The input file used for this calculation can be downloaded :download:`here <../../tutorials/tutorial_2/si.json>`.

Wannierization
--------------
While this calculation will bear a lot of similarity to the previous tutorial, there are several differences between performing a Koopmans calculation on a bulk system vs. a molecule. One major difference is that we use Wannier functions as our variational orbitals.

What are Wannier functions?
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Most electronic-structure codes try to calculate the Bloch states :math:`\psi_{n\mathbf{k}}` of periodic systems (where :math:`n` is the band index and :math:`\mathbf{k}` the crystal momentum). However, other representations are equally valid. One such representation is the Wannier function (WF) basis. In contrast to the very delocalised Bloch states, WFs are spatially localised and as such represent a very convenient basis to work with. In our case, the fact that they are localised means that they are suitable for use as variational orbitals.

Wannier functions :math:`w_{n\mathbf{R}}(\mathbf{r})` can be written in terms of a transformation of the Bloch states:

.. math::

  \begin{equation}
  w_{n \mathbf{R}}(\mathbf{r})
  =\frac{V}{(2 \pi)^{3}}
  \int_{\mathrm{BZ}}
  \left[
    \sum_{m} U_{m n}^{(\mathbf{k})} \psi_{m \mathbf{k}}(\mathbf{r})
  \right]
  e^{-i \mathbf{k}\cdot \mathbf{R}}
  \mathrm{d} \mathbf{k}
  \end{equation}

where our Wannier functions belong to a particular lattice site :math:`\mathbf{R}`, :math:`V` is the unit cell volume, the integral is over the Brillouin zone (BZ), and :math:`U^{(\mathbf{k})}_{mn}` defines a unitary rotation that mixes the Bloch states with crystal momentum :math:`\mathbf{k}`. Crucially, this matrix :math:`U^{(\mathbf{k})}_{mn}` is not uniquely defined --- indeed, it represents a freedom of the transformation that we can exploit.

We choose our :math:`U^{(\mathbf{k})}_{mn}` that gives rise to WFs that are "maximially localised". We quantify the spread :math:`\Omega` of a WF as

.. math::

  \begin{equation}
  \Omega=
  \sum_{n}
  \left[
    \left\langle w_{n \mathbf{0}}(\mathbf{r}) \right|
    r^{2}
    \left| w_{n \mathbf{0}}(\mathbf{r})\right\rangle
    -
    \left|
    \left\langle w_{n \mathbf{0}}(\mathbf{r})\right|
    \mathbf{r}
    \left| w_{n \mathbf{0}}(\mathbf{r})\right\rangle
    \right|^{2}
    \right]
  \end{equation}

The Wannier functions that minimise this spread are called maximally localised Wannier functions (MLWFs). For more details, see Ref. :cite:`Marzari2012`.

How do I calculate Wannier functions?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MLWFs can be calculated with `Wannier90 <http://www.wannier.org/>`_, an open-source code that is distributed with ``Quantum ESPRESSO``. Performing a Wannierization with Wannier90 requires a series of calculations to be performed with ``pw.x``, ``wannier90.x``, and ``pw2wannier90.x``. This workflow is automated within ``koopmans``, as we will see in this tutorial.

.. note::

  This tutorial will not discuss in detail how perform a Wannierization with Wannier90. The Wannier90 website contains lots of excellent `tutorials <http://www.wannier.org/support/>`_.

  One important distinction to make for Koopmans calculations --- as opposed to many of the Wannier90 tutorials --- is that we need to separately Wannierize the occupied and empty manifolds.

Let's now inspect this tutorial's :download:`input file <../../tutorials/tutorial_2/si.json>`. At the top you will see that

.. literalinclude:: ../../tutorials/tutorial_2/si.json
  :lines: 2-4
  :lineno-start: 2
  :emphasize-lines: 2

which tells the code to perform a standalone Wannierization calculation. Meanwhile, near the bottom of the file there are some Wannier90-specific parameters provided in the ``w90`` block

.. literalinclude:: ../../tutorials/tutorial_2/si.json
  :lines: 45-52
  :lineno-start: 45

Here, the keywords provided in the ``emp`` subdictionary are only applied during the Wannierization of the empty manifold. The ``w90`` block format is explained more fully :ref:`here <input_file:The w90 subblock>`.

We run this calculation as per usual:

.. code-block:: bash

  koopmans si.json | tee si_wannierize.out

After the usual header, you should see something like the following:

.. literalinclude:: ../../tutorials/tutorial_2/si_wannierize.out
  :lines: 15-27

These various calculations that are required to obtain the MLWFs of bulk silicon. You can inspect the  and various output files will have been generated in a new ``wannier/`` directory.

.. collapse:: Click here for detailed descriptions of each calculation

  scf
    a ``pw.x`` self-consistent DFT calculation performed with no empty bands. This obtains the ground-state electronic density

  nscf
    a ``pw.x`` non-self-consistent DFT calculation that determines the Hamiltonian, now including some empty bands

  occ/wann_preproc
    a preprocessing ``wannier90.x`` calculation that generates some files required by ``pw2wannier90.x``

  occ/pw2wan
    a ``pw2wannier90.x`` calculation that extracts from the eariler ``pw.x`` calculations several key quantities required for generating the Wannier orbitals for the occupied manifold: namely, the overlap matrix of the cell-periodic part of the Block states (this is the ``wann.mmn`` file) and the projection of the Bloch states onto some trial localised orbitals (``wann.amn``)

  occ/wann
    the ``wannier90.x`` calculation that obtains the MLWFs for the occupied manifold

  emp/...
    the analogous calculations as those in ``occ/``, but for the empty manifold
  
  bands
    a ``pw.x`` calculation that calculates the band structure of silicon explicitly, used for verification of the Wannierization (see the next section)

|

The main output files of interest in ``wannier/`` are files ``occ/wann.wout`` and ``emp/wann.wout``, which contain the output of ``wannier90.x`` for the Wannierization of the occupied and empty manifolds. If you inspect either of these files you will be able to see a lot of starting information, and then under a heading like

.. code-block::

  *------------------------------- WANNIERISE ---------------------------------*
  +--------------------------------------------------------------------+<-- CONV
  | Iter  Delta Spread     RMS Gradient      Spread (Ang^2)      Time  |<-- CONV
  +--------------------------------------------------------------------+<-- CONV

you will then see a series of steps where you can see the Wannier functions being optimised and the spread (labelled ``SPRD``) decreasing from one step to the next. Scrolling down further you should see a statement that the Wannierization procedure has converged, alongside with a summary of the final state of the WFs.

How do I know if the Wannier functions I have calculated are "good"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Performing a Wannierization calculation is not a straightforward procedure, and requires the tweaking of the Wannier90 input parameters in order to obtain a "good" set of Wannier functions.

One check is to see if an interpolated bandstructure generated by the MLWFs resembles an explicitly-calculated band structure. (For an explanation of how one can use Wannier functions to interpolate band structures, refer to Ref. :cite:`Marzari2012`.) You might have noticed that when we ran ``koopmans`` earlier it also generated a file called ``si_wannierize_bandstructure.png``. It should look something like this:

.. figure:: ../../tutorials/tutorial_2/2x2x2/si_wannierize_bandstructure.png
  :width: 400
  :align: center
  :alt: comparison of the interpolated and explicitly-calculated band structures of silicon

  Comparing the interpolated and explicitly-calculated band structures of bulk silicon

Clearly, we can see that the interpolation is no good! (The interpolated band structure ought to lie on top of the explicitly calculated band structure.) The reason for this is that the Brillouin zone is undersampled by our :math:`2\times2\times2` :math:`k`-point grid. Try increasing the size of the k-point grid and see if the interpolated bandstructure improves.

.. tip::

  Trying different grid sizes can be very easily automated within ``python``. Here is a simple script that will run the Wannierization for three different grid sizes:

  .. literalinclude:: ../../tutorials/tutorial_2/wannierize.py

.. tip::
  You may also notice a file called ``si_wannierize_bandstructure.fig.pkl``. This is a version of the figure that you can load in python and modify as you see fit. e.g. here is a script that changes the y-axis limits and label:

  .. literalinclude:: ../../tutorials/tutorial_2/load_pickled_figure.py

The KI calculation
------------------
Having obtained a Wannierization of silicon that we are happy with, we can proceed with the KI calculation. In order to do this simply change the ``task`` keyword in the input file from ``wannierize`` to ``singlepoint``.

.. tip::

  Although we just discovered that a :math:`2\times2\times2` :math:`k`-point grid was inadequate for producing good Wannier functions, this next calculation is a lot more computationally intensive and will take a long time on most desktop computers. We therefore suggest that for the purposes of this tutorial you switch back to the small :math:`k`-point grid. (But for any proper calculations, always use high-quality Wannier functions!)

Initialization
^^^^^^^^^^^^^^
If you run this new input the output will be remarkably similar to that from the previous tutorial, with a couple of exceptions. At the start of the workflow you will see there is a Wannierization procedure, much like we had earlier when we running with the ``wannierize`` task:

.. literalinclude:: ../../tutorials/tutorial_2/si_ki.out
  :lines: 16-27
  :lineno-start: 16
  :language: text

which replaces the previous series of semi-local and PZ calculations that we used to initialize the variational orbitals for a molecule.

There is then an new "folding to supercell" subsection:

.. literalinclude:: ../../tutorials/tutorial_2/si_ki.out
  :lines: 28-31
  :lineno-start: 28
  :language: text

In order to understand what these calculations are doing, we must think ahead to the next step in our calculation, where we will calculate the screening parameters using the ΔSCF method. These calculations, where we remove/add an electron from/to the system, require us to work in a supercell. This means that we must transform the :math:`k`-dependent primitive cell results from previous calculations into equivalent :math:`\Gamma`-only supercell quantities that can be read by ``kcp``. This is precisely what the above ``wan2odd`` calculations do.

Calculating the screening parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Having transformed into a supercell, the calculation of the screening parameters proceeds as usual. The one difference to tutorial 1 that you might notice at this step is that we are skipping the calculation of screening parameters for some of the orbitals:

.. literalinclude:: ../../tutorials/tutorial_2/si_ki.out
  :lines: 35-45
  :lineno-start: 35
  :emphasize-lines: 7
  :language: text

The code is doing this because of what we provided for the ``orbital_groups`` in the input file:

.. literalinclude:: ../../tutorials/tutorial_2/si.json
  :lines: 10-12
  :lineno-start: 10
  :emphasize-lines: 2

which tells the code to use the same parameter for orbitals belonging to the same group. In this instance we are calculating a single screening parameter for all four filled orbitals, and a single screening parameter for the empty orbitals.

The final calculation and postprocessing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The final difference for the solids calculation is that there is an additional preprocessing step at the very end:

.. literalinclude:: ../../tutorials/tutorial_2/si_ki.out
  :lines: 72-
  :lineno-start: 72
  :language: text

Here, we transform back our results from the supercell sampled at :math:`\Gamma` to the primitive cell with :math:`k`-space sampling. This allows us to obtain a bandstructure. The extra Wannierization step that is being performed is to assist the interpolation of the band structure in the primitive cell, and has been performed because in the input file we specified

.. literalinclude:: ../../tutorials/tutorial_2/si.json
  :lines: 53-55
  :lineno-start: 53
  :emphasize-lines: 2

For more details on the "unfold and interpolate" procedure see :ref:`here <input_file:The ui subblock>` and Ref. :cite:`DeGennaro2022`.

Extracting the KI bandstructure and the bandgap of Si
-----------------------------------------------------
The bandstructure can be found in ``postproc/bands_interpolated.dat`` as a raw data file, but there is a more flexible way for plotting the final bandstructure using the python machinery of ``koopmans``:

.. literalinclude:: ../../tutorials/tutorial_2/plot_bandstructure.py

Running this script will generate a plot of the bandstructure (``si_bandstructure.png``) as well as printing out the band gap. You should get a result around 1.35 eV. Compare this to the PBE result of 0.68 eV and the experimental value of 1.22 eV. If we were more careful with the Wannier function generation, our result would be even closer (indeed in Ref. :cite:`Nguyen2018` the KI band gap was found to be 1.22 eV!) 
