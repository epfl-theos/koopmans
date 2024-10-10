.. _tutorial_1:

Tutorial 1: the ionization potential and electron affinity of ozone
===================================================================
In this tutorial, we will calculate the ionisation potential and electron affinity of ozone.

The input
---------

The input file for this calculation can be downloaded :download:`here <../../tutorials/tutorial_1/ozone.json>`. Just to briefly highlight the most important details of the workflow block

.. literalinclude:: ../../tutorials/tutorial_1/ozone.json
  :lines: 2-4
  :lineno-start: 2
  :emphasize-lines: 2

here we select the KI functional (as opposed to KIPZ),

.. literalinclude:: ../../tutorials/tutorial_1/ozone.json
  :lines: 3-5
  :lineno-start: 3
  :emphasize-lines: 2

specifies that we are going to calculate the screening parameters via a ΔSCF procedure, whereby we compute the energies of various :math:`N`, :math:`N-1`, and :math:`N+1`-electron systems (see :ref:`the theory section<theory_dscf>` for details), and

.. literalinclude:: ../../tutorials/tutorial_1/ozone.json
  :lineno-start: 4
  :lines: 4-6
  :emphasize-lines: 2

specifies that we have chosen to use the Kohn-Sham orbitals as our :ref:`variational orbitals <theory_vorbs_vs_corbs>`. This is common practice for molecules.

Meanwhile, the ``atoms`` block describes the both the cell and the atoms it contains. If you are familiar with ``Quantum ESPRESSO`` input files then most of this should look very familiar to you (albeit in JSON format).

Running the calculation
------------------------
In order to run the calculation, simply run

.. code-block:: bash

  koopmans ozone.json | tee ozone.md

.. tip::
  In order to run in parallel, set the ``PARA_PREFIX`` environment variable to ``mpirun -np 4`` (or similar)

The outline
-----------
First, let us inspect the contents of ``ozone.md``. This provides us with an outline of the workflow that we ran. It is in `markdown <https://commonmark.org/help/>`_ format, which is rendered nicely by modern code editors such as vscode.

After the header we can see there are a list of calculations that have been run by ``koopmans``. These come under three headings.

.. _tutorial_1_initialization:

Initialization
^^^^^^^^^^^^^^
The first step in any Koopmans calculation is the initialization step. In this step we initialize the density and the variational orbitals.

----

.. include:: ../_static/tutorials/tutorial_1/md_excerpts/ozone_init.md
  :parser: myst_parser.sphinx_

----

For this calculation we can see that ``koopmans`` has run four calculations. These initialize the density with the PBE density. Indeed, from this point onwards in the calculation the density will never change, because the KI functional yields the same density as the base functional. (N.B. This is not true of KIPZ.)

These PBE calculations also have provided us with our variational orbitals -- we can see that the calculation has selected the PBE Kohn-Sham orbitals as the variational orbitals (because earlier we set ``"init_orbitals": "kohn-sham"``).

But why four calculations? The reason for this is that the calculations we will perform later involve the addition/removal of a single electron, which means the density we need to generate here must correspond to a ``nspin = 2`` calculation. However, we know ozone is a closed-shell molecule and simply performing a ``nspin = 2`` PBE calculation risks introducing spin contamination (i.e. falling into a local minimum where :math:`n^\uparrow(\mathbf{r}) \neq n^\downarrow(\mathbf{r})`).

This sequence of four calculations is designed to avoid this; we first optimize the density constrained to be spin-unpolarized, and only once that density has been minimized do we lift this restriction. This additional step can be disabled by adding ``"fix_spin_contamination": false`` to the ``workflow`` block of ``ozone.json``.

The input and output ``Quantum ESPRESSO`` files can all be found in the corresponding directory ``01-koopmans-dscf/01-initialization/``.

.. _tutorial_1_screening_parameters:

Calculating the screening parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The second step in the calculation involves the calculation of the screening parameters:

----

.. include:: ../_static/tutorials/tutorial_1/md_excerpts/ozone_alpha.md
  :parser: myst_parser.sphinx_

----

etc. Here, we are calculating the screening parameters using the :ref:`ΔSCF method <theory_dscf>`. For filled orbitals (orbitals 1-9 of ozone) this requires an :math:`N-1`-electron PBE calculation where we freeze the i\ :sup:`th` orbital, empty it, and allow the rest of the density to relax. This yields :math:`E_i(N-1)`. Meanwhile, :math:`E(N)`, :math:`\varepsilon_{i}^\alpha(1)`, and :math:`\varepsilon_{i}^0(1)` are all obtained during the trial KI calculation ``calc_alpha/ki``. Together with the value of our guess for the screening parameters (:math:`\alpha^0_i = 0.6`), this is sufficient to update our guess for :math:`\alpha_i` (see the :ref:`theory section <theory_dscf>` for details).

The procedure for empty orbitals is slightly different, as we can see when it comes to orbital 10:

----

.. include:: ../_static/tutorials/tutorial_1/md_excerpts/ozone_alpha_10.md
  :parser: myst_parser.sphinx_

----

where now we must call Quantum ESPRESSO several times in order to obtain :math:`E_i(N+1)`.

.. collapse:: Click here for detailed descriptions of each calculation

  ``01-dft_n+1_dummy`` and ``02-pz_print``
    preliminary calculations that generate files required by the subsequent constrained DFT calculation
  
  ``03-dft_n+1``
    a :math:`N+1`-electron PBE calculation where we freeze the 10\ :sup:`th` orbital, fill it, and allow the rest of the density to relax. This yields :math:`E_i(N+1)`

At the end of this section we can see a couple of tables:

----

.. include:: ../_static/tutorials/tutorial_1/md_excerpts/ozone_tables.md
  :parser: myst_parser.sphinx_

----

The first table lists the screening parameters :math:`\alpha_i` that we obtained -- we can see from row 0 we started with a guess of :math:`\alpha_i = 0.6` for every orbital `i`, and row 1 shows the alpha values.

The second table lists :math:`\Delta E_i - \lambda_{ii}^\alpha`. This is a measure of how well converged the alpha values are: if this value is close to zero, then the alpha values are well-converged. Note that the values listed above correspond to our starting guess of :math:`\alpha_i = 0.6`; this table does not show how well-converged the final alpha values are.

.. note::
  In order to see how well-converged our new screening parameters are, try increasing ``alpha_numsteps`` in the input file from ``1`` to ``2``. Can you make sense of the contents of the resulting tables?

The final calculation
^^^^^^^^^^^^^^^^^^^^^
Having determined the screening parameters, the final KI calculation is now run:

----

.. include:: ../_static/tutorials/tutorial_1/md_excerpts/ozone_final.md
  :parser: myst_parser.sphinx_

----

The files
---------

The markdown output provides an outline of the workflow that we ran, but for more details we can inspect the ``Quantum ESPRESSO`` input and output files themselves. These can be found in a nested directory structure that mirrors the outline we just followed. For example, the input and output files for the final KI calculation can be found in ``01-koopmans-dscf/03-ki_final/``. To get a sense of the entire directory structure, you can run the command

.. code-block:: bash

  tree -I TMP-CP 01-koopmans-dscf/

Extracting the ionization potential and electron affinity
---------------------------------------------------------
Let's now extract the KI ionization potential and electron affinity for ozone from our calculation.

The ionization potential (IP) corresponds to the negative of the energy of the HOMO (highest occupied molecular orbital). If you open ``01-koopmans-dscf/03-ki_final/ki_final.cpo`` and search near the bottom of the file you will see a section something like

.. code-block:: text

  ...
  HOMO Eigenvalue (eV)

  -12.5199

  LUMO Eigenvalue (eV)

  -1.8218

  Eigenvalues (eV), kp =   1 , spin =  1

  -40.1869  -32.9130  -24.2288  -19.6841  -19.4902  -19.2696  -13.6037  -12.7618  -12.5199

  Empty States Eigenvalues (eV), kp =   1 , spin =  1

  -1.8218

  Electronic Gap (eV) =    10.6981


  Eigenvalues (eV), kp =   1 , spin =  2

  -40.1869  -32.9130  -24.2288  -19.6841  -19.4902  -19.2696  -13.6037  -12.7618  -12.5199

  Empty States Eigenvalues (eV), kp =   1 , spin =  2

  -1.8218

  Electronic Gap (eV) =    10.6981
  ...

Very clearly we can see the HOMO eigenvalue of -12.52 eV. Thus we have a KI IP of 12.52 eV. This compares extremely well to the `experimental value <https://webbook.nist.gov/cgi/cbook.cgi?ID=C10028156&Mask=20#Ion-Energetics>`_ of ~ 12.5 eV, and is a marked improvement on the PBE result of 7.95 eV (which we can obtain from the ``HOMO Eigenvalue`` in ``01-initialization/01-dft_init_nspin2/04-dft_init_nspin2.cpo``).

Meanwhile, the electron affinity (EA) corresponds to the negative of the energy of the LUMO (lowest unoccupied molecular orbital). From the same files we can see that the KI EA is 1.82 eV (cf. ~ 2.1 eV experiment, 6.17 eV PBE)

.. tip::
  If you prefer working within ``python``, you need not write a script to parse the contents of ``.cpo`` files in order to extract the IP and EA. Instead, ``koopmans`` will have generated a python-readable file ``ozone.pkl`` containing all of the important calculation data.

  You can read these files like so:

  .. literalinclude:: ../../tutorials/tutorial_1/read.py

  Indeed, it is also possible to run the workflow from within ``python`` (rather than calling ``koopmans`` from the command line)

  .. literalinclude:: ../../tutorials/tutorial_1/run.py

  in which case you have immediate access to the workflow object ``wf`` rather than having to load a ``.pkl`` file.

