The input file
==============

``koopmans`` takes a single JSON file as an input. Here is an example, taken from :ref:`Tutorial 1 <tutorial_1>`:

.. literalinclude:: ../tutorials/tutorial_1/ozone.json

As you can see, the file is divided into three blocks: :ref:`workflow <input_file:The workflow block>`, :ref:`atoms <input_file:The atoms block>`, and :ref:`calculator_parameters <input_file:The calculator_parameters block>`. Other valid blocks are :ref:`kpoints <input_file:The kpoints block>`, :ref:`pseudopotentials <input_file:The pseudopotentials block>`, and :ref:`plotting <input_file:The plotting block>`. These are all explained below.

The workflow block
^^^^^^^^^^^^^^^^^^

The ``workflow`` block contains variables that define the details of the workflow that we are going to run.

.. toctree::

  List of valid keywords <input_file/workflow_keywords>

The atoms block
^^^^^^^^^^^^^^^

The ``atoms`` block contains details about the atomic positions and the simulation cell. It contains two subdictionaries

``atomic_positions``
  contains the atomic positions e.g.

  .. literalinclude:: ../tutorials/tutorial_1/ozone.json
    :lines: 18-25
    :dedent:

  Valid options for ``units`` are ``alat``, ``angstrom``, ``bohr``, and ``crystal``

``cell_parameters``
  describes the simulation cell. This can be specified explicitly

  .. literalinclude:: ../tutorials/tutorial_1/ozone.json
    :lines: 11-17
    :dedent:

  (with units ``angstrom``, ``bohr``, or ``alat``), or using ``ibrav`` and ``celldms`` following the conventions of Quantum ESPRESSO e.g.

  .. literalinclude:: ../tutorials/tutorial_2/si.json
    :lines: 15-19
    :dedent:

The kpoints block
^^^^^^^^^^^^^^^^^

The ``k_points`` block specifies the k-point sampling e.g.

.. literalinclude:: ../tutorials/tutorial_2/si.json
  :lines: 26-30
  :dedent:

There are five possible entries in this block

``grid``
  a list of three integers specifying the shape of the regular grid of k-points
``offset``
  a list of three integers, either ``0`` or ``1``. If ``1``, the regular k-point grid is offset by half a grid step in that dimension
``path``
  the path to be used in band structure plots, specified as a string of characters corresponding to special points of the Bravais lattice
``density``
  the number k-points per inverse Angstrom along the path
``gamma_only``
  set to ``True`` if the calculation is only sampling the gamma point


The pseudopotentials block
^^^^^^^^^^^^^^^^^^^^^^^^^^

``koopmans`` ships with several pre-installed pseudopotential libraries. To use these, select a ``pseudo_library`` in the ``workflow`` block. Valid options include ``pseudo_dojo_standard/stringent`` and ``sg15``.

Alternatively, you can provide your own pseudopotentials via the ``pseudopotentials`` block, where we specify the filenames of the pseudopotentials for each element e.g.

.. code-block:: json

  "pseudopotentials": {"O": "O.upf", "H": "H.upf"}

The directory that these pseudopotentials are contained in should be provided via the ``pseudo_directory`` keyword in the ``workflow`` block. See :ref:`here for more details on setting up pseudopotentials <running:Pseudopotentials>`.

.. warning::
  Koopmans currently only works with norm-conserving pseudopotentials


The calculator_parameters block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``calculator_parameters`` block can be used to specify code-specific codes e.g.

.. literalinclude:: ../tutorials/tutorial_2/si.json
  :lines: 31-48
  :dedent:

Note that any keyword specified outside of a subblock (e.g. ``ecutwfc`` in the above example) is applied to all calculators for which it is recognized keyword.

.. note::
  Generally, you will need to provide very few keywords in this block. The vast majority of keywords for each calculation will be generated automatically. Notable exceptions to this are ``ecutwfc`` and (if relevant) the ``Wannier90`` projection

The pw subblock
~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``pw.x`` (see the `list of valid pw.x keywords <https://www.quantum-espresso.org/Doc/INPUT_PW.html>`_). Note that they should not be grouped by namelists as they are in a ``pw.x`` input file.

The w90 subblock
~~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``wannier90.x``, which are documented `here <https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf>`_. The one keyword for which the syntax differs is the ``projections`` block, via which the user specifies the projections used to initialize the Wannier functions.

An individual projection can be specified as either a dictionary or a string.

As a dictionary
  If specifying a projection via a dictionary, the required entries for this dictionary are

  ``site``/``csite``/``fsite``
    an atom label/cartesian coordinate/fractional coordinate to be used as the projections' center. The three are mutually exclusive.

  ``ang_mtm``
    a string specifying the angular momentum states e.g. ``"l=2"``  

  The user can also optionally specify ``zaxis``, ``xaxis``, ``radial``, ``zona`` (see the `Wannier90 User Guide <https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf>`_ for details).

As a string
  If specifying a projection via a string, this string must follow the ``Wannier90`` syntax e.g. ``"f=0.25,0.25,0.25:sp3"``

These individual projections (either as dictionaries or as strings) must be provided to ``projections`` within a list of lists. This is because for Koopmans calculations, we want to perform the Wannierization in quite a particular way

  - the occupied and empty manifolds must be wannierized separately.
  - the occupied or empty manifold can consist of several well-separated blocks of bands. In this instance it is desirable to Wannierize each block separately, preventing the Wannierization procedure from mixing bands that are far apart in energy space.

We can achieve both of the above via the list-of-lists syntax. Consider the following example for the wannierization of bulk ZnO

.. literalinclude:: ../tutorials/tutorial_3/zno.json
  :lines: 47-55
  :dedent:

In ZnO, the bands form several distinct blocks. The first block of occupied bands have Zn 3s character, the next Zn 3p, then O 2s, and finally Zn 3d hybridized with O 2p. The first empty bands have Zn 4s character. You can see this reflected in the way the projections have been specified. If we were to run the workflow with this configuration, it will run five separate Wannierizations, one for each sub-list.

This means that

  - the occupied and empty manifolds will be wannierized separately, because the cumulative number of projections in the first four blocks is commensurate with the number of occupied bands in our system
  - we prevent mixing bands that belong to different sub-lists
  
See :ref:`here for a more detailed tutorial on projections <projections_blocks_explanation>`.

.. note::
  The order of the projections blocks is important: they must run from lowest-energy to highest-energy.

.. note::
  If disentanglement keywords such as ``dis_win_max`` are provided, these will only be used during the Wannierization of the final block of projections 

.. note::
  If running a spin-polarized calculation, you need to provide separately the spin-up and spin-down projections. You can specify this by splitting the ``w90`` block into an ``up`` and ``down`` subblocks, each containing the spin dependent ``projections``. See :ref:`here for an example of a spin-polarized set-up <tutorial_6>` 

The pw2wannier subblock
~~~~~~~~~~~~~~~~~~~~~~~
This subblock contains ``pw2wannier90.x`` keywords, in a single dictionary with no subdictionaries.

The kcp subblock
~~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``kcp.x``, a modified version of ``cp.x`` for performing Koopmans calculations. In addition to `the keywords associated with cp.x <https://www.quantum-espresso.org/Doc/INPUT_CP.html>`_ there are several new keywords associated with the Koopmans implementation in ``kcp.x``. Non-experts will never need to change these.


The kcw subblock
~~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``kcw.x`` (see the `list of valid kcw.x keywords <https://www.quantum-espresso.org/Doc/INPUT_KCW.html>`_). Non-experts will never need to change these keywords.


The ui subblock
~~~~~~~~~~~~~~~
This subblock controls the unfolding and interpolation procedure for generating band structures and densities of states from Γ-only supercell calculations.

.. toctree::

  List of valid keywords <input_file/ui_keywords>

The convergence block
^^^^^^^^^^^^^^^^^^^^^
This block can be used to customize a convergence calculation.

.. toctree::

  List of valid keywords <input_file/convergence_keywords>

See :ref:`here for a more detailed tutorial on performing convergence calculations <tutorials/tutorial_4:Tutorial 4: running convergence tests>`.

The plotting block
^^^^^^^^^^^^^^^^^^
This block can be used to customize the band structures and densities of states plots that ``koopmans`` generates.

.. toctree::

  List of valid keywords <input_file/plotting_keywords>

The ml block
^^^^^^^^^^^^

.. warning::
    This feature is experimental

This block controls the machine-learning of screening parameters.

.. toctree::

  List of valid keywords <input_file/ml_keywords>
