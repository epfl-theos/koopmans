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

  Valid options for ``units`` are ``angstrom``, ``crystal``, and ``alat``

``cell_parameters``
  describes the simulation cell. This can be specified explicitly

  .. literalinclude:: ../tutorials/tutorial_1/ozone.json
    :lines: 11-17
    :dedent:

  or using ``ibrav`` and ``celldms`` (following the conventions of Quantum ESPRESSO) e.g.

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
  :lines: 31-56
  :dedent:

Note that any keyword specified outside of a subblock (e.g. ``ecutwfc`` in the above example) is applied to all calculators for which it is recognized keyword.

.. note::
  Generally, you will need to provide very few keywords in this block. The vast majority of keywords for each calculation will be generated automatically. Notable exceptions to this are ``ecutwfc`` and (if relevant) the ``Wannier90`` projection

The pw subblock
~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``pw.x`` (see the `list of valid pw.x keywords <https://www.quantum-espresso.org/Doc/INPUT_PW.html>`_). Note that they should not be grouped by namelists as they are in a ``pw.x`` input file.

The w90 subblock
~~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``wannier90.x``, which are documented `here <https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf>`_. Because the occupied and empty manifolds are wannierized separately, you may want to use slightly different wannierization protocols for each. You can do this by placing keywords within ``occ`` and ``emp`` sub-dictionaries as in the above example. In this case both the occupied and empty manifolds will use ``sp3`` projections, but only the empty manifold will use the provided ``dis_froz_max`` etc.

Projections
  The projections can be specified as a list of dictionaries, where each dictionary corresponds to a single projection. The required entries for this dictionary are
  
    ``site``/``csite``/``fsite``
      an atom label/cartesian coordinate/fractional coordinate to be used as the projections' center. The three are mutually exclusive.
  
    ``ang_mtm``
      a string specifying the angular momentum states e.g. ``"l=2"``  
  
  The user can also optionally specify ``zaxis``, ``xaxis``, ``radial``, ``zona`` (see the `Wannier90 User Guide <https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf>`_ for details).
  
  Alternatively, the user can simply specify each projection as a single string using Wannier90 syntax e.g.
  
  .. code:: json
    
    "projections": [
      "f=0.25,0.25,0.25:sp3"
    ]

Projections blocks
  Often, the occupied or empty manifold consist of several well-separated blocks of bands. In this instance it is desirable to Wannierize each block separately, preventing the Wannierization procedure from mixing bands that are far apart in energy space.

  To do this, in place of ``projections`` we provide ``projections_blocks`` e.g.

  .. literalinclude:: ../tutorials/tutorial_3/zno.json
    :lines: 47-55
    :dedent:

  This corresponds to bulk ZnO where the first block of bands have Zn 3s character, the next Zn 3p, then O 2s, and finally Zn 3d hybridized with O 2p. Note that this syntax is practically identical to the ``projections`` field, except now we are providing a list of lists of dictionaries, not just a list of dictionaries. Note that the order of the projections blocks is important: they must run from lowest-energy to highest-energy. See :ref:`here for a tutorial on projections blocks <projections_blocks_explanation>`.

The pw2wannier subblock
~~~~~~~~~~~~~~~~~~~~~~~
This subblock contains ``pw2wannier.x`` keywords, in a single dictionary with no subdictionaries.

The kcp subblock
~~~~~~~~~~~~~~~~
This subblock contains keywords specific to ``kcp.x``, a modified version of ``cp.x`` for performing Koopmans calculations. In addition to `the keywords associated with cp.x <https://www.quantum-espresso.org/Doc/INPUT_CP.html>`_ there are several new keywords associated with the Koopmans implementation in ``kcp.x``. Non-experts will never need to change these.

The ui subblock
~~~~~~~~~~~~~~~
This subblock controls the unfolding and interpolation procedure for generating band structures and densities of states from Γ-only supercell calculations.

.. toctree::

  List of valid keywords <input_file/ui_keywords>


The plotting block
^^^^^^^^^^^^^^^^^^
This block can be used to customize the band structures and densities of states plots that ``koopmans`` generates.

.. toctree::

  List of valid keywords <input_file/plotting_keywords>
