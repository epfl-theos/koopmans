Tutorial 4: calculating electron affinities for small anions
============================================================
.. tip:: To run this tutorial, you will need a version of ``pw.x`` with the `environ <https://environ.readthedocs.io/en/latest/>`_ patch installed.

The ``koopmans`` package can also calculate the PBE electron affinities of small molecules using the method of Nattino *et al.*. These anions are typically unbound (wrongly) by PBE, which means we cannot perform a standard ΔSCF calculation. Instead, the molecule is embedded within a cavity and finite difference calculations are performed with increasingly small values of :math:`\varepsilon_\infty`. See :cite:`Nattino2019` for a more detailed description.

Running these calculations is enabled with the ``environ_dscf`` task, and ``eps_cavity`` is a list of the trial values of :math:`\varepsilon_\infty` to use e.g.

.. literalinclude:: ../../tutorials/tutorial_4/o2_environ_dscf.json
  :lines: 2-5
  :linenos:

The full input file can be downloaded `here <https://raw.githubusercontent.com/elinscott/koopmans_docs/main/o2_environ_dscf.json>`_. When you run this calculation, the output will be as follows:

.. code-block:: text

  $ koopmans o2_environ_dscf.json
  PBE ΔSCF WORKFLOW

  Performing neutral calculations...
  Running neutral/8/pbe... done
  Running neutral/6/pbe... done
  Running neutral/4/pbe... done
  Running neutral/3.5/pbe... done
  Running neutral/3/pbe... done
  Running neutral/2.5/pbe... done
  Running neutral/2/pbe... done
  Running neutral/1/pbe... done

  Performing charged calculations...
  Running charged/8/pbe... done
  Running charged/6/pbe... done
  Running charged/4/pbe... done
  Running charged/3.5/pbe... done
  Running charged/3/pbe... done
  Running charged/2.5/pbe... done
  Running charged/2/pbe... failed to converge

  WORKFLOW COMPLETE

so we can see that for :math:`\varepsilon_\infty = 2` the anion became unstable, as expected. If we perform a quartic fit to the energies (following the example of Nattino *et al.*) we can extrapolate back to :math:`\varepsilon_\infty = 1` to obtain the electron affinity of 1.30 eV.

.. image:: ../../tutorials/tutorial_4/o2_dscf_ea_result.png
  :width: 800
  :alt: Quartic fit to embedded energies of O2 to calculate its vertical electron affinity 
  :align: center

.. warning::
  The `koopmans` implementation of this workflow differs from Nattino *et al.* in one major aspect: we use the same atomic positions for the anion as the neutral molecule. This means that we obtain *vertical* rather than *adiabatic* electron affinities. The reason for this choice is to be consistent with Koopmans spectral functionals, whose LUMO energies correspond to vertical electron affinities.
