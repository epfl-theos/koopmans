.. _tutorial_5:

Tutorial 5: using machine learning to predict the screening parameters of water molecules
=========================================================================================
In this tutorial, we will train a machine learning model to predict the screening parameters of water molecules directly from their orbital densities. 
To create a trajectory with 20 different atomic configurations, we run a 
:download:`python script <../../tutorials/tutorial_5/perturb_positions.py>` that applies random noise to the unperturbed atomic positions of a water molecule. The resulting atomic positions are saved in a :download:`xyz-file <../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.json>` and visualized below

.. figure:: ../../tutorials/tutorial_5/snapshots.gif
   :width: 400
   :align: center
   :alt: The 20 snapshots generated with perturb_positions.py 

Our goal in this tutorial is to perform 
Koopmans calculations on each of these 20 snapshots using a machine learning model to predict the screening parameters.  


Running a convergence analysis
------------------------------

When running the machine learning workflow for the first time on a completely new system, it can be sensible to run a convergence analysis with respect to the number of training data.
If you are already confident in the model, you can skip this section.

The input file for the convergence analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The input file for running a convergence analysis for the machine learning model can be downloaded :download:`here <../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.json>`. 

First, we have to specify the corresponding task in the ``workflow`` block:

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.json
  :lines: 2-4
  :linenos:
  :emphasize-lines: 2

For this task, we don't provide the atomic ``"positions"`` directly to the input file since we don't want to perform a Koopmans calculation on a single snapshot but on many snapshots. Instead, we provide the xyz-file containing all the atomic positions of each snapshot that we would like to simulate

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.json
  :lines: 23-27
  :linenos:
  :emphasize-lines: 5

Finally, we have to provide a ``ml`` block with keywords specific to the machine learning model

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.json
  :lines: 12-22
  :linenos:

To predict the screening parameters from the orbital densities, we have to translate the orbital densities into input vectors for the machine learning model. To do so, we decompose the orbital densities into radial basis functions :math:`g_{nl}(r)` and angular basis functions :math:`Y_{ml}(\theta,\phi)`. 
This decomposition has the following four hyperparameters that we have to provide in the input file:

* :math:`n_{max}` determines the number of radial basis functions
* :math:`l_{max}` determines the number of angular basis functions
* :math:`r_{min}` determines the smallest cutoff radius for the radial basis functions
* :math:`r_{max}` determines the smallest cutoff radius for the radial basis functions


For the ``convergence_ml`` task, setting ``"number_of_training_snapshots": 10`` means that we will perform the convergence analysis with respect to 1,2,.., and 10 training snapshots and use the remaining snapshots (in this case snapshots 11 to 20) for testing. 

The ``"quantities_of_interest"`` is the list of parameters with respect to which we would like to perform the convergence analysis. In addition to performing it only with respect to the screening parameters ``"alphas"``, we could also perform it with respect to the eigenvalues if we had added ``"evs"`` to this list.
However, this would require an additional ``final calculation`` for each snapshot and therefore takes longer to run. 

In anticipation that the machine learning model will be most useful in extended systems (liquids or solids), we apply periodic boundary conditions and use maximally localized Wannier functions as our variational orbitals.

The output file for the convergence analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should see that the workflow first computes the screening parameters ab-initio for the last 10 snapshots.

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.out
  :lines: 18-20
  :language: text
  :lineno-start: 18

Next, snapshot 1 is added to the training data. 

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.out
  :lines: 733-734
  :language: text
  :lineno-start: 733

After having trained the machine learning model on the orbitals of the first snapshot we use the trained model to predict the screening parameters of the last 10 snapshots and compare our results to the results from the ab initio computation.

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.out
  :lines: 862-863
  :language: text
  :lineno-start: 862

Next, we add snapshot 2 to the training data. 

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5b/h2o_convergence_ml.out
  :lines: 1625-1626
  :language: text
  :lineno-start: 1625

Our model is now trained on the orbitals of 2 snapshots. We use this model again to predict the screening parameters of the last 10 snapshots and compare the results to the ab initio computation. 
We repeat this procedure until we have added all 10 snapshots to the training data. Then we can have a look at the convergence of the mean absolute error of the predicted screening parameters:

.. figure:: ../../tutorials/tutorial_5/tutorial_5b/spin_0_alphas_MAE_convergence.png
   :width: 1000
   :align: center
   :alt: The 20 snapshots generated with perturb_positions.py 

and (if we additionally had specified ``"evs"`` in the ``"quantities_of_interest"`` list) the convergence of the mean absolute error of the predicted orbital energies:

.. figure:: ../../tutorials/tutorial_5/tutorial_5b/spin_0_evs_MAE_convergence.png
   :width: 1000
   :align: center
   :alt: The 20 snapshots generated with perturb_positions.py 

We can see that we converged to a reasonable accuracy after about 5 training snapshots (which corresponds to 20 occupied and 10 empty orbitals).


Running a machine learning workflow
------------------------------------------

From the previous section (or from the results of similar systems) we know now, that 5 training snapshots (or 30 training orbitals) are sufficient to train our model for this system. 
Normally, we want to save computational time by computing the screening parameters only for the first 5 snapshots ab initio and not for all of them. This is the goal of this section.

The input file for the machine learning workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The corresponding :download:`input file <../../tutorials/tutorial_5/tutorial_5a/h2o_trajectory_ml.json>` differs from the previous input file only in the ``"task"`` keyword

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5a/h2o_trajectory_ml.json
  :lines: 2-4
  :linenos:
  :emphasize-lines: 2

and the ``"number_of_training_snapshots"``

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5a/h2o_trajectory_ml.json
  :lines: 12-22
  :linenos:
  :emphasize-lines: 8

The output file for the machine learning workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You should see that we compute the screening parameters of the first 5 snapshots ab initio and add the results to our training data

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5a/h2o_trajectory_ml.out
  :lines: 49-60
  :language: text
  :lineno-start: 49

Then we use the trained model to predict the screening parameters of the remaining snapshots 

.. literalinclude:: ../../tutorials/tutorial_5/tutorial_5a/h2o_trajectory_ml.out
  :lines: 680-683
  :language: text
  :lineno-start: 680