.. _Ratel_I43_09_elastoplastic:

###########################################################
Ratel I43_09 elastoplastic cylinder - Quasi-static Implicit
###########################################################

This study investigates micromorphic upscaling of a heterogeneous
material composed of grains and a binder.
The DNS is an inelastic, implicit finite element
(FE) model assuming quasi-static conditions conducted in Ratel.
The results of this DNS are prepared for the Micromorphic Filter
for a variety of filtering domains.
Filter output is then used to calibrate an elastic micromorphic material
model which is implemented for simulation in Tardigrade-MOOSE.
Unique material calibrations are produced for each filtering domain
which are prescribed to their corresponding elements in the
macroscale simulations.

A variety of simulation variables are provided in the :code:`I43_damage_coarse_finetime`
dictionary in the :code:`model_package/DNS_Ratel/simulation_variables_nominal.py` file.
This dictionary is loaded into the workflow as the :code:`params` dictionary.

.. literalinclude:: DNS_Ratel_simulation_variables_nominal.py
   :lines: 176-219
   :linenos:

.. warning::

   The DNS files needed for this upscaling study can only be access by users
   with a CU Boulder identikey and access to the PetaLibrary!

DNS files are copied from the PetaLibrary into a local directory
"peta_data_copy" using the following command. This command usually only needs
to be used once. Files are copied using the secure copy protocal (SCP).
A user will be asked for their identikey, password, and two factor authentication
(2FA).

   .. code-block:: bash

      $ scons --peta-data-copy

The analysis is executed in individual stages by the
:code:`Ratel_I41_02_elastic_multi_domain` SConscript.
First, the existing
DNS results are processed into the required XDMF file format for
the Micromorphic Filter using the following command:

   .. code-block:: bash

      $ scons Ratel_I43_09_multi_domain

Next the homogenization is performed. Macroscale meshes with
1, 24, 48, 264, 768, and 2160 elements are considered for the default configuration.
The micromorphic filter is parallelized to run on 8 cpus. This value may be
modified by changing the value of "filter_parallel" in the :code:`I43_damage_coarse_finetime`
parameter dictionary.

   .. code-block:: bash

      $ scons Ratel_I43_09_multi_domain --filter

Calibration is then performed for each macroscale element (i.e. filtering domain).
This process can be rather expensive depending on what model is being calibrated,
however, WAVES provides a simple way to parallelize the independent calibration
processes using the :code:`--jobs=N` command line option.
Calibration may be
performed using a single process or multiple (10 in this example) with
either of the following two options.
Parameter sets 0 through 5 correspond to the macroscale meshes with 1, 24,
48, 264, 768, and 2160 elements.
One may choose to only calibrate parameter sets 0, 1, and 2 using the
:code:`--selected-parameter-sets` command line option.


   .. code-block:: bash

      $ scons Ratel_I43_09_multi_domain --calibrate
      $ scons Ratel_I43_09_multi_domain --calibrate --jobs=10
      $ scons Ratel_I43_09_multi_domain --calibrate --jobs=10 --selected-parameter-sets='0 1 2'

Once calibration is completed, Tardigrade-MOOSE simulations may be performed.
A Tardigrade-MOOSE simulation is performed for each case of filtering domains.
Macroscale simulations may be parallelized using the :code:`--solve-cpu=N`
command line option.
Macroscale simulations may be performed using a single process or multiple
(12 in the example) below.
There are three different types of macroscale simulations that may specified:

* The :code:`--macro` command line option will run simulations with clamped boundary conditions
* The :code:`--macro-platen` command line option will run simulations with loading platens
* The :code:`--macro-damage` command line option will run simulations using gradient-enhanced damage plasticity if Tardigrade-MOOSE is configured properly.

Finally, the results across filtering domains may be summarized into several
plots and csv files using the "summary" command:

   .. code-block:: bash

        $ scons Ratel_I41_02_elastic_multi_domain --summary

***************************
DNS Description and Results
***************************

Discussion coming soon

******************
Filter Preparation
******************

Discussion coming soon

**************
Filter Results
**************

Discussion coming soon

*******************************************
Micromorphic Constitutive Model Calibration
*******************************************

Discussion coming soon

*********************
Macroscale Simulation
*********************

Discussion coming soon
