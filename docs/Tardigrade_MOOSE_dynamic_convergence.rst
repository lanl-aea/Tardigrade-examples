.. _tardigrade_moose_dynamic_convergence:

***************************************
Dynamic Convergence of Tardigrade-MOOSE
***************************************

Similar to the quasi-static convergence study presented in :ref:`tardigrade_moose_convergence`,
a simple study is conducted to investigate the convergence behavior of
Tardigrade-MOOSE dynamic implicit simulations.
See :ref:`Dynamic_verification` for description of the dynamically loaded
bar problem and solution by Meirovitch.
A cylinder with 5 mm diameter and 5 mm height is loaded with a force of 25 N.
The modulus of elasticity is 250 MPa and the density is 2000 kg/m^3.
To agree with the Meirovitch solution, a Poisson ratio of 0.0 is selected
to prevent lateral expansion/contraction to satisfy the 1D assumption.
For this geometry and material properties, the first natural frequency of vibration
is calculated to be 17.678 kHz. To simulate 3 cycles of vibration, the simulation
duration is approximately 0.170 millseconds.
It was found that 240 fixed timesteps produces satisfactory results
without requiring prohibitively long run times.
For the expected sawtooth displacement profile with 3 periods of vibration,
selecting timesteps in multiples of 6 is necessary to capture the maximum
and minimum displacements.

Four meshes are considered for simulation in Tardigrade-MOOSE with
192, 960, 2520, and 7680 elements.
See the :code:`Tardigrade_dynamic_convergence` Sconscript for details of how this
workflow is setup, as well as the :code:`dynamic_elastic_cylinder` dictionary in the
:code:`model_package/Tardigrade_MOOSE/simulation_variables_nominal.py` file
where input parameters are defined.

The mesh convergence analysis is setup as a parametric study. After simulations
are initialized, exectuted, and post-processed, the displacement versus time
results are collected for each case.

The analysis may be exectuted using the following command:

   .. code-block:: bash

      $ scons Tardigrade_dynamic_convergence

The :code:`--jobs=4` argument may be included to run each simulation concurrently.

The final force versus displacement plot is shown in
Figure :numref:`{number} <Tardigrade_dynamic_convergence_all_force_displacements>`.

.. figure:: Tardigrade_dynamic_convergence_all_force_displacements.png
   :name: Tardigrade_dynamic_convergence_all_force_displacements
   :align: center
   :width: 70%

   Displacement vs time results of dynamic Tardigrade-MOOSE convergence study

There is decent agreement between the analytical solution and Tardigrade-MOOSE
simulations. Using 60 and 120 timesteps provided poor results. It is expected
that using more timesteps (such as 480 or 600) would produce better agreement
at the expense of longer run times.