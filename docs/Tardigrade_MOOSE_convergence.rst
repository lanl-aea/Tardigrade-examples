.. _tardigrade_moose_convergence:

***************************************************
Convergence of Tardigrade-MOOSE for uniaxial stress
***************************************************

Tardigrade-MOOSE is a finite deformation code. The linear elastic micromorphic
constitutive model is geometrically nonlinear which may be considered an extension
of the Saint Venant-Kirchhoff model :cite:`ogden_non-linear_1997`.
As such, simulation results should be compared to analytical solutions in the
finite deformation regime.
A simple study is conducted that investigates the convergence behavior of
Tardigrade-MOOSE for a set of classical linear elastic material inputs with
similar loading and boundary conditions as described in the previous sections.
A cylinder with 5 mm diameter and 5 mm height is displaced in the z-direction
to produce a nominal strain of -1% with boundary conditions that result in a
uniaxial stress state.

Refer to :ref:`finite_strain_solution` for details of the analytical solution.
A user may find the :py:mod:`model_package.Tardigrade_MOOSE.finite_stVK_calculation`
script helpful for reproducing the relevant calculations.
It can be shown that the relevant stresses and forces for this problem are

.. math::

    S_{33} &= -2.48750 MPa

    \sigma_{33} &= -2.45286 MPa

    F_3 &= -48.3535 N.

Four meshes are considered for simulation in Tardigrade-MOOSE with
192, 960, 2520, and 7680 elements.
See the :code:`Tardigrade_convergence_study` Sconscript for details of how this
workflow is setup, as well as the :code:`elastic_cylinder` dictionary in the
:code:`model_package/Tardigrade_MOOSE/simulation_variables_nominal.py` file
where input parameters are defined.

The mesh convergence analysis is setup as a parametric study. After simulations
are initialized, exectuted, and post-processed, the force versus displacement
results for each case are collected.

The analysis may be exectuted using the following command:

   .. code-block:: bash

      $ scons Tardigrade_convergence

The :code:`--jobs=4` argument may be included to run each simulation concurrently.

The final force versus displacement plot is shown in
Figure :numref:`{number} <Tardigrade_convergence_all_force_displacements>`.

.. figure:: Tardigrade_convergence_all_force_displacements.jpeg
   :name: Tardigrade_convergence_all_force_displacements
   :align: center
   :width: 70%

   Absolute value of force versus displacement results of Tardigrade-MOOSE convergence study

The final force results are -47.1203, -48.0434, -48.1549, and -48.2759 N for the
meshes with 192, 960, 2520, and 7680 elements, respectively.
These results indicate that Tardigrade-MOOSE appears to converge to the
analytical force value of -48.3535 N.
See the :code:`docs/Tardigrade_convergence_all_force_displacements.csv`
file for the full force-displacement profiles.

.. note::

    Tardigrade-MOOSE does not currently provide stress quantities evaluated in the current configuration
    including Cauchy stress.