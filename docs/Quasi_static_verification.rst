.. _Quasi_static_verification:

###################################
Introduction & Analytical Solutions
###################################

A simple study is considered to verify that the micromorphic upscaling workflow is
functioning properly. The goal is to conduct a linear elastic simulation in a DNS code,
homogenize the DNS results via the Micromorphic Filter, calibrate the simplified "two
parameter" micromorphic linear elasticity (i.e. Saint Venant-Kirchhoff) model which should
recover the input material properties of the DNS, and run a macroscale simulation in
Tardigrade-MOOSE.

A uni-axial stress state is considered for a cylindrical geometry with a height, :math:`h_0`,
and diameter, :math:`d_0`, of 5 mm.
A coordinate system is assumed in which the x-, y-, and z-directions correspond to
unit vectors identified with subscripts 1, 2, and 3, respectively. The axis of the cylinder
is aligned with the z-direction. A compressive displacement, :math:`\Delta h`,
of -0.05 mm is applied to the
top surface of the cylinder corresponding to a nominal strain, :math:`\varepsilon_{33}`, of
-1%. The bottom surface is restricted from motion in the z-direction. For this study, nominal
elastic properties for the FK-800 binder are chosen with an elastic modulus, :math:`E^*`, of
165.0 MPa and a Poisson's ratio, :math:`\nu^*`, of 0.39
:cite:`allard_material_2022,3m_specialty_materials_fk-800_nodate,clements_kel-f_2007`.

.. include:: small_strain_solution.txt

.. include:: finite_strain_solution.txt