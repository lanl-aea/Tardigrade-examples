.. _clamped_elastic_cylinder:

################################################
Clamped elastic cylinder - Quasi-static Implicit
################################################

Additional worklows are available to investigate the effect
of boundary conditions on upscaling. As it pertains to the discussion
for the heterogeneous Ratel I41.02 study (see :ref:`Ratel_I41_02_elastic`),
"clamped" boundary conditions introduce concentrated stresses in a DNS
at the top and bottom boundaries since these surfaces are not allowed
to move laterally. These studies attempt to isolate the boundary condition
effects for homogeneous cylinders. The only difference between these
studies and those described in :ref:`Ratel_elastic_cylinder` and
:ref:`Abaqus_elastic_cylinder` are the boundary conditions applied
in DNS and macroscale simulations, as well as the handling of
calibration results.

The Ratel study may be executed using SCons with the
:code:`Ratel_elastic_cylinder_clamped` or
:code:`Ratel_elastic_cylinder_clamepd_multi_domain` SConscripts.
For example:

   .. code-block:: bash

      $ scons Ratel_elastic_cylinder_clamped
      $ scons Ratel_elastic_cylinder_clamped --filter
      $ scons Ratel_elastic_cylinder_clamped --calibrate
      $ scons Ratel_elastic_cylinder_clamped --macro
      $ scons Ratel_elastic_cylinder_clamped --summary

Similarly, the Abaqus study may be executed using SCons with the
:code:`Abaqus_elastic_cylinder_clamped` or
:code:`Abaqus_elastic_cylinder_clamepd_multi_domain` SConscripts.
For example:

   .. code-block:: bash

      $ scons Abaqus_elastic_cylinder_clamped
      $ scons Abaqus_elastic_cylinder_clamped --filter
      $ scons Abaqus_elastic_cylinder_clamped --calibrate
      $ scons Abaqus_elastic_cylinder_clamped --macro
      $ scons Abaqus_elastic_cylinder_clamped --summary

Documentation currently has not been created for these studies.

Future work will also consider the more experimentally realistic
treatment of boundary conditions using frictional contact for loading
platens.