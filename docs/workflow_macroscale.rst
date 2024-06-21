.. _workflow_macroscale:

#####################
Macroscale Simulation
#####################

Tardigrade-MOOSE is a micromorphic, finite element tool built as a
custom application in
multiphysics object-oriented simulation environment (`MOOSE`_).

Tardigrade-MOOSE simulations require a mesh and an input file.
Currently, all macroscale simulations use simple cylindrical meshes
which may be created using the
:py:mod:`model_package.Tardigrade_MOOSE.cylinder_from_bounds` script.

The Tardigrade-MOOSE input deck may be created using the
:py:mod:`model_package.Tardigrade_MOOSE.build_Tardigrade_input_deck` script,
which will create a simple, quasi-static simulation with displacement
applied to the top surface in the z-direction with the bottom surface
fixed in the z-direction.
This script accepts several command line arguments including:

* the name of the Tardigrade-MOOSE input (.i) file to be generated,
* the name of the mesh file,
* the name(s) material input YAML file(s), either:

  * the path to a single material file to be assigned to the entire mesh,
  * the paths to multiple material files to be assigned to individual elements,

* an optional calibration map file which specifies which material parameter
  file will be associated with each individual macroscale element,
* the boundary conditions type, either:

  * "slip" which allows for the top and bottom surfaces to expand laterally, or
  * "clamp" which prevents the top and bottom surfaces from expanding laterally,

* the displacement to apply to the top surface, and
* the simulation duration.

See the API documentation and source code for :py:mod:`model_package.Tardigrade_MOOSE.build_Tardigrade_input_deck`
for futher details. A similar script,
:py:mod:`model_package.Tardigrade_MOOSE.build_dynamic Tardigrade_input_deck`
will construct an input file for dynamic simulations.

The input file may then be run in Tardigrade-MOOSE as described in the
:ref:`macroscale_command_line` or :ref:`macroscale_WAVES` sections.

Future efforts will consider expanding Tardigrade-MOOSE simulation capabilities
to include other loading conditions, more precise control over time
incrementation parameters, other material models, application of
micro-displacement degrees of freedom, frictional boundary conditions,
other types of meshes, and much more.

Most upscaling studies assume that the macroscale mesh for Tardigrade-MOOSE
is identical to the filtering domains used for the Micromorphic Filter
and subsequent material calibration (referred to as "direct" or "1 to 1"
upscaling).
Some studies explore other configurations in which a homogeneous macroscale
is considered with a material calibration informed from a statistical
sample of the calibration results.
Future efforts will explore heterogeneous macroscale simulations in which
material properties vary spatially by sampling from a distribution of
possible parameters.

To better understand the details associated with setting up a Tardigrade-MOOSE simulation,
consider inspecting the :code:`tardigrade_moose.scons` SConscript located in
the :code:`model_package/workflows` directory.