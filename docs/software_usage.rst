.. _software_usage:

#####################
Software Usage
#####################

This section describes how software is used by the `WAVES`_ workflows, which use the
`SCons`_ automated build system.
It is assumed that the relevant software packages have been successfully installed
and linked in the :code:`config.yml` file using :code:`scons --config-software`.
See the :ref:`user_manual` section for description of basic `WAVES`_ and `SCons`_ commands.
For a thorough understanding of how to use the `WAVES`_ computation workflow tool, consider
exercising their provided tutorials.
A non-exhaustive discussion is provided to understand some aspects of the workflow and
software configuration.

The :code:`SConstruct` file is the main input that `SCons`_ reads to control workflows.
The software described in :ref:`software_installation` is linked using "BUILDERS":

.. literalinclude:: SConstruct.py
   :lines: 255-267

Individual workflows are created in SConscript files. These SConscript files are
stored in the :code:`model_package/workflows/` directory.
They are linked as simulation targets in the main :code:`SConstruct` file.
The current collection of workflows are as follows:

.. literalinclude:: SConstruct.py
   :lines: 294-314

Unless a user is adding another piece of software or creating a custom workflow,
these files do not need to be modified.

***********
General DNS
***********

In general, all DNS are setup and run in a similar way.
A mesh is generated along with a simulation input file.
The simulation is run and results are post-processed.
For each DNS, results are converted to the XDMF file format
required by the Micromorphic Filter.
Details vary for each DNS code.

**********
Abaqus FEM
**********

`WAVES`_ has a number of built in Abaqus tools including:

* :code:`waves.scons_extensions.journal()`
* :code:`waves.scons_extensions.solver()`
* :code:`waves.scons_extensions.abaqus_extract()`

`WAVES`_ includes many tutorials that describe how to use these Abaqus tools.

The :code:`waves.scons_extensions.journal()` tool uses Abaqus's built-in Python
interpretter to generate simulation input files.
Abaqus has a rhobust meshing capability which may be accessed through the journal
WAVES extension. Although meshes may also be generated using Cubit, the Abaqus workflows in
this repository use the built-in Abaqus meshing capabilities for simplicity.
All simulation inputs (parameters for geometry, mesh, material properties,
boundary and loading condition,
simulation duration and incrementation, etc.) are specified as arguments to an Abaqus
journal file.

After generating an Abaqus input deck, the :code:`waves.scons_extensions.solver()` tool
executes the simulation.

Abaqus simulation results are stored in a proprietary output database (ODB) format, which
is often hard to parse. However, the :code:`waves.scons_extensions.abaqus_extract()` tool
provides a method to extract an ODB file to an HDF5 file, which is easy to use for
post-processing to extract force vs displacement plots, stress vs strain plots, and much more.

Finally, the HDF5 file containing Abaqus results are converted to the XDMF file required for
the Micromorphic Filter using the :py:mod:`model_package.DNS_Abaqus.ODBextract_to_XDMF`
Python script.

To better understand the details associated with setting up an Abaqus simulation,
consider inspecting the :code:`Abaqus_elastic_cylinder.scons` SConscript located in
the :code:`model_package/workflows` directory.

.. _ratel_builder:

*********
Ratel FEM
*********

Meshes for simple Ratel DNS are generated using a Cubit Python script.

An options file is then generated using the
:py:mod:`model_package/DNS_Ratel/build_options_file` Python script.
All other simulation inputs (material properties, boundary and loading condition,
simulation duration and incrementation, etc.) are specified here.

To use the Ratel solver, a custom SCons Builder is setup as shown in the :code:`SConstruct`
file here:

.. literalinclude:: SConstruct.py
   :lines: 140-152

The options and mesh files are passsed to this Builder as input arguments.
Output arguments are also specified for the VTK monitor file, which will contain the simulation
results, and CSV force file, which will contain the reaction force histories.
The "-diagnostic_order 1" argument guarantees that the simulation output will be associated
with a first order finite element space. This is important because Ratel uses a nodal FEM, so
all outputs are provided at the nodes. The definition of nodal volumes (required for the
Micromorphic Filter) are only well-posed for linear elements.

Finally, the VTK file(s) containing the Ratel results are converted to the XDMF file required for
the Micromorphic Filter using the :py:mod:`model_package.DNS_Ratel.vtk_to_xdmf`
Python script.

To better understand the details associated with setting up a Ratel simulation,
consider inspecting the :code:`Ratel_elastic_cylinder.scons` SConscript located in
the :code:`model_package/workflows` directory.

********
GEOS MPM
********

..
   TODO: Describe how to build and link GEOS MPM

Coming soon!

*****
Cubit
*****

Cubit is used to generate meshes in the Exodus format through Cubit/Python scripts.
These meshes may be used as input to simple DNS, as macroscale (filtering domain)
meshes for the Micromorphic Filter, or as macroscale meshes for simulation in Tardigrade-MOOSE.
An example script is :py:mod:`model_package.Tardigrade_MOOSE.cylinder_from_bounds`.
Exodus meshes may be used for Tardigrade-MOOSE, but are converted to the XDMF format
required by the Micromorphic Filter using the `meshio`_ library.

*******************
Micromorphic Filter
*******************

Refer to section :ref:`workflow_homogenization` for discussion of how the Micromorphic Filter is used.

****************
Tardigrade-MOOSE
****************

Tardigrade-MOOSE may be used directly from the command line or through the included WAVES workflow.

.. _macroscale_command_line:

Using Tardigrade-MOOSE from Command Line
========================================

Similar to the ctest mentioned in :ref:`software_installation`, a user will need to
add the micromorphic element shared libraries to the LD_LIBRARY_PATH
(see :ref:`LD_PATH_NOTE`).
Specify the path to the compiled "tardigrade-opt" program and the input file.
A simulation may be run with multiple threads using the "--n-threads={n}" option.
For other parallelization options, see the relevant MOOSE documentation `MOOSE_parallel`.

   .. code-block:: console

      $ /path/to/tardigrade/build/tardigrade-opt -i {input_file}

Alternatively, the LD_LIBRARY_PATH may be specified on the command line when executing
Tardigrade-MOOSE.

   .. code-block:: console

      $ LD_LIBRARY_PATH=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp /path/to/tardigrade/build/tardigrade-opt -i {input_file}

One may set up a simple alias in their ~/.bashrc similar to this:
`alias run_tardigrade="LD_LIBRARY_PATH=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp /path/to/tardigrade/build/tardigrade-opt"`

which can be used as:

   .. code-block:: console

      $ run_tardigrade -i {input_file} --n-threads={#}`

.. _macroscale_WAVES:

Using Tardigrade-MOOSE with WAVES
=================================

It will be easiest to use Tardigrade-MOOSE through the WAVES workflow.
A custom SCons Builder is setup as shown in the :code:`SConstruct` file here:

.. literalinclude:: SConstruct.py
   :lines: 166-175

Setting up Tardigrade-MOOSE simulations
=======================================

Refer to section :ref:`workflow_macroscale` for discussion of how Tardigrade-MOOSE
is used.

*****************************
Micromorphic Calibration Tool
*****************************

The Micromorphic Calibration Tool provides the capability to evaluate micromorphic material
models pointwise. Calibration workflows evaluate stress measures for a given set of material
parameters and strains. These stresses are then compared to homogenized stresses
to inform the objective function for minimization to find the "best" micromorphic parameters.
See the :py:mod:`model_package.Calibrate.calibrate_element` script to see how this tool
is used for calibration tasks.

Refer to section :ref:`workflow_calibration` for discussion of how the micromorphic
calibration tool is used.

***************************************
Micromorphic Linear Elastic Constraints
***************************************

In order to ensure that the Helmholtz free energy function is positive definite,
the micromorphic linear elastic parameters must satisfy the contrainst discussed
in :ref:`linear_elastic_constraints`.

Calibration tasks evaluate these constraints by importing the
:code:`linear_elastic_parameter_constraint_equations` function from the
:code:`linear_elastic_parameter_constraint_equations.py` script.
The constrains are evaluated for a trial set of micromorphic linear elastic parameters
when evaluating the objective function in :py:mod:`model_package.Calibrate.calibrate_element`.
If any of the constraint evaluations are below 0, then the objective function
returns an error of infinity and that trial set of parameters is rejected.

Although not used by the WAVES workflows, the
:py:mod:`model_package.Calibrate.return_minimum_smith_conditions` utility script
may be used to evaluate these constraints for a set of parameters.