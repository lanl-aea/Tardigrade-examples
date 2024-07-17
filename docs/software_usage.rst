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
   :lines: 327-340

Individual workflows are created in SConscript files. These SConscript files are
stored in the :code:`model_package/workflows/` directory.
They are linked as simulation targets in the main :code:`SConstruct` file.
The current collection of workflows are as follows:

.. literalinclude:: SConstruct.py
   :lines: 366-388

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

Meshes for simple Ratel DNS are either generated using a Cubit/Python script
or pulled from the :code:`model_package/meshes` directory.

An options file is then generated using the
:py:mod:`model_package.DNS_Ratel.build_options_file` Python script.
All other simulation inputs (material properties, boundary and loading condition,
simulation duration and incrementation, etc.) are specified here.

It is assumed that Ratel simulations will be executed through WAVES workflows.
Refer to the `RATEL`_ documentation for more details.

To run Ratel simulations, a custom SCons Builder is setup as shown in the :code:`SConstruct`
file below. This is the basic builder for serial simulations. To use Ratel for parallel
simulations, see the discussion in :ref:`serial_vs_parallel`.

.. literalinclude:: SConstruct.py
   :lines: 160-167

The options and mesh files are passsed to this Builder as input arguments.
Output arguments are also specified for the VTK monitor file, which will contain the simulation
results, and CSV force file, which will contain the reaction force histories.
The "-diagnostic_order 1" argument guarantees that the simulation output will be associated
with a first order finite element space. This is important because Ratel uses a nodal FEM, so
all outputs are provided at the nodes. The definition of nodal volumes (required for the
Micromorphic Filter) is only well-posed for linear elements.

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
An example script is :py:mod:`model_package.Tardigrade_MOOSE.cylinder_from_bounds`.
These meshes may be used as input to simple DNS, as macroscale (filtering domain)
meshes for the Micromorphic Filter, or as macroscale meshes for simulation in Tardigrade-MOOSE.
The Micromorphic Filter uses the XDMF file format, so Exodus meshes are convert using the `meshio`_ library.

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
For other parallelization options, see the relevant MOOSE documentation `MOOSE_parallel`_.

   .. code-block:: console

      $ /path/to/tardigrade/build/tardigrade-opt -i {input_file}

Alternatively, the LD_LIBRARY_PATH may be specified on the command line when executing
Tardigrade-MOOSE.

   .. code-block:: console

      $ LD_LIBRARY_PATH=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp /path/to/tardigrade/build/tardigrade-opt -i {input_file}

One may set up a simple alias in their ~/.bashrc similar to:
:code:`alias run_tardigrade="LD_LIBRARY_PATH=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp /path/to/tardigrade/build/tardigrade-opt"`

which can be used as:

   .. code-block:: console

      $ run_tardigrade -i {input_file} --n-threads={#}`

.. _macroscale_WAVES:

Using Tardigrade-MOOSE with WAVES
=================================

It will be easiest to use Tardigrade-MOOSE through WAVES workflows.
A custom SCons Builder is setup as shown in the :code:`SConstruct` file below.
This is the basic builder for serial simulations with multi-threading.
To use Tardigrade-MOOSE for parallel simulations, see the discussion in
:ref:`serial_vs_parallel`.

.. literalinclude:: SConstruct.py
   :lines: 214-221

Setting up Tardigrade-MOOSE simulations
=======================================

Refer to section :ref:`workflow_macroscale` for discussion of how Tardigrade-MOOSE
is used.

*****************************
Micromorphic Calibration Tool
*****************************

The Micromorphic Calibration Tool provides the capability to evaluate micromorphic material
models pointwise. Calibration workflows evaluate stress measures for a given set of material
parameters and deformations. These stresses are then compared to homogenized stresses
to inform the objective function for minimization to find the "best" micromorphic parameters.
See the :py:mod:`model_package.Calibrate.calibrate_element` script to see how this tool
is used for calibration tasks.

Refer to section :ref:`workflow_calibration` for discussion of how the micromorphic
calibration tool is used.

***************************************
Micromorphic Linear Elastic Constraints
***************************************

In order to ensure that the Helmholtz free energy function is positive definite,
the micromorphic linear elastic parameters must satisfy the contraints discussed
in :ref:`linear_elastic_constraints`.

Calibration tasks evaluate these constraints by importing the
:code:`linear_elastic_parameter_constraint_equations` function from the
:code:`linear_elastic_parameter_constraint_equations.py` script.
The constrains are evaluated for a trial set of micromorphic linear elastic parameters
when evaluating the objective function in :py:mod:`model_package.Calibrate.calibrate_element`.
If any of the constraint evaluations are below 0, then the objective function
returns an error of infinity and that trial set of parameters is rejected.

Although not used by the WAVES workflows, the
:py:mod:`model_package.Calibrate.return_minimum_smith_constraint` utility script
may be used to evaluate these constraints for a set of parameters.

.. _serial_vs_parallel:

*******************
Serial vs. Parallel
*******************

Workflows in this repository chain together a considerable number of
computational utilities. Some workflow stages, especially pre/post-processing scripts,
run quickly on a single process. Other simulation stages may require parallelization.
By default, workflows will be executed in serial, however, some capabilities
to parallelize workflows have been implemented.

.. warning::

   Workflow parallelization techniques are under active deverlopment, but no general
   solution has been found. Different workflow stages will have different computational
   requirements, so a user is encouraged to carefully plan based on the relevant
   analysis for a specific system!

Homogenization activities via the Micromorphic Filter may be parallelized
using Scipy's "workers" argument. This value of this argument is currently specified in the
Python dictionary by the "filter_parallel" key for individual analyses.
A value of 8 workers seems to work well for most cases.

Ratel and Tardigrade-MOOSE may be parallelized using MPI (see :ref:`mpi`).
If a path to :code:`mpiexec` has been added to the :code:`config_software.yml`
file, then simulations will be executed for a specified number of cpus.
For some workflows, the number of cpus is specified using the
:code:`--solve-cpus=<n>` flag. For workflows such as `Tardigrade_convergence`,
the number of cpus for various meshes is hard coded in the parameteric study
parameter definition. Abaqus DNS are similarly parallelizable using the
:code:`--solve-cpus=<n>` flag without the need to configure MPI.

Additionally, Ratel and Tardigrade-MOOSE may be parallelized on HPC systems
with SLURM.

Note that the "BUILDER" definition (see top of this page) for Ratel and Tardigrade are
:code:`"RatelSolver": ratel_builder_select()` and
:code:`"TardigradeSolver": tardigrade_builder_select()`, respectively.
**Focusing on Tardigrade-MOOSE**, the :code:`tardigrade_builder_select` function provides
a method to configure the Tardigrade-MOOSE solver depending on the available
resources, see here:

.. literalinclude:: SConstruct.py
   :lines: 241-247

**Auto-batch**

If SCons detects that the user is on an HPC system with SLURM, the
:code:`env['sbatch']` environment variable is set.
If the user also specifies the :code:`--use-sbatch` flag, then the :code:`Tardigrade_solver`
will be set to the :code:`tardigrade_solver_sbatch` BUILDER.
This builder will automatically submit SLURM jobs.
This means that a user may run a workflow from a login node and Tardigrade-MOOSE simulations
will be automatically submitted to a compute node with :code:`--solve-cpus=<n>` cores
without requiring the user to write any SLURM scripts.
The workflow will wait for the job to finish
and then continue with the rest of the analysis.

Unfortunately, the :code:`env['sbatch']` configuration may cause workflows
to run slowly as individual workflow stages must wait through a SLURM queue.
It may be more efficient to submit an entire workflow execution directly
through a SLURM script as described in :ref:`WAVES_with_SLURM`.

**MPI**

If the path to :code:`mpiexec` is provided in the :code:`config_software.yml` file,
the :code:`Tardigrade_solver` will be set to the :code:`tardigrade_solver_mpi` BUILDER.
Tardigrade-MOOSE simulations will be parallelized using MPI on :code:`--solve-cpus=<n>` cores.
This configuration is primarily intended for workflows being submitted on standalone machines
without a job scheduler.

**Serial**

Finally, if a user does not (1) provide the :code:`--use-sbatch` flag with
:code:`env['sbatch']` detected, or (2) specify the "mpi" entry in
:code:`config_software.yml` file, then jobs will be run in serial using the
basic :code:`tardigrade_solver` BUILDER.

Identical "auto-batch", "MPI", and "Serial" configurations are provided for Ratel simulations.

.. _WAVES_with_SLURM:

**********************
Using WAVES with SLURM
**********************

As mentioned, several operations of Tardigrade-examples may be computationally
expensive. These operations may include running large direct numerical simulations,
running the Micromorphic Filter for large DNS on multiple cores,
running many calibrations simultaneously, or running large Tardigrade-MOOSE simulations.
Some WAVES configuration options have been discussed in the previous section, but
some situations using HPC systems may warrant more manual control using SLURM scripts
for specific workflow stages.
A template script is provided in the root repository with the name
:code:`sbatch_template.sh`, as shown below:

.. literalinclude:: sbatch_template.sh
   :linenos:

A user may modify this `SLURM`_ script as necessary for the intended system,
especially the "partition" and "qos" SBATCH directives.
The "node", "ntasks", "time", "job-name", and "output" directives should be modified
for the specific workflow stage to be run.
The "LD_LIBRARY_PATH" should be changed according to the discussion provided
in :ref:`LD_PATH_NOTE` especially for tasks involving the micromorphic
calibration tool or Tardigrade-MOOSE.
Similarly, the :code:`/path/to/conda/bin/activate <tardigrade-examples-env>`
line should be modified to point to the correct conda environment and
associated directory where conda is installed.

The specific SCons task to be executed is passed in from the command line.
For example, the :code:`Ratel_elastic_cylinder` workflow can be executed
with:

   .. code-block:: console

      $ sbatch sbatch_template.sh Ratel_elastic_cylinder

A more complex, custom command may be specified by wrapping a command in quotations.
The following example would run the calibration Stage of the :code:`Ratel_elastic_cylinder_multi_domain`
workflow on 28 CPUs (notes that the "ntasks" SLURM directive would also need to be set to 28):

   .. code-block:: console

      $ sbatch sbatch_template.sh "Ratel_elastic_cylinder_multi_domain --calibrate --keep-going --jobs=28"

This script is simply a template and a user is encouraged to adapt it as needed.