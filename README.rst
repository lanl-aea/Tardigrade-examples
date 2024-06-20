.. target-start-do-not-remove

.. _Abaqus: https://www.3ds.com/products/simulia/abaqus
.. _Conda: https://docs.conda.io/en/latest/
.. _Conda installation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _Conda environment management: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _SCons: https://scons.org/
.. _SCons documentation: https://scons.org/documentation.html
.. _SCons manpage: https://scons.org/doc/production/HTML/scons-man.html
.. _WAVES: https://lanl.github.io/waves/index.html
.. _WAVES-EABM: https://github.com/lanl/waves/tree/main/modsim_template
.. _`Thomas Allard`: tea@lanl.gov

.. _PetaLibrary: https://www.colorado.edu/rc/resources/petalibrary
.. _ResearchComputing: https://www.colorado.edu/rc/
.. _meshio: https://github.com/nschloe/meshio
.. _MOOSE: https://mooseframework.inl.gov/index.html
.. _MOOSE_parallel: https://mooseframework.inl.gov/getting_started/examples_and_tutorials/tutorial01_app_development/step07_parallel.html
.. _PSAAP: https://micromorph.gitlab.io

.. target-end-do-not-remove

###################
Tardigrade-Examples
###################

.. inclusion-marker-do-not-remove

.. project-description-start-do-not-remove

Tardigrade-examples (LANL code O4735) is a repository of computational workflows that
exercise the Tardigrade software package.

The Tardigrade software package is an implementation of Eringen's
micromorphic continuum theory with capabilities to support multiscale
material modeling workflows. These capabailities include homogenization
through the Micromorphic Filter, calibration of micromorphic materials
models, and macroscale simulation in Tardigrade-MOOSE.

This repository investigates continuum upscaling of various direct numerical
simulations (DNS) conducted in Abaqus finite element (FE), Ratel FE, and
GEOS material point method (MPM) software. Verfication of the upscaling
workflow is first investigated by considering DNS of trivial stress states
for homogeneous materials, results of which indicate that classical
continuum behavior is recovered as expected. DNS of heterogeneous materials
are then considered.

***********
Description
***********

Tardigrade-examples includes upscaling analyses of a variety of DNS.
The `WAVES`_ computational workflow tool is used to automate these studies
in a reproducible manner primarily through a collection of Python scripts.

A full upscaling study considers four main analysis stages:

#. **DNS**
#. **Filter**
#. **Calibrate**
#. **Tardigrade-MOOSE**

The included **DNS** investigate simple compression of a homogenoeous cylinder.
Alternatively, externally performed DNS results may be upscaled that
consider more interesting materials and stress states. DNS results are
post-processed into the XDMF file format required by the Micromorphic
Filter. This repository currently includes Python scripts for XDMF conversion
of the Abaqus output database (ODB) and Ratel visualization toolkit (VTK) formats,
with GEOS MPM VTK coming soon.

The **Micromorphic Filter** homogenizes DNS results onto a macroscale mesh where
each element is considered a "filtering domain." Meshes are generated using
Cubit. For users without access to Cubit, a collection of example meshes are
included, but functionality is limited.
Upscaling workflows consider macroscale meshes of varying element size
using WAVES parametric study tools to analyze how homogenization changes
as the number of filtering domains increase.

The homogenized, micromorphic stress and deformation measures output by the
Micromorphic Filter are used to **calibrate** various material models. A unique
calibration for each filtering domain (macroscale element) is produced. For
verification studies of simple DNS, classical inear elasticity is calibrated,
while several forms of micromorphic linear elasticity are
calibrated for heterogeneous DNS. Calibration of elasto-plasticity will be
included in a future release of this repository.

Finally, material calibration results are applied to macroscale simulations in
**Tardigrade-MOOSE** with similar boundary and loading conditions as the original
DNS. Python scripts are included that collect and summarize results across
parametric studies for different analysis stages.

Documentation is included which details the upscaling methodology, steps for building
and connecting required software, relevant micromorphic theory, how the
Micromorphic Filter works, constitutive theories, and results of several
upscaling studies.

.. project-description-end-do-not-remove

#############
Documentation
#############

Documentation for this repository is included here:

* Production version (``main`` branch): WIP, provide link
* Development version (``dev`` branch): https://aea.re-pages.lanl.gov/models-and-simulations/tardigrade-examples

******************
Upscaling Workflow
******************

.. upscaling-workflow-description-start-do-not-remove

..
   TODO: decide if I'll use rst to copy over content from the README

.. upscaling-workflow-description-end-do-not-remove

*********************
Software Requirements 
*********************

.. software-requirements-description-start-do-not-remove

This repository uses the `WAVES`_ computational workflow tool (which utilizes the `SCons`_
automate build system) to configure and exectute upscaling workflows.
Workflows contained in this repository use a wide
array of software. As such, basic description is provided for how to
properly install and link required software.
 
Currently, the following DNS codes are considered:

* Abaqus implicit finite element (FE) and
* Ratel implicit FE.

Abaqus explicit FE and GEOS explicit material point method (MPM) DNS will be added in the future.

Additionally, Cubit is used for a number of meshing operations. For users without access to Cubit,
several example meshes are contained in :code:`model_package/meshes/`, however, functionality
of workflows will be limited.

Tardigrade software products include:

* Micromorphic Filter
* Tardigrade-MOOSE
* Micromorphic Calibration Tool
* Micromorphic Linear Elastic Constraints

The Tardigrade software must be properly configured to perform any upscaling
analyses. The DNS software, on the other hand, is not necessarily required.

..
   TODO: finish describing software requirements!

.. software-requirements-description-end-do-not-remove

Developers
==========

* `Thomas Allard`_

****************
Copyright Notice
****************

.. copyright-start-do-not-remove

Copyright (c) 2024, Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL),
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
copies to the public, perform publicly and display publicly, and to permit others to do so.

.. copyright-end-do-not-remove

********************
Activate Environment
********************

.. env-start-do-not-remove

AEA server environments
=======================

A shared `AEA Compute environment`_ is maintained on AEA servers. See the `AEA Compute environment`_ documentation for
the official use and activation instructions. A minimal activation description is included below for convenience.

1. Add the AEA modulefiles directory

   .. code-block::

      $ module use /projects/aea_compute/modulefiles

2. Load the shared environment modulefile

   .. code-block::

      $ module load tardigrade-examples-env

Local environments
==================

For users external to LANL systems, an environment to run workflows in this repository can be installed in a
`Conda`_ environment with the `Conda`_ package manager.
See the `Conda installation`_ and `Conda environment management`_ documentation
for more details about using `Conda`_.

1. Create the environment if it doesn't exist

   .. code-block::

      $ conda env create --name tardigrade-examples-env --file environment.yml

2. Activate the environment

   .. code-block::

      $ conda activate tardigrade-examples-env

.. env-end-do-not-remove

******************
SCons Build System
******************

.. build-start-do-not-remove

The `SCons`_ automated build system is used to execute workflows.
This section will discuss some common build operations. An abbreviated
options description can be displayed with ``scons -H``. For a full list of `SCons`_ command line options and target
build behavior, see the `SCons manpage`_. The `SCons manpage`_ is also installed with `Scons`_ in the environment and
can be opened from the command line as ``man scons`` in the `AEA Compute environment`_. In local environments, the
manpage may not be in the ``man`` program's search path, ``MANPATH``. You can find the manpage file and make them
available with something similar to any of the following, in increasing order of required background knowledge.

.. code-block::

   # Activate the environment
   conda activate tardigrade-examples-env

   # Find the scons manpage file
   $ find $CONDA_PREFIX -name scons.1
   /path/to/tardigrade-examples-env/scons.1

   # Open manpage directly
   $ man $CONDA_PREFIX/scons.1

   # Link SCons manpage to expected path and update MANPATH
   $ ln -s $CONDA_PREFIX/scons.1 $CONDA_PREFIX/man/man1/scons.1
   $ export MANPATH=$MANPATH:$CONDA_PREFIX/man
   $ man scons

- View project specific command line options

  .. code-block::

     $ scons -h
     ...

This project limits the default target list to the documentation with the `SCons`_ ``Default`` command. Simulation
targets must be specified directly on the command line. The `SCons`_ "all targets" character, ``.``, may also be
specified to build every target in the repository, including *all* simulation targets. Simulation targets may be
specified by output file name or by target alias, which is set to match the parent directory for the target
configuration, e.g. ``Abaqus_elastic_cylinder``.

- View the default targets and target aliases

  .. code-block::

     $ scons -h
     ...

- Build default targets

  .. code-block::

     $ scons

- Build *all* targets

  .. code-block::

     $ scons .

- Build a specific target

  .. code-block::

     $ scons <target name>

- Remove *all* build target artifacts

  .. code-block::

     $ scons . --clean

.. build-end-do-not-remove

*********************
PetaLibrary Data Copy
*********************

.. peta-start-do-not-remove

Several WAVES workflows upscale DNS run by others from the CU Boulder PSAAP project
and stored on the `PetaLibrary`_.
These DNS results may be copied using the following command:

  .. code-block::

     $ scons --peta-data-copy

A user will be asked for their identikey, password, and a dual authentication requrest
before the secure copy (SCP) transfers files.

.. note::
    This data may only be accessed for users with a Colorado `ResearchComputing`_ account
    with an allocation to the appropriate PSAAP user group.

.. peta-end-do-not-remove


************************************
Configure paths to required software
************************************

.. config-paths-start-do-not-remove

Paths to required software are specified by modifying the contents of the
:code:`config.yml` file in the root directory.
By default, these paths are empty so they must be configured. Upon using
``scons -h``, a user may see a list of local options for
Upon executing the ``scons -h`` command, one may se a number of local options
including ``--config-software``. Additionally, a user may modify the contents
of :code:`config.yml` directly.

- Configure the paths to required software

  .. code-block::

     $ scons --config-software

The user will be asked if new or additional paths will be appended to the
:code:`config.yml` file. Some of these paths are paths to executable programs
(e.g. Abaqus, Ratel, and Tardigrade-MOOSE), while some are paths to Python
programs and scripts.

The :code:`config.yml` file is read into the SCons configuration file (:code:`SConstruct`).
The YAML file is parsed into a dictionary where each key corresponds to a program and
each entry is a list of program paths.
For exeuctable programs, the :code:`waves.scons_extensions.find_program()` function
is used to search the list of paths with the first executable found being set as
the program path.
For paths to importable Python objects, only the last path in the list is set.

.. config-paths-end-do-not-remove

*******
Testing
*******

.. test-start-do-not-remove

Unlike software projects, the primary model/simulation project tests are the successful completion of some subset of the
simulation targets. If the selected simulations run successfully, then the target passes. Secondary project tests will
use `SCons`_ to execute unit and integration testing for project specific scripts, such as journal files and Python
processing scripts.

- Build the required target(s). Test targets may not be part of the default target list. If so, each target will
  need to be listed explicitly or the "all targets" character, ``.``, should be used to build *all* project targets.

  .. code-block::

     $ scons <target_1_name> <target-2_name>

- Run *all* simulation and test targets. Try to run all targets even if some fail.

  .. code-block::

     scons . --keep-going

.. test-end-do-not-remove

Test Local Module
=================

.. test-local-module-start-do-not-remove

When testing CLI changes locally, the waves module must be run as a script. We must also set the ``PYTHONPATH``
in order to include the current waves module when operating on a configuration that imports waves.

Below is an example of a visualization test of an SConstruct file using the local waves module.

.. code-block::

   $ pwd
   path/to/local/git/clone/waves/
   $ PYTHONPATH=$PWD python -m waves.main visualize . --sconstruct /path/to/local/SConstruct

.. test-local-module-end-do-not-remove

*************
Documentation
*************

.. docs-start-do-not-remove

The documentation build is also automated with SCons as the ``documentation`` target alias.

- Build all documentation targets

  .. code-block::

     $ scons documentation

- Build the HTML documentation

  .. code-block::

     $ scons html

.. docs-end-do-not-remove
