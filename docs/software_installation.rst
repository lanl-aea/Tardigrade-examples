.. _software_installation:

#####################
Software Installation
#####################

*************************************
Determining what software is required
*************************************

Some users do not need to install and link every piece of software
depending on their use case.

Verification workflows require installations of :ref:`abaqus_fem_installation`
and/or :ref:`ratel_fem_installation`.
Workflows for upscaling heterogeneous DNS use simulations stored on the
PetaLibrary and can be accessed by following the instructions provided
here: :ref:`peta_library`.

All upscaling studies require the installation of the
:ref:`micromorphic_filter_installation`.
Several workflows only run Tardigrade-MOOSE simulations.

All upscaling studies perform calibration of material models which
require the installation of the
:ref:`micromorphic_calibration_tool_installation` and
:ref:`micromorphic_elastic_constraints_installation` tool.

All macroscale and upscaling studies require the installation of
:ref:`tardigrade_moose_installation`.

********************
Activate environment
********************

It is assumed that the tardigrade-examples-env environment has been installed
per the instructions provided in :ref:`build`.

   .. code-block:: console

      $ conda activate -n tardigrade-examples-env

.. _abaqus_fem_installation:

**********
Abaqus FEM
**********

Abaqus is used for several upscaling verification workflows, but is otherwise
optional. Users may refer to the `Abaqus`_ documentation for installation
instructions.

Add Abaqus to software configuration path
=========================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/abaqus` to the :code:`config_software.yml` entry for "Abaqus".

For most Windows installations,
the path to the :code:`abaqus.bat` script may be specified. The default
path is already added to the :code:`SConstruct` file assuming the program
is installed in the default location.

.. _ratel_fem_installation:

*********
Ratel FEM
*********

Ratel may be installed according to directions provided in the main repository,
https://ratel.micromorph.org/doc/intro/,
however, further instructions are provided here for clarity.
For the DNS required by workflows of this repository, the only
Ratel software prerequisites are libCEED and PETSc.
Once libCEED and PETSc are properly installed, Ratel may be built.

First, create a root directory for Ratel software in a convenient location.

   .. code-block:: console

      $ mkdir /path/to/directory/Ratel
      $ export RATEL_DIR=/path/to/directory/Ratel
      $ cd $RATEL_DIR

Build libCEED
=============

Clone and build libCEED.

   .. code-block:: console

      $ cd $RATEL_DIR
      $ git clone https://github.com/CEED/libCEED
      $ cd libCEED
      $ make

Export the build directory (``$CEED_DIR``) as an environment variable.

   .. code-block:: console

      $ export CEED_DIR=$RATEL_DIR/libCEED

Build PETSc
===========

Clone, configure, and build PETSc then export the build location as an
environment variable. The following instructions provided by PETSc may
be a helpful reference if problems arise: https://petsc.org/release/install/install_tutorial/#qqtw-quickest-quick-start-in-the-west.
Several configure options are specified here to allow Ratel simulations
to run with Exodus meshes generated using Cubit.

   .. code-block:: console

      $ cd $RATEL_DIR
      $ git clone https://gitlab.com/petsc/petsc
      $ cd petsc
      $ ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack --download-exodusii --download-hdf5 --download-netcdf --download-pnetcdf --download-zlib

After configuring PETSc, a specific ``make`` command will be provided that
lists the PETSc build architercture (``$PETSC_ARCH``) and the PETSc build
directory (``$PETSC_DIR``) which will be needed for Ratel. The following ``make``
command was provided when building on ``sstbigbird.lanl.gov``.

   .. code-block:: console

      $ make PETSC_DIR=$RATEL_DIR/petsc PETSC_ARCH=arch-linux-debug all

Export the build directory (``$PETSC_DIR``) and build architecture (``$PETSC_ARCH``)
as environment variables.

   .. code-block:: console

      export PETSC_DIR=$RATEL_DIR/petsc
      export PETSC_ARCH=arch-linux-debug

.. note::

   Make sure to use the :code:`PETSC_ARCH` specified by PETSc after the configuration step

Build Ratel
===========

Clone and build Ratel. This should work if the ``$CEED_DIR``, ``$PETSC_DIR``, and
``$PETSC_ARCH`` environment variables have been set.

   .. code-block:: console

      $ cd $RATEL_DIR
      $ git clone https://gitlab.com/micromorph/ratel
      $ cd ratel
      $ make

Test
====

The Ratel documentation includes instructions for how to test the installation
which a user is welcome to follow. Another simple test may be run using the
following commands:

   .. code-block:: console

      $ cd $RATEL_DIR
      $ ./bin/ratel-quasistatic -options_file examples/ex02-quasistatic-elasticity-linear-platen.yml

Many other examples can be found in the :code:`$RATEL_DIR/examples` directory.

Add Ratel to software configuration path
========================================

Currently, all Ratel DNS used in this repository only require the `ratel-quasistatic` program.
This executable should be located in ``$RATEL_DIR/ratel/bin/ratel-quasistatic``.
Either using :code:`scons --config-software` or manually, add
:code:`/path/to/ratel/bin/ratel-quasistatic` to the :code:`config_software.yml` entry for "Ratel".

********
GEOS MPM
********

..
   TODO: Describe how to build and link GEOS MPM

Coming soon!

*****
Cubit
*****

Cubit is used for a number of meshing operations.
Users may refer to the `Cubit`_ documentation for installation instructions.

For users without access to Cubit,
several example meshes are contained in :code:`model_package/meshes/`, however, functionality
of workflows will be limited.

Add Cubit to software configuration path
========================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/cubit` to the :code:`config_software.yml` entry for "Cubit".

.. _micromorphic_filter_installation:

*******************
Micromorphic Filter
*******************

All workflows use the Micromorphic Filter for homogenization. This software
is written entirely in Python and does not need to be compiled or built in any
capacity. Workflows using the Micromorphic Filter are already configured to
instantiate the Filter class and call relevant functions. Simply clone the
repository to a desired location.

   .. code-block:: console

      $ git clone git@github.com:UCBoulder/tardigrade_filter.git


In order to clone this repository, a user may need to configure their
GitHub account to be associated with University of Colorado Boulder's
single sign-on (SSO). For instructions, see the section titled
"Access GitHub" from the Office of Information Technology at the
following link:
https://oit.colorado.edu/services/business-services/github-enterprise

The Conda Environment for this repo includes all of the same packages
included in the Micrormophic Filter repository to guarantee that this
software functions appropriately.

Test
====

The Micromorphic Filter comes with built in tests using PyTest. To run these
tests, simply run the following commands:

   .. code-block:: console

      $ cd /path/to/tardigrade_filter
      $ pytest

Add Micromorphic Filter to software configuration path
======================================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/tardigrade_filter/src/python` to the
:code:`config_software.yml` entry for "filter".

The path to the Micromorphic Filter's :code:`src/python` directory needs to be inserted
into the Python path whenever it is to be used. This is handled automatically by
the SCons workflow.

.. _tardigrade_moose_installation:

****************
Tardigrade-MOOSE
****************

Tardigrade-MOOSE is built using CMake and requires a number of compilers and 
Python libraries which are included in the :code:`environment.txt` file included
in this repository.

.. note::

   Note that `MOOSE`_ and associated Python package update frequently,
   so the conda environment for this repository should be rebuilt each time
   Tardigrade-MOOSE is to be compiled. See the following link for more
   information: https://mooseframework.inl.gov/getting_started/new_users.html#update.

Clone Tardigrade
================

   .. code-block:: console

      $ git clone https://github.com/UCBoulder/tardigrade.git
      $ cd tardigrade

CMake
=====

   .. code-block:: console

      $ mkdir build
      $ cd build
      $ cmake .. -DTARDIGRADE_BUILD_PYTHON_BINDINGS=OFF
      $ make -j 4

.. _LD_PATH_NOTE:

Set LD_LIBRARY_PATH
===================

There is an LD_LIBRARY_PATH that needs to be specified.
A user may either:
(1) export this path as an environment variable or
(2) include this path on the command line each time a Tardigrade package is run.

For option 1, an environment variable may be set with the following command.
It is NOT recommended to include this environment variable in a ~/.bashrc as
there may be unintended consequences.

   .. code-block:: console

      $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp

For details using option 2, see the following subsection for "Test" or :ref:`macroscale_command_line`.
**Workflows that run Tardigrade-MOOSE are configured to automatically use option 2 in which the
LD_LIBRARY_PATH is prepended to the command that launches a simulation.** However, note that
other operations may still require manual intervention (such as those described in the sections
just mentioned).

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp` to the
:code:`config_software.yml` entry for "LD_PATH". This configuration will ensure that
Tardigrade-MOOSE simulations run through SCons workflows will access the appropriate shared libraries.

If one encounters :code:`error while loading shared libraries: libmicromat.so: cannot open shared object file`,
then the LD_LIBRARY_PATH is not configured correctly.

Test
====

The tests may be run using the :code:`ctest -v` command from within the Tardigrade build directory.
As discussed in :ref:`LD_PATH_NOTE`, the tests may be run with the LD_LIBRARY_PATH already
set as an environment variable with:

   .. code-block:: console

      $ cd /path/to/tardigrade/build
      $ ctest -v

or by specifying the LD_LIBRARY_PATH on the command line:

   .. code-block:: console

      $ cd /path/to/tardigrade/build
      $ LD_LIBRARY_PATH=/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-build/src/cpp ctest -v

Most or all of the tests should pass. If they do not all pass, the tests may be run again with the
:code:`--rerun-failed` and :code:`--output-on-failure` to see what tests failed. If a test fails
with the "EXODIFF" reason, then it is likely that the most recent of Tardigade produces output
that does not exactly match the "gold" results file. Otherwise, if tests fail because a specific
library is not found (e.g., :code:`libmicromat.so`) then Tardigrade is configured improperly and/or
the LD_LIBRARY_PATH has not been specified correctly.

Add Tardigrade-MOOSE to software configuration path
===================================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/tardigrade/build/tardigrade-opt` to the
:code:`config_software.yml` entry for "Tardigrade".

.. _micromorphic_calibration_tool_installation:

*****************************
Micromorphic Calibration Tool
*****************************

The micromorphic calibration tool is a shared Python library that can be
built after building :code:`tardigrade_micromorphic_element`.
This tool may be installed either (1) using the repository cloned automatically
by the :ref:`tardigrade_moose_installation` CMake build, or (2) after cloning
the Micromorphic Element repository.

Option 1. Build from repository cloned by Tardigrade-MOOSE
==========================================================

If the Tardigrade-MOOSE build went smoothly, then the
calibration tool will be contained in the
:code:`/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-src/src/python`
directory.

.. note::

   It is likely that the :code:`setup.py` file will need to be modified!

Set the :code:`library_dirs` in :code:`setup.py` to the following path:

   .. code-block:: python

      library_dirs = [os.path.abspath('../../../tardigrade_micromorphic_element-build/src/cpp')]

The LD_LIBRARY_PATH must be set according to the instuctions provided in :ref:`LD_PATH_NOTE`.

The shared library may be built as follows:

   .. code-block:: console

      $ cd /path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-src/src/python
      $ python setup.py build_ext --inplace

Option 2. Clone and build from Micromorphic Element Repository
==============================================================

Building Tardigrade-MOOSE can be challenging and unnecessary if a user is only performing
workflows for homogenization and material model calibration without needing to run
macroscale simulations.

Clone Micromorphic Element
--------------------------

   .. code-block:: console

      $ git clone https://github.com/UCBoulder/tardigrade_micromorphic_element
      $ tardigrade_micromorphic_element

CMake
-----

.. note::

   It is likely that the `CMakeLists.txt` file will need to be modified!
   Make sure that :code:`"tardigrade_stress_tools" "tardigrade_solver_tools" "tardigrade_hydra" "tardigrade_constitutive_tools"`
   are added to the upstream packages list.

Perform the CMake build.

   .. code-block:: console

      $ mkdir build
      $ cd build
      $ cmake .. -DTARDIGRADE_MICROMORPHIC_BUILD_PYTHON_BINDINGS=OFF
      $ make -j 4

Build Calibration Tool
----------------------

The calibration tool will be contained in the
:code:`/path/to/tardigrade_micromorphic_element/src/python` directory.

.. note::

   It is likely that the :code:`setup.py` file will need to be modified!

Set the :code:`library_dirs` in :code:`setup.py` to the following path:

   .. code-block:: python

      library_dirs = [os.path.abspath('../../build/src/cpp')]

The LD_LIBRARY_PATH must be set according to the instuctions provided in :ref:`LD_PATH_NOTE`.

The shared library may be built as follows:

   .. code-block:: console

      $ cd /path/to/tardigrade_micromorphic_element/src/python
      $ python setup.py build_ext --inplace

Test
====

To test that the shared library is working correctly, one may start
an interactive Python session in the directory containing `setup.py`
and use :code:`import micromorphic`. Similarly, an interactive session may be run
from any directory, but the location of the micromorphic shared library must be
appended to the Python path as follows:

   .. code-block:: python

      import sys
      sys.path.append('/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-src/src/python')
      import micromorphic

Further discussion is provided in :ref:`software_usage` to show how the WAVES workflow
automatically sets these Python paths. 

Add Micromorphic Calibration Tool to software configuration path
================================================================

Either using :code:`scons --config-software` or manually, add

* :code:`/path/to/tardigrade/build/_deps/tardigrade_micromorphic_element-src/src/python`
  if using the :ref:`tardigrade_moose_installation` CMake build, or
* :code:`/path/to/tardigrade_micromorphic_element/src/python`
  if using the :ref:`micromorphic_calibration_tool_installation` (option 2) CMake build

to the :code:`config_software.yml` entry for "micromorphic".

The path to the :code:`micromorphic` shared library needs to be inserted
into the Python path whenever it is to be used. This is handled automatically by
the SCons workflow.

.. _micromorphic_elastic_constraints_installation:

***************************************
Micromorphic Linear Elastic Constraints
***************************************

Constraints of the micromorphic linear elasticity model of Eringen and Suhubi
:cite:`eringen_nonlinear_1964` must be enforced. See discussion of these
constraints in :ref:`linear_elastic_constraints`.

The calibration stage of upscaling workflows must evaluate these constraints
when determining linear elastic parameters. 
The :code:`linear_elastic_parameter_constraint_equations.py` script is provided in
the :code:`tardigrade_micromorphic_linear_elasticity` repository to
evluate these 13 constraints. This repository is automatically pulled
during the CMake builds for :ref:`tardigrade_moose_installation` or
:ref:`micromorphic_calibration_tool_installation` (option 2).

Add Micromorphic Linear Elastic Constraints to software configuration path
==========================================================================

Either using :code:`scons --config-software` or manually, add

* :code:`/path/to/tardigrade/build/_deps/tardigrade_micromorphic_linear_elasticity-src/src/python`
  if using the :ref:`tardigrade_moose_installation` CMake build, or
* :code:`/path/to/tardigrade_micromorphic_element/build/_deps/tardigrade_micromorphic_linear_elasticity-src/src/python`
  if using the :ref:`micromorphic_calibration_tool_installation` (option 2) CMake build

to the :code:`config_software.yml` entry for "constraints".

The path to the :code:`linear_elastic_parameter_constraint_equations.py` script needs to be inserted
into the Python path whenever it is to be used. This is handled automatically by
the SCons workflow.

.. _mpi:

***
MPI
***

Parallel jobs for Ratel and Tardigrade-MOOSE may be run using MPI (message passing interface).
The location of the :code:`mpiexec` utility will depend on the system being used,
however, it may have been installed when creating the conda environment for
this project (i.e. :code:`/path/to/tardigrade-examples-env/bin/mpiexec`).
One may be able to locate this utility by executing :code:`which mpiexec`
on the command line.

The mpiexec command should only be necessary for parallelizing simulations run
on systems without a job scheduler such as SLURM. For HPCs with SLURM, see the
discussion in :ref:`serial_vs_parallel`.

Add MPI to software configuration path
======================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/mpiexec`
to the :code:`config_software.yml` entry for "mpi".

*****
Neper
*****

Neper should be installed to the Python environment discussed in :ref:`build`.
One may test if neper has been installed using:

   .. code-block:: console

      $ which neper

This command will provide the path to the Neper program.
If Neper is not found, one may try installing it into the Python environment
using :code:`conda install neper`.

Add Neper to software configuration path
=========================================

Either using :code:`scons --config-software` or manually, add
:code:`/path/to/neper`
to the :code:`config_software.yml` entry for "Neper".