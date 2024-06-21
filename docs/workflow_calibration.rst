.. _workflow_calibration:

###########
Calibration
###########

After DNS results have been homogenized through the Micromorphic Filter,
constitutive models are calibrated.
Calibrated models are then applied to Tardigrade-MOOSE simulations.
Figure :numref:`{number} <micromorphic_filter_experimental_data>` illustrates the
calibration methodology in which pairs of stress and deformation quantities output
by the Micromorphic Filter are used as the data to inform constitutive model
parameters.
Calibration is conducted for quantities in the reference configuration.
In reality, the filter outputs 45 stress components (9 for symmetric micro-stress,
9 for Second Piola-Kirchhoff stress, and 27 for higher order stress) and provides
quantities which may be used to compute 45 deformation
tensor components (9 for micro-strain, 9 for Green-Lagrange strain, and 27 for
micro-deformation gradient).

.. figure:: micromorphic_filter_experimental_data.png
   :name: micromorphic_filter_experimental_data
   :align: center
   :width: 50%

   Illustration of calibration method using homogenized stress and deformation
   quantities output by the Micromorphic Filter

As discussed in :ref:`macroscale_definition`, either single or multiple
filtering domains may be used for homogenization.
For upscaling studies in this repository,
each filtering domain is a trilinear hexahedral finite element.
All homogenized quantities (see Eq. :math:numref:`{number} <homogenized_quantities>`
and Eq. :math:numref:`{number} <deformation_measures>`) output by the Micromorphic Filter
are provided at the 8 integration points of the filter domain / element.

While calibration could be conducted using the output for each integration point,
the current setup of macroscale Tardigrade-MOOSE simulations only allows for
unique material properties to be assigned to "blocks" which contain either one or
more elements.
As such, quantities are first averaged to the element centroid prior to calibration.
In most cases, unique material inputs are assigned to each finite element in
macroscale simulations.
Figure :numref:`{number} <single_element_calibration>` illustrates this approach for
a single filter element.

.. figure:: single_element_calibration.png
   :name: single_element_calibration
   :align: center
   :width: 50%

   Micromorphic material calibration for a single macroscale finite element

The `WAVES`_ computational workflow tool provides an efficient way to parallelize
calibration efforts. While the Micromorphic Filter conducts homogenization
over the entire macroscale to ensure that kinematic quantities are
continuous across the finite elements nodes, subsequent calibration
tasks may be conducted independently.
This is convenient because there may be many macroscale elements to
be calibrated and is a rather computationally expensive process.
For example, in a multiple filtering domain study with macroscale meshes
containing 1, 10, 100, and 1000 elements, there are 1111 calibrations to perform.
This process may be parallelized using the :code:`--jobs=n` argument as follows:

   .. code-block:: console

      $ scons study_name --calibrate --jobs=n

Calibration is performed using the :py:mod:`model_package.Calibrate.calibrate_element`
script. This script accepts several command line arguments including:

* the XDMF file containing homogenization results output by the Micromorphic Filter,
* the name of an output YAML file containing calibrated parameters to be used for
  macroscale Tardigrade-MOOSE simulations,
* approximate, homogenized DNS material properties used for initial parameter estimation,
* the element number of the macroscale mesh to be calibrated,
* optional time increment to perform calibration,
* optional file name to plot comparison between DNS and calibration results for Cauchy or
  symmetric micro-stress, and
* the calibration "case" to be discussed below.

See the API documentation and source code for :py:mod:`model_package.Calibrate.calibrate_element`
for futher details.

Currently, only linear micromorphic elasticity models are considered.
Different calibration "cases" may be specified that correspond to different versions
of linear micromorphic elasticity:

* case 1: classical elasticity simplification (see Eq. :math:numref:`{number} <constitutive_case_1>`),
* case 2: 7 parameter elasticity simplification with :math:`\tau_7^*` fixed at :math:`10^{-3}
  MPa \cdot mm^{\text{2}}`,
  (see Eq. :math:numref:`{number} <constitutive_case_2&3>`),
* case 3: 8 parameter elasticity simplification
  (see Eq. :math:numref:`{number} <constitutive_case_2&3>`), and
* case 4: full 18 parameter linear micromorphic elasticity model
  (see Eq. :math:numref:`{number} <constitutive_case_4>`).

For cases 2 and 3, the objective function is constructed using the error between predicted
and homogenized values of the Second Piola-Kirchhoff and symmetric micro-stresses.
For case 4, the higher-order stress error is included.
For case 1, only the error for the :math:`33` components of Second Piola-Kirchhoff and
symmetric micro-stress is calculated. It is assumed that all DNS used for calibration
with case 1 are loaded in the direction that creates :math:`33` stresses.
To help constrain the parameter estimation for this case, a target Poisson ratio is
first calculated from the homogenized Green-Lagrange strains.
During optimization, a trial pair of Lam\'e parameters is converted to a Poisson
ratio using :math:`\nu^* = \,^{\lambda^*}\!/\!_{2 \left(\lambda^* + \mu^* \right)}`. If this trial
Poisson ratio is not within :math:`\pm 1\%` of the target, the trial parameters are rejected.

To better understand the details associated with calibration activities,
consider inspecting the :code:`calibrate_element.scons` SConscript located in
the :code:`model_package/workflows` directory.