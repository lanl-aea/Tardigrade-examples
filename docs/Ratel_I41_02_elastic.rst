.. _Ratel_I41_02_elastic:

#####################################################
Ratel I41_02 elastic cylinder - Quasi-static Implicit
#####################################################

This study investigates micromorphic upscaling of a heterogeneous
material composed of grains and a binder.
The DNS is an implicit finite element
(FE) model assuming quasi-static conditions conducted in Ratel.
The results of this DNS are prepared for the Micromorphic Filter
for a variety of filtering domains.
Filter output is then used to calibrate an elastic micromorphic material
model which is implemented for simulation in Tardigrade-MOOSE.
Unique material calibrations are produced for each filtering domain
which are prescribed to their corresponding elements in the
macroscale simulations.

A variety of simulation variables are provided in the :code:`I41_02`
dictionary in the :code:`model_package/DNS_Ratel/simulation_variables_nominal.py` file.
This dictionary is loaded into the workflow as the :code:`params` dictionary.

.. literalinclude:: DNS_Ratel_simulation_variables_nominal.py
   :lines: 95-116
   :linenos:

.. warning::

   The DNS files needed for this upscaling study can only be access by users
   with a CU Boulder identikey and access to the PetaLibrary!

DNS files are copied from the PetaLibrary into a local directory
"peta_data_copy" using the following command. This command usually only needs
to be used once. Files are copied using the secure copy protocal (SCP).
A user will be asked for their identikey, password, and two factor authentication
(2FA).

   .. code-block:: bash

      $ scons --peta-data-copy

The analysis is executed in individual stages by the
:code:`Ratel_I41_02_elastic_multi_domain` SConscript.
First, the existing
DNS results are processed into the required XDMF file format for
the Micromorphic Filter using the following command:

   .. code-block:: bash

      $ scons Ratel_I41_02_elastic_multi_domain

Next the homogenization is performed. Macroscale meshes with
1, 24, 48, 192, and 960 elements are considered for the default configuration.
The micromorphic filter is parallelized to run on 8 cpus. This value may be
modified by changing the value of "filter_parallel" in the :code:`I41_02`
parameter dictionary.

   .. code-block:: bash

      $ scons Ratel_I41_02_elastic_multi_domain --filter

Calibration is then performed for each macroscale element (i.e. filtering domain).
This process can be rather expensive depending on what model is being calibrated,
however, WAVES provides a simple way to parallelize the independent calibration
processes using the :code:`--jobs=N` command line option. Calibration may be
performed using a single process or multiple (10 in this example) with
either of the following two options, respectively:

   .. code-block:: bash

      $ scons Ratel_I41_02_elastic_multi_domain --calibrate
      $ scons Ratel_I41_02_elastic_multi_domain --calibrate --jobs=10

Once calibration is completed, Tardigrade-MOOSE simulations may be performed.
A Tardigrade-MOOSE simulation is performed for each case of filtering domains.
Macroscale simulations may be parallelized using the :code:`--solve-cpu=N`
command line option. Macroscale simulations may be performed using a single
process or multiple (12 in this example) using either of the following
two options, respectively:

   .. code-block:: bash

        $ scons Ratel_I41_02_elastic_multi_domain --macro
        $ scons Ratel_I41_02_elastic_multi_domain --macro --solve-cpus=12

Finally, the results across filtering domains may be summarized into several
plots and csv files using the "summary" command:

   .. code-block:: bash

        $ scons Ratel_I41_02_elastic_multi_domain --summary

***************************
DNS Description and Results
***************************

**Geometry and Mesh**

This DNS is created from the I41.02 CT scan by the CU Boulder `PSAAP`_ center.
The geometry is 6 mm in diameter with a height of 5.485 mm.
A coarse mesh is generated with 591,030 nodes.

**Material**

For this study, nominal elastic properties for the binder are chosen
with an elastic modulus, :math:`E^*_{binder}`, of 230.0 MPa and a Poisson's ratio,
:math:`\nu^*_{binder}`, of 0.45. The elastic properties of the grains
are chosen with an elastic modulus, :math:`E^*_{grain}` of 22000.0 MPa
and a Poisson's ratio, :math:`\nu^*_{grain}` of 0.25.
The binder density is 1910 :math:`\frac{kg}{m^3}` and the grain density
is 2001 :math:`\frac{kg}{m^3}`.
Figure :numref:`{number} <Ratel_I41_density>` shows mass density contour plots
of the Ratel I41.02 geometry where it is clear that a sparse collection of
grains are interspersed in the binder matrix.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_density
    :class-grid: outline
    :align: center
    :width: 70%

    .. image:: Ratel_I41_DNS_density_transparent.jpeg
       :alt: (a) Low densities have lower opacity

    .. image:: Ratel_I41_DNS_density_cut.jpeg
       :alt: (b) cut view

    Density contours of Ratel I41.02 DNS

..
   TODO: describe experimental modulus results and how a composite E* and nu* are calculated for the upscaling

**Boundary conditions and loading**

The bottom face of the cylinder is fixed in the z-direction.
Both the top and bottom faces are restrained from lateral expansion,
which is referred to as "clamped" boundary conditions.
The top face of the cylinder is compressed 0.2729 mm, corresponding
to a nominal compressive strain of nearly 5 percent.

**Results**

Figure :numref:`{number} <Ratel_I41_stress>` shows the resulting Cauchy stress
contours (33 component). It is evident that there are high stress concentrations
near the clamped boundaries on the top and bottom of the cylinder in
figure :numref:`{number} <Ratel_I41_stress>` (a). Similarly, stresses
concentrate near adjacent grains shown in figure
:numref:`{number} <Ratel_I41_stress>` (b). Note that these stress concentrations
are certainly in similar locations as the grains visible in
figure :numref:`{number} <Ratel_I41_density>` (b).

.. subfigure:: AA|BC
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_stress
    :class-grid: outline
    :align: center
    :width: 70%

    .. image:: Ratel_I41_stress_colormap.jpeg

    .. image:: Ratel_I41_DNS.jpeg
       :alt: (a) full geometry

    .. image:: Ratel_I41_DNS_with_cut.jpeg
       :alt: (b) cut view

    Cauchy stress contours (33 component) of Ratel I41.02 DNS

The following figure shows the resultant force versus displacement for the top of the cylinder,
with a total force of 489.9 N.

.. figure:: Ratel_I41_elastic_fd.png
   :name: Ratel_I41_elastic_fd
   :align: center
   :width: 50%

   Ratel I41.02 DNS force versus displacement results

******************
Filter Preparation
******************

DNS results are converted to the required XDMF file format in a very similar
way to what is describd in :ref:`Ratel_elastic_cylinder`.
Here, filtering domains with 1, 24, 48, 192, and 960
macroscale elements are used, however, only the results
of the 48, 192, and 960 element cases will be discussed.

**************
Filter Results
**************

Figure :numref:`{number} <Ratel_I41_filter_cauchy>` shows the homogenized
Cauchy stresses (33 component) output by the Micromorphic Filter for the
Ratel I41.02 DNS. It is clear that, as filtering domains decrease in
size, more of the DNS stress response is recovered. The finest filtering
domain case of 960 elements shows some stress concentration near the
edges of the cylinder top and bottom surfaces, though it is more
smeared out than the DNS.

.. subfigure:: AAA|BCD
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_filter_cauchy
    :class-grid: outline
    :align: center
    :width: 70%

    .. image:: Ratel_I41_stress_colormap.jpeg

    .. image:: Ratel_I41_filter_results_48.jpeg
       :alt: (a) 48 domains

    .. image:: Ratel_I41_filter_results_192.jpeg
       :alt: (b) 192 domains

    .. image:: Ratel_I41_filter_results_960.jpeg
       :alt: (c) 960 domains

    Homogenized Cauchy stress contours (33 component) output by Micromorphic Filter for the Ratel I41.02 DNS

The distribution in Cauchy stresses are further inspected using violin plots in
figure :numref:`{number} <Ratel_I41_violin_plot>`.
The minimum, mean, and maximum values are shown for each filtering domain case,
and the DNS,
as horizontal bars while a kernel density estimate (KDE) is used to illustrated the
distribution of :math:`\sigma_{33}`.
Here it is clear that variation in stress increases as the number of filtering
domains increase. The mean stress is roughly consistent across all plots.
Note that the variation in DNS stresses is extremely high with a max
compressive stress larger than 500 MPa and a max tensile stress
greater than 100 MPa.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_violin_plot
    :class-grid: outline
    :align: center
    :width: 70%

    .. image:: Ratel_I41_cauchy33_violin_plot.png
       :alt: (a) Homogenized

    .. image:: Ratel_I41_cauchy33_violin_plot_DNS.png
       :alt: (b) DNS
       :scale: 70%

    Violin plots of Cauchy stress (33 component)

*******************************************
Micromorphic Constitutive Model Calibration
*******************************************

The 8 parameter micromorphic linear elasticity model (Eq.
:math:numref:`{number} <constitutive_case_2&3>`) is calibrated for every filtering
domain. The parameter bounds for :math:`\lambda` is from 0 to 100,000 MPa, while
the bounds for :math:`\mu`, :math:`\eta`, :math:`\tau`, :math:`\kappa`,
:math:`\nu`, and :math:`\sigma` are -100,000 to 100,000 MPa, and finally
the bounds for :math:`\tau_7` are -100,000 to 100,000 MPa-mm^2.
The calibration routine may be inspected in
:py:mod:`model_package.Calibrate.calibrate_element`.
The entirety of the calibration process is summarized in the joint
probability distribution plot in figure
:numref:`{number} <Ratel_I41_joint_probability_distributions_case_3>`.

.. figure:: Ratel_I41_joint_probability_distributions_case_3.png
   :name: Ratel_I41_joint_probability_distributions_case_3
   :align: center
   :width: 90%

   Joint Probability Distribution of Ratel I41.02 calibration results

The histograms (which are stacked) shown on the diagonal of the joint plot
indicate that most of the calibrations are similar, but the off-diagonals
sub-plots show that parameter calibrations may deviate significantly from
the average. This deviation increases dramatically as the number of filtering
domains increase (e.g., the 960 domain case shown in blue).
Additionally, it is clear that most parameter calibrations are
strongly correlated (such as the plot for :math:`\sigma` compared with
:math:`\mu`). It is difficult to determine if these correlations are
"real" or just a byproduct of the approach.

The homogenized stresses contours shown
in figure :numref:`{number} <Ratel_I41_filter_cauchy>` indicate that
the stress concentrations, present in the DNS begin, are more clearly
resolved as the size of the filtering domains decrease.
It is possible that these regions
of high stress may be affecting calibration results.
Therefore, the joint probability distribution plot is generated a second
time with all of the filtering domains (macroscale elements) on the
top and bottom surface excluded, shown in figure
:numref:`{number} <Ratel_I41_joint_probability_distributions_case_3_no_BCs>`.
Note that all elements of the 24 domain case touch the boundary, so only
the 48, 192, and 960 domain cases are relevant.

.. figure:: Ratel_I41_joint_probability_distributions_case_3_no_BCs.png
   :name: Ratel_I41_joint_probability_distributions_case_3_no_BCs
   :align: center
   :width: 90%

   Joint Probability Distribution of Ratel I41.02 calibration results with boundary elements excluded.

It is interesting to see that there is still a large spread in parameter calibrations.

As it pertains to the macroscale simulations in Tardigrade-MOOSE, applying
the unique calibrations from
:numref:`{number} <Ratel_I41_joint_probability_distributions_case_3>`
was found to be problematic and Tardigrade-MOOSE would crash for the 192 and
960 domain macroscale simulations. As such, it was decided to ignore calibration
for domains located on the top and bottom boundaries. However, those elements
still need a material assigned. It was discovered that taking the peak value
for each parameter from a kernel density estimate of the calibration results
(neglecting boundary elements) produced satisfactory inputs for heterogeneous
Tardigrade-MOOSE simulations. These kernel density estimates are shown
in figure :numref:`{number} <Ratel_I41_kdes>`.

.. subfigure:: ABCD|EFGH
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_kdes
    :class-grid: outline
    :align: center
    :width: 90%

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_lambda.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_mu.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_eta.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_tau.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_kappa.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_nu.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_sigma.png

    .. image:: Ratel_I41_full_kde_case_3_no_BCS_tau_7.png

    Kernel density estimates for summarizing calibration results for all parameters for 48, 192, and 960 filtering domains (boundary elements excluded) with peak KDE value highlighted

On the other hand, it is desirable to keep the calibrations for different filtering domain
cases separate. As such, the peak KDE value is extracted for each parameter for each
filter domain cases and summarized in table
:numref:`{number} <Ratel_I41_best_KDE_table>`.

.. list-table:: Summary of peak KDE value for each filter domain case and parameter
   :name: Ratel_I41_best_KDE_table
   :widths: 20 20 20 20
   :header-rows: 1

   * - Parameter
     - 48 domains
     - 192 domains
     - 960 domains
   * - :math:`\lambda`
     - 729.8
     - 684.3
     - 645.4
   * - :math:`\mu`
     - 122.4
     - 127.2
     - 119.5
   * - :math:`\eta`
     - -13.91
     - -6.464
     - 14.54
   * - :math:`\tau`
     - -22.82
     - -25.37
     - 4.823
   * - :math:`\kappa`
     - 28.01
     - 12.69
     - 10.66
   * - :math:`\nu`
     - -37.07
     - -24.50
     - -15.86
   * - :math:`\sigma`
     - -2.714
     - -5.979
     - -4.888
   * - :math:`\tau_7`
     - 730.5
     - 309.6
     - 32.81

*********************
Macroscale Simulation
*********************

Macro-scale simulations are now conducted in Tardigrade-MOOSE using the unique calibrations for
each filtering domain applied to identical macro-scale meshes.
Elements touching the boundary have material inputs assigned based on the peak KDE
values shown in Table :numref:`{number} <Ratel_I41_best_KDE_table>` for the
different filter domain cases.
Clamped boundary conditions are applied with a compressive displacement of 0.3 mm.
This workflow stage is run using the ``--macro-ignore-BCs`` option which is configured
to apply peak KDE calibrations for boundary elements.

   .. code-block:: bash

        $ scons Ratel_I41_02_elastic_multi_domain --macro-ignore-BCs

The Second Piola-Kirchhoff stress contour (33 compoent) results
are shown in figure :numref:`{number} <Ratel_I41_macroscale_PK2>`.
A force-displacement plot is shown in figure
:numref:`{number} <Ratel_I41_all_force_displacements_case_3>` compared with
the DNS results showing good agreement.

.. subfigure:: AAA|BCD
    :gap: 8px
    :subcaptions: below
    :name: Ratel_I41_macroscale_PK2
    :class-grid: outline
    :align: center
    :width: 70%

    .. image:: Ratel_I41_pk2_colorbar.jpg

    .. image:: Ratel_I41_macroscale_elastic_48.jpeg
       :alt: (a) 48 domains

    .. image:: Ratel_I41_macroscale_elastic_192.jpeg
       :alt: (b) 192 domains

    .. image:: Ratel_I41_macroscale_elastic_960.jpeg
       :alt: (c) 960 domains

    Second Piola-Kirchhoff stress contours (33 component) for Tardigrade-MOOSE simulations

.. figure:: Ratel_I41_all_force_displacements_case_3.png
   :name: Ratel_I41_all_force_displacements_case_3
   :align: center
   :width: 50%

   Force vs. displacement results for Tardigrade-MOOSE simulations and Ratel I41.02 DNS