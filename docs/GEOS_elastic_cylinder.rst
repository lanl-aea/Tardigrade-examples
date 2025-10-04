.. _GEOS_elastic:

#####################
GEOS Elastic Cylinder
#####################

A simple simulation is considered to verify that the Micromorphic Filter
may be used to upscale GEOS DNS.

DNS files are copied from the PetaLibrary into a local directory
"peta_data_copy" using the following command. This command usually only needs
to be used once. Files are copied using the secure copy protocal (SCP).
A user will be asked for their identikey, password, and two factor authentication
(2FA).

   .. code-block:: bash

      $ scons --peta-data-copy

The analysis is executed in individual stages by the
:code:`GEOS_elastic_cylinder` SConscript.
First, the existing
DNS results are processed into the required XDMF file format for
the Micromorphic Filter using the following command:

   .. code-block:: bash

      $ scons GEOS_elastic_cylinder_multi_domain

Next the homogenization is performed. Macroscale meshes with
1, 24, 48, and 192 are considered for the default configuration.

   .. code-block:: bash

      $ scons GEOS_elastic_cylinder_multi_domain --filter

***
DNS
***

GEOS is an explicit dynamics code, so simulating quasi-static
loading is challenging.
A displacement of about 0.125 mm is applied to the top face of a cylinder over 100 microseconds.
Boundary conditions are chosen to approximate uniaxial stress.
Figure :numref:`{number} <GEOS_elastic_cylinder_displacements>` shows
the force vs displacement and displacement field of the DNS.
Orange markers in Figure :numref:`{number} <GEOS_elastic_cylinder_displacements>`
denote the frames chosen for upscaling.
The force is roughly linear, though the variations indicate the dynamic nature of the simulation.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: GEOS_elastic_cylinder_displacements
    :class-grid: outline
    :align: center

    .. image:: GEOS_elastic_cylinder_force_displacement.png
       :alt: (a) Force vs displacement with markers denoting frames chosen for upscaling

    .. image:: GEOS_elastic_cylinder_displacement.jpeg
       :alt: (b) Displacement field 

    GEOS elastic cylinder simulation force vs displacement

The Cauchy stress (ZZ component) is plotted in Figure
:numref:`{number} <GEOS_elastic_cylinder_stress>` showing that
stress is approximately uniaxial varying from -6.2 to -6.4 MPa with higher
stresses concentrated on the bottom.

.. figure:: GEOS_elastic_cylinder_stress.jpeg
   :name: GEOS_elastic_cylinder_stress
   :align: center
   :width: 50%

   GEOS elastic cylinder stress

DNS results are converted to the XDMF file format required by the Micromorphic Filter
using the :py:mod:`model_package.DNS_GEOS.vtk_to_xdmf` script.

*********
Upscaling
*********

DNS results are homogenized using the Micromorphic Filter for 1, 24, 48, and 192
filter domains.
Only the results for the 192 element macroscale are presented.
Figure :numref:`{number} <GEOS_elastic_cylinder_filter_results>`
shows the homogenized displacement and stress fields which generally
agree with the results shown in Figures :numref:`{number} <GEOS_elastic_cylinder_displacements>`
and :numref:`{number} <GEOS_elastic_cylinder_stress>`.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: GEOS_elastic_cylinder_filter_results
    :class-grid: outline
    :align: center

    .. image:: GEOS_elastic_cylinder_192_displacement.jpeg
       :alt: (a) Homogenized displacement field

    .. image:: GEOS_elastic_cylinder_192_cauchy33.jpeg
       :alt: (b) Homogenized Cauchy stress ZZ

    Micromorphic Filter results for 192 domains

This study shows that GEOS DNS may be homogenized using the Micromorphic Filter.