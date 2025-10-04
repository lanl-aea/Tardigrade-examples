.. _GEOS_I43_01:

#####################
GEOS I43_01 Upscaling
#####################

This study investigates micromorphic upscaling of a heterogeneous
material composed of grains and a binder, specifically the I43.01 specimen.
The DNS is an explicit material point method (MPM) model of a dynamic experimental
conducted in GEOS.
The results of this DNS are homogenized using the Micromorphic Filter
for a variety of filtering domains.

DNS files are copied from the PetaLibrary into a local directory
"peta_data_copy" using the following command. This command usually only needs
to be used once. Files are copied using the secure copy protocal (SCP).
A user will be asked for their identikey, password, and two factor authentication
(2FA).

   .. code-block:: bash

      $ scons --peta-data-copy

The analysis is executed in individual stages by the
:code:`GEOS_I43_01_sim38` SConscript.
First, the existing
DNS results are processed into the required XDMF file format for
the Micromorphic Filter using the following command:

   .. code-block:: bash

      $ scons GEOS_I43_01_sim38

Next the homogenization is performed. Macroscale meshes with
1, 24, 48, and 192 are considered for the default configuration.

   .. code-block:: bash

      $ scons GEOS_I43_01_sim38 --filter

***
DNS
***

A displacement of about 0.2 mm is applied over 121 microseconds.
The output files of this DNS are very large and difficult to visualize on
a personal machine.
The following figures show the extracted point fields input to the
Micromorphic Filter containing the raw DNS data.
Figure :numref:`{number} <GEOS_I43_01_cylinder_displacements>` shows
the force vs displacement and displacement field of the DNS.
Orange markers in Figure :numref:`{number} <GEOS_I43_01_cylinder_displacements>`
denote the frames chosen for upscaling.
The force-displacement indicates that the specimen has started to fail by the end
of the simulation.
The displacement field also shows that the top of the cylinder is not perfectly
which may be a byproduct of the CT segmentation and DNS geometry generation.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: GEOS_I43_01_cylinder_displacements
    :class-grid: outline
    :align: center

    .. image:: GEOS_I43_01_sim38_force_displacement.png
       :alt: (a) Force vs displacement with markers denoting frames chosen for upscaling

    .. image:: GEOS_I43_01_sim38_displacement.jpeg
       :alt: (b) Displacement field 

    GEOS I43.01 simulation force vs displacement

The Cauchy stress (ZZ component) and scalar damage field is plotted in Figure
:numref:`{number} <GEOS_I43_01_cylinder_stress_damage>`.
It is clear that there is a high concentration of stress on the top of the cylinder.
The damage field shows regions, resembling grains, that have reached full
state of damage.
These same regions of high damage correspond to regions where Cauchy stress is zero.

.. subfigure:: AB
    :gap: 8px
    :subcaptions: below
    :name: GEOS_I43_01_cylinder_stress_damage
    :class-grid: outline
    :align: center

    .. image:: GEOS_I43_01_sim38_stress.jpeg
       :alt: (a) Force vs displacement with markers denoting frames chosen for upscaling

    .. image:: GEOS_I43_01_sim38_damage.jpeg
       :alt: (b) Displacement field 

    GEOS I43.01 simulation force vs displacement

*********
Upscaling
*********

DNS results are homogenized using the Micromorphic Filter for 1, 24, 48, 264, 768, 1680, and 6912
filter domains.
Only the results for the 6912 element macroscale are presented.
Figure :numref:`{number} <GEOS_I43_01_filter_results>`
shows the homogenized displacement, stress, and damage fields.
These results generally agree with the DNS fields shown in Figures
:numref:`{number} <GEOS_I43_01_cylinder_displacements>`
and :numref:`{number} <GEOS_I43_01_cylinder_stress_damage>`, except
here the DNS results are clearly "smeared out."


 .. subfigure:: AA|BC
    :gap: 8px
    :subcaptions: below
    :name: GEOS_I43_01_filter_results
    :class-grid: outline
    :align: center

    .. image:: GEOS_I43_01_sim38_6912_displacement.jpeg
       :alt: (a) Homogenized displacement field

    .. image:: GEOS_I43_01_sim38_6912_cauchy33.jpeg
       :alt: (b) Homogenized Cauchy stress ZZ

    .. image:: GEOS_I43_01_sim38_6912_damage.jpeg
       :alt: (c) Homogenized damage

    Micromorphic Filter results for 6912 domains