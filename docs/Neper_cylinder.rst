.. _neper_cylinder:

####################
Neper Cylinder Study
####################

.. figure:: neper_cylinder_grain_geometry.jpg
   :name: neper_cyliner_grain_geometry
   :align: center
   :width: 40%

   Example of meshed Abaqus cylinder geometry generated from grain growth tesselation results from Neper

To run this study, use the following commands.

   .. code-block:: bash

      $ scons neper_cylinder --solve-cpus=8
      $ scons neper_cylinder --filter

Modifications to the Micromorphic Filter are needed to perform this analysis.
Documentation currently has not been created for this study or required
modifications.
