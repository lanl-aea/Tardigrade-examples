.. _neper_cube:

################
Neper Cube Study
################

.. figure:: neper_cube_grain_geometry.jpg
   :name: neper_cube_grain_geometry
   :align: center
   :width: 50%

   Example of meshed Abaqus cube geometry generated from grain growth tesselation results from Neper

To run this study, use the following commands.

   .. code-block:: bash

      $ scons neper_cube --solve-cpus=8
      $ scons neper_cube --filter

Modifications to the Micromorphic Filter are needed to perform this analysis.
Documentation currently has not been created for this study or required
modifications.
