.. _workflow_overview:

#################
Workflow Overview
#################

.. include:: README.txt
   :start-after: upscaling-workflow-description-start-do-not-remove
   :end-before: upscaling-workflow-description-end-do-not-remove

Figure :numref:`{number} <upscaling-flowchart>` depicts the upscaling workflow.
First, DNS of the "microscale" is performed, processed, and homogenized
through the Micromorphic Filter. Next, homogenized quantities are calibrated to
constitutive models (either classical or micromorphic). Finally, the
calibrated constitutive model(s) are applied to a macroscale simulation
in Tardigarde-MOOSE.

.. figure:: upscaling_flowchart.png
   :name: upscaling-flowchart
   :align: center
   :width: 50%

   Upscaling workflow

Upscaling workflows are typically partitioned into 4 separate stages
(DNS, filter, calibrate, and macroscale simulation) using the following command
line arguments:

   .. code-block:: console

      $ scons study_name
      $ scons study_name --filter
      $ scons study_name --calibrate
      $ scons study_name --macro

Further details and options for controlling stages of upscaling analyses are provided in
the following sections and documentation of relevant upscaling studies.