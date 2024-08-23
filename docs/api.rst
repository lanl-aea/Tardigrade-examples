.. _sphinx_api:

#############
|PROJECT| API
#############

config_software.py
==================

.. automodule:: model_package.config_software
    :members:
    :show-inheritance:
    :synopsis: Configure software paths in a YAML file

peta.py
=======

.. automodule:: model_package.peta
    :members:
    :show-inheritance:
    :synopsis: Copy DNS results from the CU Peta library to the output directory

xdmf_reader_tools.py
====================

.. automodule:: model_package.xdmf_reader_tools
    :members:
    :show-inheritance:
    :synopsis: Functions for reading XDMF files

DNS_Abaqus.build_dynamic_elastic_cylinder.py
============================================

.. automodule:: model_package.DNS_Abaqus.build_dynamic_elastic_cylinder
    :members:
    :show-inheritance:
    :synopsis: Create an Abaqus model of an elastic cylinder under dynamic compression

DNS_Abaqus.build_elastic_cylinder.py
====================================

.. automodule:: model_package.DNS_Abaqus.build_elastic_cylinder
    :members:
    :show-inheritance:
    :synopsis: Create an Abaqus model of an elastic cylinder under static compression

DNS_Abaqus.dynamic_analytical_comparison.py
===========================================

.. automodule:: model_package.DNS_Abaqus.dynamic_analytical_comparison
    :members:
    :show-inheritance:
    :synopsis: Plot dynamic Abaqus results against an analytical solution

DNS_Abaqus.extract_history.py
=============================

.. automodule:: model_package.DNS_Abaqus.extract_history
    :members:
    :show-inheritance:
    :synopsis: Plot Abaqus history output for force versus displacement

DNS_Abaqus.extract_frames.py
============================

.. automodule:: model_package.DNS_Abaqus.extract_frames
    :members:
    :show-inheritance:
    :synopsis: Extracts 3D field output from a completed Abaqus simulation to save as 2D image

DNS_Abaqus.modify_input.py
==========================

.. automodule:: model_package.DNS_Abaqus.modify_input
    :members:
    :show-inheritance:
    :synopsis: Modify Abaqus input file to output 'COORD' at integration points

DNS_Abaqus.ODBextract_to_XDMF.py
================================

.. automodule:: model_package.DNS_Abaqus.ODBextract_to_XDMF
    :members:
    :show-inheritance:
    :synopsis: Convert Abaqus DNS results to XDMF format

DNS_Ratel.build_options_file.py
===============================

.. automodule:: model_package.DNS_Ratel.build_options_file
    :members:
    :show-inheritance:
    :synopsis: Write Ratel options file

DNS_Ratel.plot_force_displacement.py
====================================

.. automodule:: model_package.DNS_Ratel.plot_force_displacement
    :members:
    :show-inheritance:
    :synopsis: Process force-displacement from Ratel DNS results

DNS_Ratel.vtk_to_xdmf.py
========================

.. automodule:: model_package.DNS_Ratel.vtk_to_xdmf
    :members:
    :show-inheritance:
    :synopsis: Convert Ratel DNS results to XDMF format

Filter.bounds_from_DNS.py
=========================

.. automodule:: model_package.Filter.bounds_from_DNS
    :members:
    :show-inheritance:
    :synopsis: Create a csv containing the extents of a DNS file

Filter.build_filter_config.py
=============================

.. automodule:: model_package.Filter.build_filter_config
    :members:
    :show-inheritance:
    :synopsis: Write the configuration file for the Micromorphic Filter

Filter.collect_multi_domain_errors.py
=====================================

.. automodule:: model_package.Filter.collect_multi_domain_errors
    :members:
    :show-inheritance:
    :synopsis: Collect balance equation errors across filter domain studies

Filter.collect_multi_domain_stats.py
====================================

.. automodule:: model_package.Filter.collect_multi_domain_stats
    :members:
    :show-inheritance:
    :synopsis: Collect statistics of a homogenized micromorphic quantity across filter domain studies

Filter.force_bounds.py
======================

.. automodule:: model_package.Filter.force_bounds
    :members:
    :show-inheritance:
    :synopsis: Create a csv file containing information for a bounding box encompassing all DNS points

Filter.parse_balance_errors.py
==============================

.. automodule:: model_package.Filter.parse_balance_errors
    :members:
    :show-inheritance:
    :synopsis: Parse balance equation errors from Micromorphic Filter standard output

Filter.run_micromorphic_filter.py
=================================

.. automodule:: model_package.Filter.run_micromorphic_filter
    :members:
    :show-inheritance:
    :synopsis: Run the Micromorphic Filter

Filter.single_macroscale.py
===========================

.. automodule:: model_package.Filter.single_macroscale
    :members:
    :show-inheritance:
    :synopsis: Write a single macroscale domain file for the Micromorphic Filter

Filter.visualize_results.py
================================

.. automodule:: model_package.Filter.visualize_results
    :members:
    :show-inheritance:
    :synopsis: Post-process Micromorphic Filter output

Filter.xdmf_local_paths.py
==========================

.. automodule:: model_package.Filter.xdmf_local_paths
    :members:
    :show-inheritance:
    :synopsis: Create a copy of an XDMF file with absolute H5 paths replaced with relative paths

Filter.xdmf_tomfoolery.py
=========================

.. automodule:: model_package.Filter.xdmf_tomfoolery
    :members:
    :show-inheritance:
    :synopsis: Modify an XDMF file by combining elements from separate 'blocks'

Calibrate.build_calibration_map.py
==================================

.. automodule:: model_package.Calibrate.build_calibration_map
    :members:
    :show-inheritance:
    :synopsis: Create a yaml file to map calibration results

Calibrate.calibrate_element.py
==============================

.. automodule:: model_package.Calibrate.calibrate_element
    :members:
    :show-inheritance:
    :synopsis: Calibrate micromorphic linear elasticity for a single filter element

Calibrate.joint_probability_distributions.py
============================================

.. automodule:: model_package.Calibrate.joint_probability_distributions
    :members:
    :show-inheritance:
    :synopsis: Create a joint probability distribution plot to summarize calibration results

Calibrate.return_minimum_smith_constraint.py
============================================

.. automodule:: model_package.Calibrate.return_minimum_smith_constraint
    :members:
    :show-inheritance:
    :synopsis: Evaluate the 13 Smith constraints for micromorphic linear elasticity and return the minimum value

Calibrate.summarize_calibration_results.py
==========================================

.. automodule:: model_package.Calibrate.summarize_calibration_results
    :members:
    :show-inheritance:
    :synopsis: Summarize results of parameter calibration

Calibrate.summarize_calibration_results_ignore_boundary.py
==========================================================

.. automodule:: model_package.Calibrate.summarize_calibration_results_ignore_boundary
    :members:
    :show-inheritance:
    :synopsis: Summarize results of parameter calibration while ignoring elements on the z-boundary"

Tardigrade_MOOSE.add_element_blocks_to_mesh.py
==============================================

.. automodule:: model_package.Tardigrade_MOOSE.add_element_blocks_to_mesh
    :members:
    :show-inheritance:
    :synopsis: Take an existing exodus mesh, add element blocks for each element, save with new name

Tardigrade_MOOSE.build_dynamic_Tardigrade_input_deck.py
=======================================================

.. automodule:: model_package.Tardigrade_MOOSE.build_dynamic_Tardigrade_input_deck
    :members:
    :show-inheritance:
    :synopsis: Write a Tardigrade-MOOSE input file for dynamic simulation

Tardigrade_MOOSE.build_Tardigrade_input_deck.py
===============================================

.. automodule:: model_package.Tardigrade_MOOSE.build_Tardigrade_input_deck
    :members:
    :show-inheritance:
    :synopsis: Write Tardigrade-MOOSE input file

Tardigrade_MOOSE.cylinder_from_bounds.py
========================================

.. automodule:: model_package.Tardigrade_MOOSE.cylinder_from_bounds
    :members:
    :show-inheritance:
    :synopsis: Create a cylinder mesh from the bounds of a DNS file

Tardigrade_MOOSE.finite_stVK_calculation.py
============================================

.. automodule:: model_package.Tardigrade_MOOSE.finite_stVK_calculation
    :members:
    :show-inheritance:
    :synopsis: Solution for uniaxial stress of a cylinder for finite deformation using the St. Venant-Kirchhoff elasticity model

Tardigrade_MOOSE.plot_dynamic_displacement.py
=============================================

.. automodule:: model_package.Tardigrade_MOOSE.plot_dynamic_displacement
    :members:
    :show-inheritance:
    :synopsis: Process displacement vs time from Tardigrade-MOOSE results

Tardigrade_MOOSE.plot_force_displacement.py
===========================================

.. automodule:: model_package.Tardigrade_MOOSE.plot_force_displacement
    :members:
    :show-inheritance:
    :synopsis: Process force-displacement from Tardigrade-MOOSE results

Tardigrade_MOOSE.summarize_micro_macro_force_displacements.py
=============================================================

.. automodule:: model_package.Tardigrade_MOOSE.summarize_micro_macro_force_displacements
    :members:
    :show-inheritance:
    :synopsis: Process force-displacement from Tardigrade-MOOSE results

Tardigrade_MOOSE.write_elastic_material_card.py
===============================================

.. automodule:: model_package.Tardigrade_MOOSE.write_elastic_material_card
    :members:
    :show-inheritance:
    :synopsis: Write elastic Tardigrade-MOOSE input card (.yml)