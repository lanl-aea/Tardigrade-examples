.. _changelog:

#########
Changelog
#########

******************
0.2.0 (unreleased)
******************

******************
0.1.0 (2024-06-17)
******************

New Features
============
- Port over existing work for elastic cylinder. Includes hodgepodge of workflows for
  RatelF83 upscaling that needs to be improved and properly incorporated
  (:issue:`1`, :merge:`1`). By `Thomas Allard`_.
- For Abaqus elastic cylinder, build implicit dynamic workflow and fix aspects of the
  quasi-static workflow (:issue:`2`, :issue:`6`, :merge:`2`). By `Thomas Allard`_.
- Create a new "tard-ex-env" environment hosted on aea servers using the environment.yml
  file. This environment will only work with the html target for now and future
  implementation of workflows using the new micromorphic filter. Existing workflows
  may be run using a local conda environment built with the environment.txt file
  (:issue:`14`, :merge:`3`). By `Thomas Allard`_.
- Duplicate workflows for Abaqus elastic cylinder simulations (quasi-static and dynamic)
  and modify to use new micromorphic filter implementation. It is likely that previous
  workflows using old filter will be removed (:issue:`8`, :merge:`5`).
  By `Thomas Allard`_.
- Implement multidomain upscaling for Abaqus elastic cylinder DNS studies
  (:issue:`15`, :merge:`6`). By `Thomas Allard`_.
- Implement method of interpolating Abaqus nodal fields to integration points. New method
  compares well with previous (for displacements) and now velocity and acceleration
  fields may be upscaled (:issue:`17`, :merge:`7`). By `Thomas Allard`_.
- Post-processing filter results refactored to build specific plot and csv targets.
  Plots for deviatoric stress norms added. Csv files now generated to summarize statistics
  (mean, min, max, standard deviation) of various stress and deformation measures.
  Stress statistics for multi-domain workflows collected and summarized across
  parameter studies (:issue:`16`, :merge:`8`). By `Thomas Allard`_.
- Build Ratel locally and hook into SConstruct (:issue:`3`, :merge:`9`).
  By `Thomas Allard`_.
- Implement upscaling workflows for Ratel quasi-static elastic cylinder DNS
  through Micromorphic Filter including single and multiple filter domains
  (:issue:`4`, :merge:`10`). By `Thomas Allard`_.
- Implement upscaling workflow for Ratel F83 heterogeneous DNS through 
  Micromorphic Filter for multiple filter domains (including single)
  (:issue:`18`, :merge:`11`). By `Thomas Allard`_.
- Create new calibration script (calibrate_element.py) to calibrate micromorphic
  linear elasticity using averaged fields only for a specified element of the
  macroscale mesh filter domain. Updated multi domain workflows to use this script
  using a second nested parameter study to loop through each element.
  (:issue:`19`, :merge:`12`). By `Thomas Allard`_.
- Create script for parsing balance equation errors from Micromorphic Filter standard
  output which creates csv and plot files and added to workflows. Additional script
  for collecting output across multiple filtering domain studies
  (:issue:`23`, :merge:`13`). By `Thomas Allard`_.
- Added Ratel I41.02 elastic upscaling workflow (:issue:`26`, :merge:`14`).
  By `Thomas Allard`_.
- Added Tardigrade-MOOSE simulations to all workflows (:issue:`10`, :merge:`15`).
  By `Thomas Allard`_.
- Implemented better SConscript strategy to allow workflows to make use of common
  filter, calibration, and macroscale simulation steps (:issue:`20`, :merge:`17`).
  By `Thomas Allard`_.
- Add new joint probability distrbution plotting script and improve use of common
  SConscripts for upscaling workflows. Cleanup old meshes and DNS files. Add
  config file to specify program locations to be read by SConstruct
  (:merge:`19`). By `Thomas Allard`_.
- Migrate all steps for summarizing multi domain studies into a dedicated
  SConscript. Add CLI option to run this task (:issue:`33`, :merge:`20`).
  By `Thomas Allard`_.
- Migrate all "old" workflow associated with old Micromorphic Filter and rename
  all "new" workflows and scripts (:issue:`24`, :merge:`21`).
  By `Thomas Allard`_.
- Replace "options" argument in calibration and visualization scripts with explicit
  arguments for plotting, averaging, and calibration case. Remove unused "datacheck"
  target from all workflows (:issue:`26`, :merge:`22`). By `Thomas Allard`_.
- Migrate all SConscripts and workflows to a dedicated directory to declutter root
  (:issue:`34`, :merge:`23`). By `Thomas Allard`_.
- Implement --peta-data-copy local option to copy DNS files from the CU Peta library
  using the peta.py script throug SCP (:issue:`30`, :merge:`24`).
  By `Thomas Allard`_.
- Generate template meshes for users without access to Cubit and update workflows
  to handle this option (:issue:`29` and :merge:`25`). By `Thomas Allard`_.
- Add new script to make a copy of Micromorphic Filter XDMF results file where
  absolute paths are replaced with local paths to allow results to be visualized
  by Paraview without crashing (:issue:`35`, :merge:`27`). By `Thomas Allard`_.
- Add new studies for Abaqus and Ratel with clamped boundary conditions
  (:issue:`36`, :merge:`28`). By `Thomas Allard`_.
- Add new scripts and associated workflows to plot Ratel DNS and Tardigrade-MOOSE
  force vs displacement results. Additional summary script to summarize
  force vs displacement results together for multi domain workflows
  (:issue:`37`, :merge:`29`). By `Thomas Allard`_.
- Add new script to be used with "--config-software" to configure paths to
  various software. Replace previous config.yml file with a template.
  (:issue:`31`, :merge:`33`). By `Thomas Allard`_.
- Add new study for single filter domains "RVE" study for Ratel I41_02 elastic
  DNS (:issue:`9`, :merge:`35`). By `Thomas Allard`_.
- Add capability to apply "best" calibration results to Tardigrade-MOOSE simulations
  determined from peak values from a kernel density estimate to elements located on
  the boundary (:issue:`42`, :merge:`35`). By `Thomas Allard`_.
- Update Abaqus dynamic cylinder workflow to run basic macroscale simulation in
  Tardigrade-MOOSe (:issue:`46`, :merge:`38`). By `Thomas Allard`_.
- Add github workflow to deploy static documentation content to Pages for upcoming
  relesae (:issue:`44`, :merge:`41`). By `Thomas Allard`_.

Documentation
=============
- Port over existing documentation (:issue:`1`, :merge:`1`). By `Thomas Allard`_.
- Overhaul documentation for entire repository and document Abaqus elastic cylinder
  dynamic implicit and quasi-static workflows through direct numerical simulation
  (:issue:`2`, :issue:`6`, :merge:`2`). By `Thomas Allard`_.
- Update README.rst with environment activation instructions for local and AEA
  usage of new environment (:issue:`14`, :merge:`3`). By `Thomas Allard`_.
- Improve documentation for Abaqus dynamic elastic cylinder and include calculation
  of series convergence term to prescribe a load resulting in 1% strain
  (:issue:`8`, :merge:`5`). By `Thomas Allard`_.
- Provide basic instructions for building and using Ratel in WAVES worklow
  (:issue:`3`, :merge:`9`). By `Thomas Allard`_.
- Add docstrings for all relevant scripts used in current workflows and add to 
  API/CLI (:issue:`12`, :merge:`18`). By `Thomas Allard`_.
- Updating all documentation: improve uniaxial stress solutions, add pictures for
  Abaqus and Ratel elastic cylinder homogenization, improve formatting for software
  requirements and workflow overview, begin sections for all other upscaling
  studies (:issue:`13`, :merge:`26`). By `Thomas Allard`_.
- Document setup and use of all linked software. Add BSD-3 license file and add
  copyright and LANL code number O4375 to README
  (:issue:`9`, :merge:`35`). By `Thomas Allard`_.
- Document calibration workflow (:issue:`47`, :merge:`36`). By `Thomas Allard`_.
- Document macroscale simulation workflow (:issue:`48`, :merge:`37`).
  By `Thomas Allard`_.
- Create new image depicting the micromorphic reference and current configuration spaces.
  Document relevant micromorphic theory in appendix (:issue:`11`, :merge:`39`).
  By `Thomas Allard`_.
- Add and/or update all documentation for quasi-static verification, dynamic verification,
  and Ratel I41.02 upscaling studies (:issue:`43`, :merge:`40`). By `Thomas Allard`_.
