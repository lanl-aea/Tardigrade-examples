#########
|PROJECT|
#########

.. include:: project_brief.txt

.. raw:: latex

   \part{Project Setup}

.. toctree::
   :maxdepth: 2
   :caption: Micromorphic Upscaling Workflow

   workflow_overview
   workflow_micromorphic
   workflow_homogenization
   workflow_calibration
   workflow_macroscale

.. toctree::
   :maxdepth: 1
   :caption: Software Requirements

   software_overview
   software_installation
   software_usage

.. raw:: latex

   \part{Analyses}

.. toctree::
   :maxdepth: 1
   :caption: Quasi-static Verification

   Quasi_static_verification
   Ratel_elastic_cylinder
   Abaqus_elastic_cylinder
   Clamped_elastic_cylinder
   Tardigrade_MOOSE_convergence

.. toctree::
   :maxdepth: 1
   :caption: Dynamic Verification

   Dynamic_verification
   Abaqus_elastic_cylinder_dynamic_imp
   Tardigrade_MOOSE_dynamic_convergence

.. toctree::
   :maxdepth: 1
   :caption: Heterogeneous Upscaling

   Ratel_I41_02_elastic
   Ratel_I43_09_elastoplastic

.. toctree::
   :maxdepth: 1
   :caption: Micro-averaging Domain Studies

   Neper_cube
   Neper_cylinder

.. raw:: latex

   \part{Model Repository}

.. toctree::
   :maxdepth: 2
   :caption: Model Repository

   user
   api
   cli
   devops

.. raw:: latex

   \part{Help and Reference}

.. toctree::
   :caption: Help & Reference
   :maxdepth: 1

   release_philosophy
   changelog
   zreferences
   glossary

.. raw:: latex

   \part{Appendices}

.. toctree::
   :maxdepth: 2
   :caption: Appendix - Detailed Theory

   micromorphic_theory
   micromorphic_theory_kinematics
   micromorphic_theory_balance
   micromorphic_theory_constitutive
   micromorphic_theory_filter

.. only:: html

   ##################
   Indices and tables
   ##################

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
