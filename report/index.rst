#####################################
|PROJECT|: Tardigrade-examples Report
#####################################

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
   :maxdepth: 2
   :caption: Quasi-static Verification

   Quasi_static_verification
   Ratel_elastic_cylinder
   Abaqus_elastic_cylinder
   Clamped_elastic_cylinder
   Tardigrade_MOOSE_convergence

.. toctree::
   :maxdepth: 2
   :caption: Dynamic Verification
   
   Dynamic_verification
   Abaqus_elastic_cylinder_dynamic_imp
   Tardigrade_MOOSE_dynamic_convergence

.. toctree::
   :maxdepth: 2
   :caption: Heterogeneous Upscaling

   Ratel_I41_02_elastic

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

   \part{Appendices}

.. toctree::
   :maxdepth: 2
   :caption: Appendix - Detailed Theory

   micromorphic_theory
   micromorphic_theory_kinematics
   micromorphic_theory_balance
   micromorphic_theory_constitutive
   micromorphic_theory_filter
