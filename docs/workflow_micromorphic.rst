.. _workflow_micromorphic:

###################
Micromorphic Theory
###################

This section describes a small portion of micromorphic theory
This section describes a minimal portion of micromorphic theory as relevant
to understanding the inputs and outputs of the Micromorphic Filter
and Tardigarde-MOOSE. Further discussion is provided in :ref:`micromorphic_theory`.

.. _workflow_theory_filter:

**************************************
Homogenization via Micromorphic Filter
**************************************

The Micromorphic Filter calculates a variety of homogenized, macroscale
quantities using volume and surface integrals of DNS data over selected
micro-averaing domains. The following equations
(Eq. :math:numref:`{number} <homogenized_quantities>`)
define the macroscale
density, Cauchy stress, force, acceleration, couple stress, symmetric
micro stress, body force couple, and micro spin inertia. Further details
of the Micromorphic Filter and micromorphic quantities are provided by
Miller 2021 :cite:`miller_micromorphic_2021`, Miller et al. 2022
:cite:`miller_micromorphic_2022`, and a variety of other resources.

.. math::
    :label: homogenized_quantities

    \rho dv &\stackrel{\text{def}}{=} \int_{da}\rho^{\left(\alpha\right)}\,{dv^{\left(\alpha\right)}}

    \sigma_{ij} n_i da &\stackrel{\text{def}}{=} \int_{da} \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)} }

    \rho f_{i} dv &\stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} f_{i}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)} }

    \rho a_{i} dv &\stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} a_{i}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)} }

    m_{ijk} n_i da &\stackrel{\text{def}}{=} \int_{da} \sigma_{ij}^{\left(\alpha\right)} \xi_k n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)} }

    s_{ij} dv &\stackrel{\text{def}}{=} \int_{dv} \sigma_{ij}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)} }

    \rho l_{ij} dv &\stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} f_{i}^{\left(\alpha\right)} \xi_j\,{dv^{\left(\alpha\right)} }

    \rho \omega_{ij} dv &\stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \ddot{\xi_i} \xi_j \,{dv^{\left(\alpha\right)}}

..
   TODO: update description for new version of Micromorphic Filter (no surface integration?) provide description of how the new Micromorphic Filter calculates Cauchy and Higher Order Stress

The Micromorphic Filter also determines the macroscale deformation gradient,
:math:`F_{iI}`, the micro-deformation tensor, :math:`\chi_{iJ}`, and the
gradient of the micro-deformation tensor, :math:`\chi_{iJ,K}`. With these
terms, deformation measures may be determined including the Green-Lagrange
strain (:math:`E_{ij}`), Euler-Almansi strain (:math:`e_{ij}`), micro-strain
(:math:`\mathcal{E}_{IJ}`), and micro-deformation gradient (:math:`\Gamma_{IJK}`)
shown in equation :math:numref:`{number} <deformation_measures>`.

.. math::
    :label: deformation_measures

    E_{IJ} &= \frac{1}{2} \left( F_{iI} F_{iJ} - \delta_{IJ}\right)

    e_{ij} &= F_{Ii}^{-1} E_{IJ} F_{Jj}^{-1}

    \mathcal{E}_{IJ} &= F_{iI} \chi_{iJ} - \delta_{IJ}

    \Gamma_{IJK} &= F_{iI} \chi_{iJ,K}

****************
Tardigrade-MOOSE
****************

The relevant balance equations (in the current configuration) to describe
a determinant system with 12 unknowns may be defined for the balance of
linear momentum and the balance of the first moment of momentum.

.. math::
    :label: balance_equations

	\sigma_{lk,l} + \rho \left( f_k - a_k \right) &= 0

	\sigma_{mk} - s_{mk} + m_{lkm,l} + \rho \left( l_{mk} - \omega_{mk}\right) &= 0

Tardigrade-MOOSE solves these equations to find the unknown displacements,
:math:`\mathbf{u}`, and micro-displacements, :math:`\mathbf{\Phi}`, where
:math:`\mathbf{\chi} = \mathbf{I} + \mathbf{\Phi}`.

********************************
Micromorphic Constitutive Models
********************************

.. include:: workflow_constitutive.txt