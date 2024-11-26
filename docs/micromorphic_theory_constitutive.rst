.. _micromorphic_theory_constitutive:

######################
Constitutive Equations
######################

******************************
Micromorphic Linear Elasticity
******************************

..
   TODO: write this section

Eringen and Suhubi :cite:`eringen_nonlinear_1964` proposed unique sets of elastic deformation measures as functions
of the macroscopic deformation gradient :math:`\mathbf{F}` and the micro-deformation tensor :math:`\mathbf{\chi}`.
These deformation measures include the left Cauchy-Green deformation tensor
:math:`\mathbf{\mathcal{C}}`, a micro-deformation tensor :math:`\mathbf{\Psi}`,
and the micro-deformation gradient :math:`\mathbf{\Gamma}`, as calculated in
Eq. :math:numref:`{number} <deformation_measures_1>`.
From the deformation measures introduced in Eq. :math:numref:`{number} <deformation_measures_1>`,
two strain measures are defined:
the Green-Lagrange strain, :math:`\mathbf{E}`, and the micro-strain, :math:`\mathbf{\mathcal{E}}`,
calculated according to equations :math:numref:`{number} <GL_strain>` and
:math:numref:`{number} <micro_strain>`, respectively.

A quadratic form for the Helmholtz free energy is introduced in equation :math:numref:`{number} <elastic_helmholtz>`
using :math:`\mathbf{E}`, :math:`\mathbf{\mathcal{E}}`, and :math:`\mathbf{\Gamma}`,

.. math::
   :label: elastic_helmholtz

   \left(\rho_{0}\psi\right)^e &= \frac{1}{2}E_{KL}A_{KLMN}E_{MN} + \frac{1}{2}\mathcal{E}_{KL}B_{KLMN}\mathcal{E}_{MN}

   &+ \frac{1}{2}\Gamma_{KLM}C_{KLMNPQ}\Gamma_{NPQ} + E_{KL}D_{KLMN}\mathcal{E}_{MN},

where :math:`\mathbf{A}`, :math:`\mathbf{B}`, :math:`\mathbf{C}`, and :math:`\mathbf{D}`
are the elastic material moduli tensors, which may be calculated according to equations
:math:numref:`{number} <A>`, :math:numref:`{number} <B>`, :math:numref:`{number} <C>`, and
:math:numref:`{number} <D>`, respectively.

.. math::
   :label: A

   A_{KLMN} = \lambda\delta_{KL}\delta_{MN} + \mu \left(\delta_{KM}\delta_{LN} + \delta_{KN}\delta_{LM} \right)

.. math::
   :label: B

   B_{KMLN} = \left(\eta - \tau\right) \delta_{KL}\delta_{MN} + \left(\kappa - \sigma\right) \delta_{KM}\delta_{LN} + \left(\nu - \sigma \right) \delta_{KN}\delta_{LM}

.. math::
   :label: C

   C_{KLMNPQ} &=
      \tau_1 \left(\delta_{KL}\delta_{MN}\delta_{PQ} + \delta_{KQ}\delta_{LM}\delta_{NP}\right) +
      \tau_2 \left(\delta_{KL}\delta_{MP}\delta_{NQ} + \delta_{KM}\delta_{LQ}\delta_{NP}\right)

   &+
      \tau_3 \delta_{KL}\delta_{MQ}\delta_{NP} +
      \tau_4 \delta_{KN}\delta_{LM}\delta_{PQ} +
      \tau_5 \left(\delta_{KM}\delta_{LN}\delta_{PQ} + \delta_{KP}\delta_{LM}\delta_{NQ}\right)

   &+
      \tau_6 \delta_{KM}\delta_{LP}\delta_{NQ} +
      \tau_7 \delta_{KN}\delta_{LP}\delta_{MQ} +
      \tau_8 \left(\delta_{KP}\delta_{LQ}\delta_{MN} + \delta_{KQ}\delta_{LN}\delta_{MP}\right)

   &+
      \tau_9 \delta_{KN}\delta_{LQ}\delta_{MP} +
      \tau_10 \delta_{KP}\delta_{LN}\delta_{MQ} +
      \tau_{11} \delta_{KQ}\delta_{LP}\delta_{MN}

.. math::
   :label: D

   D_{KLMN} = \tau\delta_{KL}\delta_{MN} + \sigma \left(\delta_{KN}\delta_{LM} + \delta_{LN}\delta_{KM}\right)

Equations :math:numref:`{number} <A>`, :math:numref:`{number} <B>`, :math:numref:`{number} <C>`,
and :math:numref:`{number} <D>` introduce the 18 parameters for the linear elastic micromorphic constitutive model.
Calibration will seek to determine an admissible set of parameters that best describes the homogenized DNS response.
The 18 parameters are identified as :math:`\lambda`, :math:`\mu`, :math:`\eta`, :math:`\tau`, :math:`\kappa`,
:math:`\nu`, :math:`\sigma`, and :math:`\tau_{1}` through :math:`\tau_{11}`.
It should be noted that parameters :math:`\tau_{1}` through :math:`\tau_{11}` are only present in equation
:math:numref:`{number} <C>` which will be used to relate the micro-deformation gradient
:math:`\mathbf{\Gamma}` to higher order stress effects.

The stresses (second Piola-Kirchhoff stress :math:`\mathbf{S}`, symmetric micro-stress :math:`\mathbf{\Sigma}`,
and higher order stress :math:`\mathbf{M}`) may be derived from the Helmholtz free energy
as follows:

.. math::
   :label: PK2_1

   S_{IK} = 2 \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \mathcal{C}_{IJ}}
      + \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \Psi_{IQ}} \Psi_{KQ} \mathcal{C}_{JK}^{-1}
      + \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{IQK}} \Gamma_{SQK} \mathcal{C}_{JS}^{-1}

.. math::
   :label: SIGMA_1

   \Sigma_{IJ} = 2 \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \mathcal{C}_{IJ}}
      + 2 symm\left[ \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \Psi_{IQ}} \Psi_{KQ} \mathcal{C}_{JK}^{-1}
      + \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{IQK}} \Gamma_{SQK} \mathcal{C}_{JS}^{-1} \right]

.. math::
   :label: M_1

   M_{IJK} = \frac{\partial\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{JKI}}.

By taking the relevant partial derivatives of the elastic Helmholtz free energy function,
equations :math:numref:`{number} <PK2_1>`, :math:numref:`{number} <SIGMA_1>`, and :math:numref:`{number} <M_1>`
may be evaluated as follows:

.. math::
   :label: PK2_2

   S_{IJ} =& A_{IJKL}E_{KL} + D_{IJKL} \mathcal{E}_{KL} + \left\{B_{IQKL}\mathcal{E}_{KL}
      + E_{KL}D_{KLIQ}\right\}\left(\mathcal{E}_{RQ} + \delta_{RQ}\right)
      \left(\mathcal{C}_{RJ}\right)^{-1}

   &+ C_{IQRLMN} \Gamma_{LMN}  \left(\mathcal{C}_{SJ}\right)^{-1} \Gamma_{SQR}

.. math::
   :label: SIGMA_2

   \Sigma_{IJ} =&  A_{IJKL}E_{KL} + D_{IJKL} \mathcal{E}_{KL}

   &+ 2symm \left( \left\{B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right\} \left( \mathcal{E}_{RQ}
      + \delta_{RQ}\right) \left(\mathcal{C}_{RJ}\right)^{-1}\right)

   &+ 2symm \left(C_{IQRLMN} \Gamma_{LMN} \Gamma_{SQR} \left(\mathcal{C}_{SJ}\right)^{-1} \right)

.. math::
   :label: M_2

   M_{IJK} = C_{JKILMN} \Gamma_{LMN}.

Finaly, the elastic moduli tensors (equations :math:numref:`{number} <A>`, :math:numref:`{number} <B>`,
:math:numref:`{number} <C>`, and :math:numref:`{number} <D>`) may be substituted into
equations :math:numref:`{number} <PK2_2>`, :math:numref:`{number} <SIGMA_2>`, and :math:numref:`{number} <M_2>`
to express the stresses as functions of the 18 elasticity parameters, resulting in:

.. math::
   :label: PK2

   S_{IJ} = \left(\lambda^* + \tau^*\right) E_{MM} \delta_{IJ}
      + 2\left(\mu^* + \sigma^*\right) E_{IJ}
      + \eta^* \mathcal{E}_{MM} \delta_{IJ}
      + \kappa^* \mathcal{E}_{IJ}
      + \nu^* \mathcal{E}_{JI}

.. math::
   :label: SIGMA

   \Sigma_{IJ} &= \left(\lambda^* + 2\tau^*\right) E_{MM} \delta_{IJ}
      + 2\left(\mu^* + 2\sigma^*\right) E_{IJ}
      + \left(2\eta^* - \tau^*\right) \mathcal{E}_{MM} \delta_{IJ}

   &+ \left(\nu^* + \kappa^* - \sigma\right)
      \left(\mathcal{E}_{IJ} + \mathcal{E}_{JI}\right)

.. math::
   :label: M

   M_{IJK} &= \tau_1^* \left(\delta_{JK}\Gamma_{IPP} + \delta_{KI} \Gamma_{PPJ}\right)
      + \tau_2^*  \left(\delta_{JK}\Gamma_{NIN} + \delta_{JI} \Gamma_{PPK}\right)

   &+ \tau_3^* \delta_{JK} \Gamma_{NNI}
      + \tau_4^* \delta_{KI} \Gamma_{JPP}
      + \tau_5^* \left(\delta_{JI}\Gamma_{KPP} + \delta_{KI} \Gamma_{NJN}\right)

   &+ \tau_6^* \delta_{JI} \Gamma_{NKN}
      + \tau_7^* \Gamma_{JKI}
      + \tau_8^* \left(\Gamma_{IJK} \Gamma_{KIJ}\right)
      + \tau_9^* \Gamma_{JIK}

      &+ \tau_{10}^* \Gamma_{KJI}
      + \tau_{11}^* \Gamma_{IJK},

which are the same as shown in equation :math:numref:`{number} <constitutive_case_4>`
(although here different indices are used) discussed in
the :ref:`workflow_constitutive_linear_elasticity` section while describing the
micromorphic upscaling workflow.

.. _smith_conditions:

Smith Conditions
================

Refer to the :ref:`linear_elastic_constraints` section for the discussion of the Smith
conditions.

..
   TODO: Describe the Smith conditions

******************************
Micromorphic Elasto-plasticity
******************************

Kinematics
==========

A multiplicative decomposition of the macro deformation gradient and micro deformation
tensor is used to separate elastic and plastic effects.

.. math::

   F_{iI} &= F^{e}_{i\bar{I}} F^p_{\bar{I}I}

   \chi_{iI} &= \chi^{e}_{i\bar{I}} \chi^p_{\bar{I}I}

Figure :numref:`{number} <FeFp_configurations>` (borrowed from Miller 2021 :cite:`miller_micromorphic_2021`
Figure 3.1) shows the effect of the multiplicative decomposition and resulting
configuration spaces.
As before, the reference configuration maps to the current configuration through
:math:`\mathbf{F}` and :math:`\mathbf{\chi}`.
A stress free, intermediate configuration (:math:`\bar{B}`) is introduced which maps from the reference
configuration through :math:`\mathbf{F^p}` and :math:`\mathbf{\chi^p}`.
The overbar notation :math:`\bar{.}` denotes quantities in the intermediate configuration.
Figure :numref:`{number} <FeFp_configurations>` show the independent actions of macro
deformation gradients and micro deformation tensors.

..
   TODO: Make my own version of the folowing figure

.. figure:: FeFp_configurations.jpg
   :name: FeFp_configurations
   :align: center
   :width: 80%

   Configuration spaces of the multiplicative decomposition of :math:`\mathbf{F}` and :math:`\mathbf{\chi}`

The elastic deformation measures from equations :math:numref:`{number} <deformation_measures_1>`,
:math:numref:`{number} <GL_strain>`, and :math:numref:`{number} <micro_strain>` are now defined in the
intermediate configuration as

.. math::
   :label: deformation_measures_intermediate

   \bar{\mathcal{C}}^e_{\bar{I}\bar{J}} &= F^e_{i\bar{I}} F^e_{i\bar{J}}

   \bar{\Psi}^e_{\bar{I}\bar{J}} &= F^e_{i\bar{I}} \xi^e_{i\bar{J}}

   \bar{\Gamma}^e_{\bar{I}\bar{J}K} &= F^e_{i\bar{I}} \xi^e_{i\bar{J},\bar{K}}

   \bar{E}^e_{\bar{I}\bar{J}} &= \frac{1}{2} \left( \mathcal{C}^e_{\bar{I}\bar{J}} - \delta_{\bar{I}\bar{J}} \right)

   \bar{\mathcal{E}}^e_{\bar{I}\bar{J}} &= \Psi^e_{\bar{I}\bar{J}} - \delta_{\bar{I}\bar{J}}

Further details of micromorphic elastoplasticity may be found in Chapter 3 of Miller 2021 :cite:`miller_micromorphic_2021`.

Deviatoric Stress Measures
==========================

The deviatoric parts of the Cauchy, symmetric micro-, and higher-order stresses may be
defined in the current configuration as

.. math::

   dev\left(\sigma_{ij}\right) & \stackrel{\text{def}}{=} \sigma_{ij} - \frac{1}{3}\sigma_{kk}\delta_{ij}

   dev\left(s_{ij}\right) & \stackrel{\text{def}}{=} s_{ij} - \frac{1}{3} s_{kk}\delta_{ij}

   dev\left(m_{ijk}\right) & \stackrel{\text{def}}{=} m_{ijk} - \frac{1}{3} m_{llk}\delta_{ij}.


The Cauchy, micro-, and higher-order pressures may be defined in the current configuration as

.. math::

   p^u &\stackrel{\text{def}}{=} \frac{1}{3}\sigma_{kk}
        = \frac{1}{3 J^e} F^e_{i\bar{I}}\bar{S}_{\bar{I}\bar{J}}F^e_{i\bar{J}}
		= \frac{1}{3 J^e}\bar{C}^e_{\bar{I}\bar{J}}\bar{S}_{\bar{I}\bar{J}}

   p^{\chi} &= \frac{1}{3} s_{kk}
		= \frac{1}{3 J^e}\bar{C}^e_{\bar{I}\bar{J}}\bar{\Sigma}_{\bar{I}\bar{J}}

   p^{\nabla\chi}_k &= \frac{1}{3} m_{llk}
		= \frac{1}{3 J^e}\bar{M}_{\bar{I}\bar{J}\bar{K}}\bar{C}^e_{\bar{I}\bar{J}}\chi^e_{k\bar{K}}

where :math:`\bar{C}^e_{\bar{I}\bar{J}}=F^e_{i\bar{I}}F^e_{j\bar{J}}`.
These pressures may be defined in the intermediate configuration as

.. math::

   \bar{p}^u & =\frac{1}{3}\bar{C}^e_{\bar{I}\bar{J}}\bar{S}_{\bar{I}\bar{J}}

   \bar{p}^{\chi} &=\frac{1}{3}\bar{C}^e_{\bar{I}\bar{J}}\bar{\Sigma}_{\bar{I}\bar{J}}

   \bar{p}^{\nabla\chi}_{\bar{K}} &= \frac{1}{3}\bar{C}^e_{\bar{I}\bar{J}}\bar{M}_{\bar{I}\bar{J}\bar{K}}.

Using these terms, the deviatoric parts of the Second Piola-Kirchhoff, symmetric micro-
and higher-order stresses may be written in the intermediate configuration as

.. math::

   dev\left(\bar{S}_{\bar{I}\bar{J}}\right) &= \bar{S}_{\bar{I}\bar{J}}
       - \bar{p}^u\left(\bar{C}^e_{\bar{I}\bar{J}}\right)^{-1}

   dev\left(\bar{\Sigma}_{\bar{I}\bar{J}}\right) &= \bar{\Sigma}_{\bar{I}\bar{J}}
       - \bar{p}^{\chi}\left(\bar{C}^e_{\bar{I}\bar{J}}\right)^{-1}

   dev\left(\bar{M}_{\bar{I}\bar{J}\bar{K}}\right) &= \bar{M}_{\bar{I}\bar{J}\bar{K}}
       - \bar{p}_{\bar{K}}^{\nabla\chi}\left(\bar{C}^e_{\bar{I}\bar{J}}\right)^{-1}.

Helmholtz Free Energy Function
==============================

A micromorphic, linear isotropic, Drucker-Prager elastoplasticity is considered.
The total Helmholtz free energy per unit mass in the intermediate configuration may
be expressed as the addition of the elastic free energy function
:math:`\left( \bar{\rho} \bar{\psi}\right)^e`
(introduced in Eq. :math:numref:`{number} <elastic_helmholtz>`
in the reference configuration) and plastic free energy
function :math:`\left( \bar{\rho} \bar{\psi}\right)^p`.

.. math::

   \left(\bar{\rho} \bar{\psi}\right) \stackrel{\text{def}}{=}
       \left(\bar{\rho} \bar{\psi}\right)^e + \left(\bar{\rho} \bar{\psi}\right)^p

A quadratic form of the plastic Helmholtz energy function is defined as a function
of strain-like internal state variables (ISVs) and hardening moduli.

.. math::
   :label: plastic_helmholtz

   \left(\bar{\rho} \bar{\psi}\right)^p \stackrel{\text{def}}{=}
       \frac{1}{2}\bar{H}^u \left(\bar{Z}^u\right)^2
       + \frac{1}{2}\bar{H}^{\chi} \left(\bar{Z}^{\chi}\right)^2
       + \frac{1}{2}\bar{Z}^{\chi}_{,\bar{I}} \bar{H}_{\bar{I}\bar{J}}^{\nabla\chi} \bar{Z}^{\chi}_{,\bar{J}}

where :math:`\bar{Z}^u` and :math:`\bar{Z}^{\chi}` are scalars.

As will be discussed in the proceeding section, :ref:`yield_surfaces`, three yield surfaces are
defined for macro- (:math:`u`), micro- (:math:`\chi`), and micro-gradient plasticity (:math:`\nabla\chi`)
with associated strain-like ISVs (:math:`\bar{Z}^u`, :math:`\bar{Z}^{\chi}`,
and :math:`\bar{Z}^{\nabla \chi}`) and hardening moduli
(:math:`\bar{H}^u`, :math:`\bar{H}^{\chi}`, and :math:`\bar{H}^{\nabla \chi}`).
Stress-like ISVs are defined as

.. math::
   :label: micromorphic_ISVs

   \bar{Q}^u &\stackrel{\text{def}}{=}
       \frac{\partial \left( \bar{\rho} \bar{\psi} \right)}{\partial \bar{Z}^u_{\bar{I}}}
       = \bar{H}^u\bar{Z}^u

   \bar{Q}^{\chi} &\stackrel{\text{def}}{=}
       \frac{\partial \left( \bar{\rho} \bar{\psi} \right)}{\partial \bar{Z}^{\chi}_{\bar{I}}}
       = \bar{H}^{\chi}\bar{Z}^{\chi}

   \bar{Q}^{\nabla\chi}_{\bar{I}} &\stackrel{\text{def}}{=}
       \frac{\partial \left( \bar{\rho} \bar{\psi} \right)}{\partial \bar{Z}^{\nabla\chi}_{\bar{I},\bar{J}}}
       = \bar{H}^{\nabla\chi}_{\bar{I}\bar{J}}\bar{Z}^{\chi}_{,\bar{J}}.

For isotropy, the micro-gradient hardening moduli are simplified to
:math:`\bar{H}^{\nabla\chi}_{\bar{I}\bar{J}}=\bar{H}^{\nabla\chi} \delta_{\bar{I}\bar{J}}`
which reduces the stress-like ISV to :math:`\bar{Q}^{\chi}_{\bar{I}}=\bar{H}^{\nabla\chi}\bar{Z}^{\chi}_{,\bar{I}}`.

.. _yield_surfaces:

Yield Surfaces
==============

The Drucker-Prager yield functions may be defined for macro- (:math:`u`),
micro- (:math:`\chi`), and micro-gradient plasticity (:math:`\nabla\chi`) as

.. math::

   \bar{F}^u\left(\bar{\mathbf{S}},\bar{c}^u\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{S}}\right)|| - \left(A^{u,\phi}\bar{c}^u
       - B^{u,\phi} \bar{p}^u\right) \leq 0

   \bar{F}^{\chi}\left(\bar{\mathbf{\Sigma}},\bar{c}^{\chi}\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{\Sigma}}\right)|| - \left(A^{{\chi},\phi}\bar{c}^{\chi}
       - B^{{\chi},\phi} \bar{p}^{\chi}\right) \leq 0

   \bar{F}_{\bar{K}}^{\nabla\chi}\left(\bar{\mathbf{M}},\bar{\mathbf{c}}^{\nabla\chi}\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{M}}\right)||_{\bar{K}} - \left(A^{{\nabla\chi},\phi}\bar{c}_{\bar{K}}^{\nabla\chi}
       - B^{{\nabla\chi},\phi} \bar{p}_{\bar{K}}^{\nabla\chi}\right) \leq 0

where :math:`\bar{c}^u`, :math:`\bar{c}^{\chi}`, and :math:`\bar{c}^{\nabla\chi}` are the cohesion terms
for each of the respective yield surface.
Note that there are 3 yield surface for micro-gradient plasticity.
The invariants of the stress measures are defined as


.. math::

   ||dev\left(\bar{\mathbf{S}}\right)|| &= \sqrt{\left(dev \left( \bar{S}_{\bar{I}\bar{J}}\right)\right)
       : \left(dev \left(\bar{S}_{\bar{I}\bar{J}}\right) \right)}

   ||dev\left(\bar{\mathbf{\Sigma}}\right)|| &= \sqrt{\left(dev \left( \bar{\Sigma}_{\bar{I}\bar{J}}\right) \right)
       : \left(dev \left(\bar{\Sigma}_{\bar{I}\bar{J}}\right) \right)}

   ||dev\left(\bar{\mathbf{M}}\right)||_{\bar{K}} &= \sqrt{\left(dev \left( \bar{M}_{\bar{I}\bar{J}\left(\bar{K}\right)}\right) \right)
       : \left(dev \left(\bar{M}_{\bar{I}\bar{J}\left(\bar{K}\right)}\right) \right)}.

The notation :math:`\left(\bar{K}\right)` indicates that the :math:`\bar{K}` index is free so
:math:`||dev\left(\bar{\mathbf{M}}\right)||_{\bar{K}}` represents a
vector of invariants of the higher-order stress.
These yield functions depend on the friction angles (:math:`\phi^u`,
:math:`\phi^{\chi}`, and :math:`\phi^{\nabla\chi}`) through the functions

.. math::
   :label: micromorphic_friction_angle_functions

   &A^{u,\phi}=\beta^{u,\phi}\cos\left(\phi^u\right)
       , B^{u,\phi}=\beta^{u,\phi}\sin\left(\phi^u\right)
       , \beta^{u,\phi} = \frac{2\sqrt{6}}{3+\tilde{\beta}^{u,\phi} \sin\left(\phi^u\right)}

   &A^{{\chi},\phi}=\beta^{{\chi},\phi}\cos\left(\phi^{\chi}\right)
       , B^{{\chi},\phi}=\beta^{{\chi},\phi}\sin\left(\phi^{\chi}\right)
       , \beta^{{\chi},\phi} = \frac{2\sqrt{6}}{3+\tilde{\beta}^{\chi,\phi} \sin\left(\phi^{\chi}\right)}

   &A^{{\nabla\chi},\phi}=\beta^{{\nabla\chi},\phi}\cos\left(\phi^{\nabla\chi}\right)
       , B^{{\nabla\chi},\phi}=\beta^{{\nabla\chi},\phi}\sin\left(\phi^{\nabla\chi}\right)
       , \beta^{{\nabla\chi},\phi} = \frac{2\sqrt{6}}{3+\tilde{\beta}^{\nabla\chi,\phi} \sin\left(\phi^{\nabla\chi}\right)}.


The parameters :math:`\tilde{\beta}^{u,\phi}`, :math:`\tilde{\beta}^{\chi,\phi}`,
and :math:`\tilde{\beta}^{\nabla\chi,\phi}` are parameters ranging from
:math:`+/-` 1 to control how the yield function relates to the Mohr-Coloumb yield surface.
In a similar manner as the yield functions, the plastic potential functions for macro-
(:math:`u`), micro- (:math:`\chi`), and micro-gradient plasticity (:math:`\nabla\chi`)
are defined as

.. math::
   :label: micromorphic_plastic_potentials

   \bar{G}^u\left(\bar{\mathbf{S}},\bar{c}^u\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{S}}\right)|| - \left(A^{u,\psi}\bar{c}^u - B^{u,\psi} \bar{p}^u\right)

   \bar{G}^{\chi}\left(\bar{\mathbf{\Sigma}},\bar{c}^{\chi}\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{\Sigma}}\right)|| - \left(A^{{\chi},\psi}\bar{c}^{\chi} - B^{{\chi},\psi} \bar{p}^{\chi}\right)

   \bar{G}_{\bar{K}}^{\nabla\chi}\left(\bar{\mathbf{M}},\bar{\mathbf{c}}^{\nabla\chi}\right) &\stackrel{\text{def}}{=}
       ||dev\left(\bar{\mathbf{M}}\right)||_{\bar{K}} - \left(A^{{\nabla\chi},\psi}\bar{c}_{\bar{K}}^{\nabla\chi}
       - B^{{\nabla\chi},\psi} \bar{p}_{\bar{K}}^{\nabla\chi}\right).

The plastic potential functions depend on the dilation angles (:math:`\psi^u`,
:math:`\psi^{\chi}`, and :math:`\psi^{\nabla\chi}`) through functions of identical
form as shown in Eq. :math:numref:`{number} <micromorphic_friction_angle_functions>`
except the :math:`\phi` terms are replaced with :math:`\psi`.
Similarly, parameters :math:`\tilde{\beta}^{u,\psi}`, :math:`\tilde{\beta}^{\chi,\psi}`,
and :math:`\tilde{\beta}^{\nabla\chi,\psi}` range from :math:`+/-` 1 to control how the
potential functions relate to the Mohr-Coloumb yield surface.

The evolution equations are presented in full detail in :cite:`isbuga2017, miller_micromorphic_2021`.
For the present discussion, the evolution of the strain-like ISVs are

.. math::
   :label: micromorphic_ISV_evolution

   \dot{\bar{Z}}^u &= -\dot{\bar{\gamma}}^u \frac{\partial\bar{G}^u}{\partial\bar{c}^u} = A^{u,\psi} \dot{\bar{\gamma}}^u

   \dot{\bar{Z}}^{\chi} &= -\dot{\bar{\gamma}}^{\chi} \frac{\partial\bar{G}^{\chi}}{\partial\bar{c}^{\chi}}
      = A^{\chi,\psi} \dot{\bar{\gamma}}^{\chi}

   \dot{\bar{Z}}^{\chi}_{,\bar{I}} &= -\dot{\bar{\gamma}}^{\nabla\chi}_{\bar{J}}
        \frac{\partial\bar{G}^{\nabla\chi}_{\bar{J}}}{\partial\bar{c}^{\nabla\chi}_{\bar{I}}}
        = A^{\nabla\chi,\psi} \dot{\bar{\gamma}}^{\nabla\chi}_{\bar{J}}\delta_{\bar{I}\bar{J}}

where :math:`\dot{\bar{\gamma}}^u`, :math:`\dot{\bar{\gamma}}^{\chi}`, and :math:`\dot{\bar{\gamma}}^{\nabla\chi}_{\bar{J}}`
are plastic multipliers.
We also introduce the evolution equations for the plastic deformation which arise from the dissipation inequality as

.. math::

   \bar{H}_{\bar{I}\bar{J}}\left(\bar{\mathbf{S}}, \bar{\mathcal{Q}}\right) &\stackrel{\text{def}}{=}
       \dot{\bar{\gamma}} \frac{\partial \bar{G}}{\partial \bar{S}_{\bar{I}\bar{J}}}

   \bar{H}_{\bar{I}\bar{J}}^{\chi}\left(\bar{\mathbf{\Sigma}}, \bar{\mathcal{Q}}^{\chi}\right)
       &\stackrel{\text{def}}{=} \dot{\bar{\gamma}}^{\chi} \frac{\partial \bar{G}^{\chi}}{\partial
       \bar{\Sigma}_{\bar{I}\bar{J}}}

   \bar{H}_{\bar{I}\bar{J}\bar{K}}^{\nabla \chi}\left(\bar{\mathbf{M}}, \bar{\mathbf{\mathcal{Q}}}^{\nabla \chi}\right)
       &\stackrel{\text{def}}{=} \dot{\bar{\gamma}}^{\nabla \chi}_{\bar{L}} \frac{\partial
       \bar{G}^{\nabla\chi}_{\bar{L}}}{\partial \bar{M}_{\bar{I}\bar{J}\bar{K}}}

where

.. math::

   \bar{H}_{\bar{I}\bar{J}}\left(\bar{\mathbf{S}}, \bar{\mathcal{Q}}\right) &\stackrel{\text{def}}{=}
       \bar{L}_{\bar{K}\bar{I}}^p \bar{C}_{\bar{K}\bar{J}}^{e} -  \bar{\Psi}_{\bar{J}\bar{M}}^e
       \bar{L}_{\bar{M}\bar{N}}^{\chi,p} \chi_{\bar{N}k}^{e,-1} F_{k\bar{I}}^e

   \bar{H}_{\bar{I}\bar{J}}^{\chi}\left(\bar{\mathbf{\Sigma}}, \bar{\mathcal{Q}}^{\chi}\right)
       &\stackrel{\text{def}}{=} \bar{\Psi}_{\bar{J}\bar{M}}^e \bar{L}_{\bar{M}\bar{N}}^{\chi,p}
       \chi_{\bar{N}k}^{e,-1} F_{k\bar{I}}^e

   \bar{H}_{\bar{I}\bar{J}\bar{K}}^{\nabla \chi}\left(\bar{\mathbf{M}}, \bar{\mathbf{\mathcal{Q}}}^{\nabla
       \chi}\right) &\stackrel{\text{def}}{=} \left( \left( \chi_{j\bar{N},\bar{I}}^e \dot{\chi}_{\bar{N}I}^p
       + \chi_{j\bar{N}}^e \dot{\chi}_{\bar{N}I,\bar{I}}^p - \chi_{j\bar{N}}^e \bar{L}_{\bar{N}\bar{M}}^{\chi}
       \chi_{\bar{M}I,\bar{I}}^p \right) \chi_{I\bar{K}}^{p,-1} - \nu_{jm}^p \chi_{m\bar{K},\bar{I}}^{e}
       \right) F_{j\bar{J}}^e.

In this expression, :math:`\bar{\bf{L}}^p` is the macro-plastic velocity gradient,
:math:`\bar{\bf{L}}^{\chi,p}` is the micro-plastic velocity gradient,
:math:`\bar{\Psi}^e_{\bar{I}\bar{J}} = F_{i\bar{I}}^e \chi_{i\bar{J}}^e`,
and :math:`\nu_{ij}^p = \chi_{i\bar{I}}^e \dot{\chi}_{\bar{I}I}^p \chi_{I\bar{J}}^{p,-1} \chi_{\bar{J}j}^e`.
Finally, the cohesion terms will evolve as

.. math::
   :label: micromorphic_cohesion_evolution

   \bar{c}^u &= \bar{H}^u\bar{Z}^u

   \bar{c}^{\chi} &= \bar{H}^{\chi}\bar{Z}^{\chi}

   \bar{c}^{\nabla\chi}_{\bar{I}} &=\bar{H}^{\nabla\chi}\bar{Z}^{\chi}_{,\bar{I}}

and initialize as

.. math::
   :label: micromorphic_cohesion_initialization

   \bar{c}^u &= \bar{c}^{u,0} \text{ for } \bar{Z}^u = 0

   \bar{c}^{\chi} &= \bar{c}^{\chi, 0} \text{ for } \bar{Z}^{\chi} = 0

   \bar{c}^{\nabla\chi}_{\bar{I}} &=\bar{c}^{\nabla\chi,0} \overrightarrow{1}
       \text{ for } \bar{Z}^{\chi}_{,\bar{I}} = 0

where :math:`\overrightarrow{1}` is a vector of ones.
For this form of micromorphic elastoplasticity, a user must specify 18 parameters in addition to
the 18 linear elasticity parameters:

* initial cohesion value: :math:`c^{u,0}`, :math:`c^{\chi,0}`, and :math:`c^{\nabla\chi,0}`
* hardening (or softening if negative) moduli: :math:`\bar{H}^u`, :math:`\bar{H}^{\chi}`, and :math:`\bar{H}^{\nabla\chi}`
* friction angles: :math:`\phi^u`, :math:`\phi^{\chi}`, and :math:`\phi^{\nabla\chi}`
* yield :math:`\tilde{\beta}` parameters: :math:`\tilde{\beta}^{u,\phi}`, :math:`\tilde{\beta}^{\chi,\phi}`, and :math:`\tilde{\beta}^{\nabla\chi,\phi}`
* dilation angles: (:math:`\psi^u`, :math:`\psi^{\chi}`, and :math:`\psi^{\nabla\chi}`)
* flow potential :math:`\tilde{\beta}` parameters: :math:`\tilde{\beta}^{u,\psi}`, :math:`\tilde{\beta}^{\chi,\psi}`, and :math:`\tilde{\beta}^{\nabla\chi,\psi}`

Choosing separate definition of the yield and flow potential parameters allows for
non-associative plasticity to be modeled.

This model may be simplified by setting all friction and dilation angles to zero
resulting associative, pressure insensitive, hardening plasticity with yield functions in the form of

.. math::
   :label: micromorphic_yield_functions_simplified

   \bar{F}^u\left(\bar{\mathbf{S}},\bar{c}^u\right) = ||dev\left(\bar{\mathbf{S}}\right)||
       - 2 \sqrt{\frac{2}{3}}\bar{c}^u - \leq 0

   \bar{F}^{\chi}\left(\bar{\mathbf{\Sigma}},\bar{c}^{\chi}\right) = ||dev\left(\bar{\mathbf{\Sigma}}\right)||
       - 2 \sqrt{\frac{2}{3}}\bar{c}^{\chi} \leq 0

   \bar{F}_{\bar{K}}^{\nabla\chi}\left(\bar{\mathbf{M}},\bar{\mathbf{c}}^{\nabla\chi}\right)=
       ||dev\left(\bar{\mathbf{M}}\right)||_{\bar{K}} - 2 \sqrt{\frac{2}{3}}\bar{c}_{\bar{K}}^{\nabla\chi} \leq 0.

