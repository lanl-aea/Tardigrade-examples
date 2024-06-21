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

   S_{IK} = 2 \frac{\left(\rho_{0}\psi\right)^e}{\partial \mathcal{C}_{IJ}}
      + \frac{\left(\rho_{0}\psi\right)^e}{\partial \Psi_{IQ}} \Psi_{KQ} \mathcal{C}_{JK}^{-1}
      + \frac{\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{IQK}} \Gamma_{SQK} \mathcal{C}_{JS}^{-1}

.. math::
   :label: SIGMA_1

   \Sigma_{IJ} = 2 \frac{\left(\rho_{0}\psi\right)^e}{\partial \mathcal{C}_{IJ}}
      + 2 symm\left[ \frac{\left(\rho_{0}\psi\right)^e}{\partial \Psi_{IQ}} \Psi_{KQ} \mathcal{C}_{JK}^{-1}
      + \frac{\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{IQK}} \Gamma_{SQK} \mathcal{C}_{JS}^{-1} \right]

.. math::
   :label: M_1

   M_{IJK} = \frac{\left(\rho_{0}\psi\right)^e}{\partial \Gamma_{JKI}}.

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

****************
Smith Conditions
****************

Refer to the :ref:`linear_elastic_constraints` section for the discussion of the Smith
conditions.

..
   TODO: Describe the Smith conditions

******************************
Micromorphic Elasto-plasticity
******************************

..
   TODO: Write this section

..
   TODO: Add in discussion of elasto-plastic kinematics

.. note::

    Discussion coming soon!