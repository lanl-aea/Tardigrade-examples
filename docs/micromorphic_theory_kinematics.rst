.. _micromorphic_theory_kinematics:

##########
Kinematics
##########

The micromorphic continuum theory developed by Eringen :cite:`eringen_nonlinear_1964, suhubl_nonlinear_1964`
accounts for behavior of a micro-structure coupled with a macro-scale.
The approach of classical nonlinear continuum mechanics is enriched by embedding
micro-position vectors corresponding to micro-scale differential elements.

Figure :numref:`{number} <micromorphic_configuration_out_no_beta>` by a body
undergoing deformation from the reference
configuration (:math:`B_0`, :math:`\partial B_0`) to the current configuration (:math:`B`, :math:`\partial B`),
with differential elements and their encompassing areas (:math:`dV`, :math:`dA`) and (:math:`dv`, :math:`da`),
respectively.
Micro-scale quantities are identified with a superscript alpha :math:`(\bullet)^{\left(\alpha\right)}`,
and the same capitalized notation distinguishes between the reference (uppercase)
and current (lowercase) configurations.

.. figure:: micromorphic_configuration_out_no_beta.svg
   :name: micromorphic_configuration_out_no_beta
   :align: center
   :width: 50%

   The reference and current configurations with detail of differential elements

The micro-position vectors in the reference Eq. (:math:`\mathbf{X}^{\left(\alpha\right)}`) and current configurations
(:math:`\mathbf{x}^{\left(\alpha\right)}`), respectively, are shown in :math:numref:`{number} <micro_positions>`.
Vectors :math:`\mathbf{X}` and :math:`\mathbf{x}` point to the centers of mass of the differential macro-elements
and :math:`\mathbf{\Xi}` and :math:`\bf{\xi}` point to the centers of mass of the micro-elements relative to the
centers of mass of the differential elements.
It is assume that the magnitudes of micro-position vectors are significantly small, e.g.,
:math:`||\mathbf{\Xi}|| << 1` and :math:`||\mathbf{\xi}|| << 1`.

.. math::
   :label: micro_positions

   X_i^{\left(\alpha\right)} &= X_i + \Xi_i

   x_i^{\left(\alpha\right)} &= x_i + \xi_i

********************
Deformation Gradient
********************

Equation :math:numref:`{number} <ximap>` shows the smooth map between differential micro-elements in the
reference and current configurations through the the micro-deformation tensor :math:`\mathbf{\chi}`.
While the macro-position vectors (:math:`\mathbf{X}` and :math:`\mathbf{x}`) are mapped to their respective
differential elements (:math:`d\mathbf{X}` and :math:`d\mathbf{x}`) through the deformation gradient
(:math:`\mathbf{F}=\mathbf{I}+\partial \mathbf{u} / \partial \mathbf{X}`), the micro-scale relative position
vectors (:math:`\mathbf{\Xi}` and :math:`\mathbf{\xi}`) are assumed to map in their entirety through the
linear transformation :math:`\mathbf{\chi}`.

.. math::
   :label: dxmap

   dx_i = F_{iI} dX_I

.. math::
   :label: ximap

   \xi_i = \chi_{iI} \Xi_I

The micro-deformation tensor :math:`\mathbf{\chi}` is related to the micro-displacement tensor
:math:`\mathbf{\Phi}` as :math:`\mathbf{\chi}= \mathbf{I} + \mathbf{\Phi}`.

The micro-position deformation gradient :math:`\mathbf{F}^{\left(\alpha\right)}` must now be computed to
relate the differential micro-elements positions :math:`d\mathbf{X}^{\left(\alpha\right)}` and
:math:`d\mathbf{x}^{\left(\alpha\right)}`.

.. math::
   :label: Fprime_1

   F_{iI}^{\left(\alpha\right)} &= \frac{\partial}{\partial X_I^{\left(\alpha\right)}} (x_i + \xi_i)

   &= \frac{\partial x_i}{\partial X_L} \frac{\partial X_L}{\partial X_{I}^{\left(\alpha\right)}}
      + \frac{\partial}{\partial X_{I}^{\left(\alpha\right)}} (\chi_{iK} \Xi_{K})

   &= F_{iL} \frac{\partial X_L}{\partial X_I^{\left(\alpha\right)}}
      + \frac{\partial \chi_{iK}}{\partial X_L} \frac{\partial X_L}{\partial X_I^{\left(\alpha\right)}}\Xi_{K}
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L}\frac{\partial X_L}{\partial X_I^{\left(\alpha\right)}}

   &= \left( F_{iL} + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \frac{\partial X_L}{\partial X_I^{\left(\alpha\right)}}

Since :math:`X_L^{\left(\alpha\right)} = X_L + \Xi_L`, then :math:`X_L = X_L^{\left(\alpha\right)} - \Xi_L`, so

.. math::
   :label: XdelXprime

	\frac{\partial X_L}{\partial X_I^{\left(\alpha\right)}} = \frac{\partial X_L^{\left(\alpha\right)}}{\partial X_I^{\left(\alpha\right)}}
       - \frac{\partial \Xi_L}{\partial X_I^{\left(\alpha\right)}} = \delta_{LI} - \frac{\partial \Xi_L}{\partial X_I^{\left(\alpha\right)}}

Similarly,

.. math::
   :label: XprimedelX

   \frac{\partial X_I^{\left(\alpha\right)}}{\partial X_L} = \frac{\partial X_I}{\partial X_L}
      + \frac{\partial \Xi_I}{\partial X_L} = \delta_{IL} + \frac{\partial \Xi_I}{\partial X_L}

Substitute Eq. :math:numref:`{number} <XdelXprime>` back into :math:`\mathbf{F}^{\left(\alpha\right)}`
in Eq. :math:numref:`{number} <Fprime_1>` and manipulate.

.. math::
   :label: F_prime_2

   F_{iI}^{\left(\alpha\right)} &= \left( F_{iL} + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \left( \delta_{LI}
      - \frac{\partial \Xi_L}{\partial X_I^{\left(\alpha\right)}} \right)

   &= F_{iL}\delta_{LI}  + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K \delta_{LI}
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L}\delta_{LI} - \left( F_{iL}
      + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \frac{\partial \Xi_L}{\partial X_I^{\left(\alpha\right)}}

   &= F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K + \chi_{iK}\frac{\partial \Xi_K}{\partial X_I}
      - \left( F_{iL} + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \frac{\partial \Xi_L}{\partial X_I^{\left(\alpha\right)}}

We may inspect the :math:`\frac{\partial \mathbf{\Xi}}{\partial \mathbf{X}^{\left(\alpha\right)}}` term and
substitute the expression from Eq. :math:numref:`{number} <XprimedelX>`.

.. math::
   :label: XidelXprime

   \frac{\partial \Xi_{L}}{\partial X_{I}^{\left(\alpha\right)}} &= \frac{\partial \Xi_{L}}{\partial X_{M}} \frac{\partial X_{M}}{\partial X_{I}^{\left(\alpha\right)}}
      = \frac{\partial \Xi_{L}}{\partial X_{M}} \left( \frac{\partial X_{I}^{\left(\alpha\right)}}{\partial X_{M}} \right)^{-1}

   &= \frac{\partial \Xi_{L}}{\partial X_{M}} \left( \delta_{IM}
      + \frac{\partial \Xi_I}{\partial X_M}\right)^{-1}

Substitute back into :math:`\mathbf{F}^{\left(\alpha\right)}`. Eq. :math:numref:`{number} <F_prime_general>` shows the
\textit{general} form of deformation gradient!

.. math::
   :label: F_prime_general

	\boxed{F_{iI}^{\left(\alpha\right)} = F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K
       + \chi_{iK}\frac{\partial \Xi_K}{\partial X_I} - \left( F_{iL} + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
       + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \frac{\partial \Xi_{L}}{\partial X_{M}} \left( \delta_{IM}
       + \frac{\partial \Xi_I}{\partial X_M}\right)^{-1}}

=============================================
Simplification for small variation in density
=============================================

The micro-position deformation gradient may be further simplified for different cases. One such case is when
:math:`\frac{\partial \mathbf{\Xi}}{\partial \mathbf{X}} << 1` which indicates that the distribution of
:math:`\mathbf{\Xi}` is nearly the same at every location within the body. Thus, the variation in the mass
distribution is small between different :math:`dV`, however, the variation in mass itself is not necessarily small.

If :math:`\frac{\partial \mathbf{\Xi}}{\partial \mathbf{X}} << 1`, then

.. math::

	\left( \delta_{IM} + \frac{\partial \Xi_I}{\partial X_M}\right)^{-1} \approx \left( \delta_{IM} - \frac{\partial \Xi_I}{\partial X_M}\right)


which we may substitute into Eq. :math:numref:`{number} <XidelXprime>`

.. math::
   :label: XidelXprime2

   \frac{\partial \Xi_{L}}{\partial X_{I}^{\left(\alpha\right)}} &= \frac{\partial \Xi_{L}}{\partial X_{M}} \left( \delta_{IM}
      + \frac{\partial \Xi_I}{\partial X_M}\right)^{-1} \approx \frac{\partial \Xi_{L}}{\partial X_{M}} \left( \delta_{IM}
      - \frac{\partial \Xi_I}{\partial X_M}\right)

   &\approx \frac{\partial \Xi_{L}}{\partial X_{M}} \delta_{IM}
      + \cancelto{0}{\frac{\partial \Xi_{L}}{\partial X_{M}}\frac{\partial \Xi_{I}}{\partial X_{M}}}

   &\approx \frac{\partial \Xi_{L}}{\partial X_{I}}

One may insert the results of Eq. :math:numref:`{number} <XidelXprime2>` into Eq. :math:numref:`{number} <F_prime_2>` to provide

.. math::
   :label: F_prime_3

   F_{iI}^{\left(\alpha\right)} = F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_I} - \left( F_{iL} + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K
      + \chi_{iK}\frac{\partial \Xi_K}{\partial X_L} \right) \frac{\partial \Xi_L}{\partial X_I}

This expression may be expanded and a single term may be canceled.

.. math::

	F_{iI}^{\left(\alpha\right)} = F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K + \chi_{iK}\frac{\partial \Xi_K}{\partial X_I}
       - F_{iL}\frac{\partial \Xi_L}{\partial X_I} - \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K\frac{\partial \Xi_L}{\partial X_I}
       - \cancelto{0}{\chi_{iK}\frac{\partial \Xi_K}{\partial X_L}\frac{\partial \Xi_L}{\partial X_I}}

The dummy indices on the third term may be switched from :math:`K` to :math:`L`.

.. math::

	F_{iI}^{\left(\alpha\right)} = F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K
       + \chi_{iL}\frac{\partial \Xi_L}{\partial X_I} - F_{iL}\frac{\partial \Xi_L}{\partial X_I}
       - \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K\frac{\partial \Xi_L}{\partial X_I}

Finally, recollect terms to give:

.. math::
   :label: final_F

	F_{iI}^{\left(\alpha\right)} = F_{iI} + \frac{\partial \chi_{iK}}{\partial X_I}\Xi_K - \left( F_{iL}
       + \frac{\partial \chi_{iK}}{\partial X_L}\Xi_K - \chi_{iL} \right) \frac{\partial \Xi_L}{\partial X_I}

Materials with a large variation in density can include functionally graded materials, additively manufactured lattice structures,
granular materials with particles spreading out, and some foams (with spatially varying density). Micromorphic descriptions of
these materials need to use the micro-position deformation gradient from Eq. :math:numref:`{number} <F_prime_general>`.

********************
Deformation Measures
********************

We now investigate several deformation measures.
For now, we will assume only elastic deformations and will introduce more specific kinematics
for an elastoplastic model using a multiplicative decomposition of the deformation gradient
and micro deformation tensor.
A full treatment of deformation measures will be added in the future, but for now we may start
with the following definitions.

.. math::
   :label: deformation_measures_1

   \mathcal{C}_{IJ} &= F_{iI} F_{iJ}

   \Psi_{IJ} &= F_{iI} \xi_{iJ}

   \Gamma_{IJK} &= F_{iI} \xi_{iJ,K}

These measures may be used to define the Green-Lagrange strain (Eq. :math:numref:`{number} <GL_strain>`)
and the Micro strain (Eq. :math:numref:`{number} <micro_strain>`) which will be used for an elastic
constitutive model along with the Micro-deformation gradient :math:`\Gamma_{IJK}`.

.. math::
   :label: GL_strain

   E_{IJ} = \frac{1}{2} \left( \mathcal{C}_{IJ} - \delta_{IJ} \right)

.. math::
   :label: micro_strain

   \mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}

..
   TODO: Fill out this section with more details (e.g., differential line segments)