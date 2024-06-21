.. _micromorphic_theory_balance:

#################
Balance Equations
#################

..
   TODO: write this section: alp = ^{\left(\alpha\right)}

..
   TODO: check that punctuation makes sense across sentences with equations

For the micromorphic formulation considered herein, we insist that the Cauchy continuum balance
equations are satisfied on the micro-scale and that macro-scale balance equations arise from this constraint.

***************
Balance of Mass
***************

The mass of a differential element :math:`dv` can be calculated as the volume average of micro-element densities $\rho^{\left(\alpha\right)}$.

.. math::

	m \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}

The conservation of mass for a non-reacting material with no additional mass sources becomes

.. math::

	\frac{Dm}{Dt} &= \frac{D}{Dt} \int_{dv} \rho^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}} = 0

	&=  \int_{dV} \frac{D}{Dt} \rho_{0}^{\left(\alpha\right)} \,{dV}^{\left(\alpha\right)} = 0

After localizing the integral, the following is true in the reference configuration

.. math::

   \frac{\partial \rho_{0}^{\left(\alpha\right)}}{\partial t} = 0

Similarly, we may define the balance of mass in the current configuration

.. math::

   \frac{Dm}{Dt} &= \int_{dV} \frac{D}{Dt} \left( \rho^{\left(\alpha\right)} J^{\left(\alpha\right)}\right) \,{dV^{\left(\alpha\right)}} = 0

   &=  \int_{dV} \left( \frac{\partial \rho^{\left(\alpha\right)}}{\partial t} J^{\left(\alpha\right)}
      + \rho^{\left(\alpha\right)} \frac{\partial J^{\left(\alpha\right)}}{\partial t}\right) \,{dV^{\left(\alpha\right)}} = 0

Recall from equation 2.177 of Holzapfel :cite:`holzapfel_nonlinear_2000`

.. math::

	\frac{\partial J}{\partial t} &= J \frac{\partial v_{i}}{\partial x_{i}}

	\Rightarrow \frac{\partial J^{\left(\alpha\right)}}{\partial t} &= J
       \frac{\partial v_{i}^{\left(\alpha\right)}}{\partial x_{i}^{\left(\alpha\right)}}

So

.. math::

	\frac{Dm}{Dt} &=  \int_{dV} \left( \frac{\partial \rho^{\left(\alpha\right)}}{\partial t}
       + \rho^{\left(\alpha\right)} \frac{\partial v_{i}^{\left(\alpha\right)}}{\partial x_{i}^{\left(\alpha\right)}}\right)
       J^{\left(\alpha\right)} \,{dV^{\left(\alpha\right)}} = 0

	&=  \int_{dv} \left( \frac{\partial \rho^{\left(\alpha\right)}}{\partial t}
       + \rho^{\left(\alpha\right)} \frac{\partial v_{i}^{\left(\alpha\right)}}{\partial x_{i}^{\left(\alpha\right)}}\right)
       \,{dV^{\left(\alpha\right)}} = 0

Localizing

.. math::

	\frac{\partial \rho^{\left(\alpha\right)}}{\partial t} + \rho^{\left(\alpha\right)}
       \frac{\partial v_{i}^{\left(\alpha\right)}}{\partial x_{i}^{\left(\alpha\right)}} = 0


If we define the average density over the differential element as

.. math::

	\rho dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}

then the total mass of the body is

.. math::

   M = \int_{B} \rho \,{dv} = \int_{B} \left[ \int_{dv} \rho^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}\right] \,{dv}
      = \int_{B^0} \left[ \int_{dV} \rho^{\left(\alpha\right)} J^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}\right] J \,{dV}

The rate of change of the total mass

.. math::

   \frac{DM}{Dt} &= \int_{B^0} \left[ \int_{dV} \frac{D \left(\rho^{\left(\alpha\right)} J^{\left(\alpha\right)}\right)}{Dt}
      \,{dv^{\left(\alpha\right)}}\right] J \,{dV}

   &= \int_{B^0} \left[ \frac{Dm}{Dt} = 0 \right] J \,{dV}

   &= 0

Which means

.. math::

   \frac{DM}{Dt} &= 0 = \frac{D}{Dt} \int_B \rho \,{dv} = \int_{B^0} \frac{D \left(\rho J \right)}{Dt} \,{dV}

   &= \int_B \left[\frac{\partial \rho}{\partial t} + \rho \frac{\partial v_i}{\partial x_i}\right] \,{dv}


Finally, upon localizing the integral, the classical balance of mass in the current configuration is recovered.

.. math::

   \frac{\partial \rho}{\partial t} + \rho \frac{\partial v_i}{\partial x_i} = 0

*************************
Balance of Micro Inertias
*************************

..
   TODO: Find that one equation for relating micro-spin inertia to moment of inertia!

Because :math`\mathbf{\Xi}` is the relative location of the mass center of :math:`dV^{\left(\alpha\right)}`
with respect to :math:`dV` we can write

.. math::

	\int_{dV} \rho_{0}^{\left(\alpha\right)} \Xi_I \,{dv^{\left(\alpha\right)}} = 0_I.

We can define the micro-moment of inertia in the reference configuration as

.. math::

	\rho_0 I_{IJ}dV \stackrel{\text{def}}{=} \int_{dV} \rho_{0}^{\left(\alpha\right)} \Xi_I \Xi_J \,{dv^{\left(\alpha\right)}}

and in the current configuration the form is

.. math::

	\rho i_{ij} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \xi_i \xi_j \,{dv^{\left(\alpha\right)}}
        = \int_{dV} \rho_{0}^{\left(\alpha\right)} \chi_{iI} \Xi_I \chi_{jJ} \Xi_J \,{dv^{\left(\alpha\right)}}

Because the deformation of :math:`dv` is affine
(i.e. :math:`\mathbf{\chi} \neq \mathbf{\chi} \left(\mathbf{x^{\left(\alpha\right)}}\right)`) which means

.. math::

   \rho i_{ij} dv &= \chi_{iI} \chi_{jJ} \int_{dv} \rho^{\left(\alpha\right)} \Xi_I \Xi_J \,{dv^{\left(\alpha\right)}}
      = \chi_{iI} \chi_{jJ} \int_{dV} \rho_{0}^{\left(\alpha\right)} \Xi_I \Xi_J \,{dv^{\left(\alpha\right)}}

   &= \rho_0 \chi_{iI} \chi_{jJ} I_{IJ} dV

Now we can relate :math:`i_{ij}` to :math:`I_{IJ}` through

.. math::
   :label: micro_inertias

   i_{ij} &= \chi_{iI} \chi_{jJ} I_{IJ} \frac{\rho_0}{\rho} \frac{dV}{dv}
      = \chi_{iI} \chi_{jJ} I_{IJ} \frac{J \rho}{\rho} \frac{dV}{JdV}

   &= \chi_{iI} \chi_{jJ} I_{IJ}

Which may be written as 

.. math::

   \mathbf{i} = \mathbf{\chi} \cdot \mathbf{I} \cdot \mathbf{\chi^T}

We now compute the balance of inertia on the micro-scale via

.. math::

   \frac{D}{Dt} \int_{B^0} \rho_0 I_{IJ} \,{dV} = \int_{B^0} \rho_{0} \frac{D I_{IJ}}{Dt} \,{dV} = 0.

From the relation between the micro-moments of inertia, we can write

.. math::

   \frac{D}{Dt} \left[\chi_{iI} \chi_{jJ} I_{IJ}\right] &= \frac{D i_{ij}}{Dt}
 
   \dot{\chi_{iI}} \chi_{jJ} I_{IJ} + \chi_{iI} \dot{\chi_{jJ}} I_{IJ}
      + \chi_{iI} \chi_{jJ} \frac{D I_{IJ}}{Dt} &= \frac{D i_{ij}}{Dt}

   \chi_{iI} \chi_{jJ} \frac{D I_{IJ}}{Dt} &= \frac{D i_{ij}}{Dt} - \dot{\chi_{iI}} \chi_{jJ} I_{IJ}
      - \chi_{iI} \dot{\chi_{jJ}} I_{IJ}

   \frac{D I_{IJ}}{Dt} &= \chi_{Li}^{-1} \chi_{Jj}^{-1} \left[\frac{D i_{ij}}{Dt}
   - \dot{\chi_{iI}} \chi_{jJ} I_{IJ} - \chi_{iI} \dot{\chi_{jJ}} I_{IJ}\right]

By substituting :math:`i_{ij}\chi_{iI}^{-1} \chi_{jJ}^{-1}` for :math:`I_{IJ}`, it can be shown that 

.. math::

	\frac{D I_{IJ}}{Dt} = \chi_{Li}^{-1} \chi_{Jj}^{-1} \left[\frac{D i_{ij}}{Dt}
       - \dot{\chi_{iI}} \chi_{Ik}^{-1} i_{kj} - \dot{\chi_{jJ}} \chi_{Jm}^{-1} i_{im}\right] 

Noting that

.. math::

   \dot{\chi_i} &= \dot{\chi_{iI}} \Xi_I = \dot{\chi_iI} \chi_{Ij}^{-1} \xi_j

   \Rightarrow \dot{\chi_{iI}} &= v_{ij} \chi_{jI}

where we have defined the micro-gyration tensor as 

.. math::

   v_{ij} \stackrel{\text{def}}{=} \dot{\chi_{iI}} \chi_{Ij}^{-1}

we can now write

.. math::

   \frac{D I_{IJ}}{Dt} = \chi_{Li}^{-1} \chi_{Jj}^{-1} \left[\frac{D i_{ij}}{Dt} - v_{ik} i_{kj} - v_{jm} i_{im}\right]

This means,

.. math::

	\int_{B^0} \rho_0 \frac{D I_{IJ}}{Dt} \,{dV} = \int_{B^0} \rho_0 \chi_{Li}^{-1} \chi_{Jj}^{-1} \left[\frac{D i_{ij}}{Dt}
       - v_{ik} i_{kj} - v_{jm} i_{im}\right] \,{dV} = 0

Or in the current configuration

.. math::

	\int_{B} \rho \frac{D I_{IJ}}{Dt} \,{dv} = \int_{B} \rho \chi_{Li}^{-1} \chi_{Jj}^{-1} \left[\frac{D i_{ij}}{Dt}
       - v_{ik} i_{kj} - v_{jm} i_{im}\right] \,{dv} = 0

If we localize the integral and, assuming for admissible deformations and materials :math:`\rho \chi_{Li}^{-1} \chi_{Jj}^{-1}`
is non-zero, we can state the balance of micro-inertia as 

.. math::

	\frac{D i_{ij}}{dt} - v_{ik} i_{kj} - v_{jk} i_{im} = 0

******************
Balance of Momenta
******************

We follow the weighted residual approach of Eringen and Suhubi to compute the balance of linear, angular, and first moment
of momentum. We postulate that the balance of linear and angular momentum are satisfied in the micro-element.

.. math::

	\sigma_{lk,l}^{\left(\alpha\right)} + \rho^{\left(\alpha\right)} \left(f_{k}^{\left(\alpha\right)} - a_{k}^{\left(\alpha\right)}\right) &= 0

	\sigma_{lk}^{\left(\alpha\right)} &= \sigma_{kl}^{\left(\alpha\right)}

where :math:`\mathbf{\sigma^{\left(\alpha\right)}}` is the micro-scale Cauchy stress, :math:`\mathbf{f^{\left(\alpha\right)}}`
is the micro-scale body force per unit mass, and :math:`\mathbf{a^{\left(\alpha\right)}}` is the micro-scale acceleration.
Following a similar approach as was done for the balance of mass we find

.. math::

   \int_B \int_{dv} \phi^{\left(\alpha\right)} \left[\sigma_{lk,l^{\left(\alpha\right)}}
      + \rho^{\left(\alpha\right)} \left(f_{k}^{\left(\alpha\right)}
      - a_{k}^{\left(\alpha\right)}\right) \right] \,{dv^{\left(\alpha\right)}} \,{dv} &= 0

   \Rightarrow \int_B \int_{dv} \left[\phi^{\left(\alpha\right)} \sigma_{lk,l}^{\left(\alpha\right)}
      + \phi^{\left(\alpha\right)} \rho^{\left(\alpha\right)} \left(f_{k}^{\left(\alpha\right)}
      - a_{k}^{\left(\alpha\right)}\right) \right] \,{dv^{\left(\alpha\right)}} \,{dv} &= 0

Where :math:`\mathbf{\phi^{\left(\alpha\right)}}` is some weighting function we will change to explore different
momentum balance conditions. Using the chain rule we find

.. math::

   \phi^{\left(\alpha\right)} \sigma_{ij,i}^{\left(\alpha\right)} = \left(\phi^{\left(\alpha\right)}
      \sigma_{ij}^{\left(\alpha\right)}\right)_{,i} - phi_{,i}^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}

Upon substitution into the integral equations we find

.. math::

   \int_B \left\{ \int_{dv} \left[\left(\phi^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}\right)_{,i}
      - \phi_{,i}^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)} + \phi^{\left(\alpha\right)}
      \rho^{\left(\alpha\right)} \left(f_{j}^{\left(\alpha\right)} - a_{j}^{\left(\alpha\right)}\right) \right]
      \,{dv^{\left(\alpha\right)}} \right\} \,{dv} &= 0

   \int_{\partial B} \left\{\int_{da} \phi^{\left(\alpha\right)}\sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)}
      \,{da}^{\left(\alpha\right)}\right\} + \int_B \left\{ \int_{dv} \left[ - \phi_{,i}^{\left(\alpha\right)}
      \sigma_{ij}^{\left(\alpha\right)} + \phi^{\left(\alpha\right)} \rho^{\left(\alpha\right)} \left(f_{j}^{\left(\alpha\right)}
      - a_{j}^{\left(\alpha\right)}\right) \right] \,{dv^{\left(\alpha\right)}} \right\} \,{dv} &= 0

==========================
Balance of Linear Momentum
==========================

We obtain the macro-balance of linear momentum by letting $\phi^{\left(\alpha\right)} = 1$ which yields

.. math::

   \int_{\partial B} \left\{\int_{da} \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)}
      \,{da^{\left(\alpha\right)}}\right\} + \int_B \left\{ \int_{dv} \left[ \rho^{\left(\alpha\right)}
      \left(f_{j}^{\left(\alpha\right)} - a_{j}^{\left(\alpha\right)}\right) \right] \,{dv^{\left(\alpha\right)}} \right\} \,{dv} = 0

We now make the following definitions

.. math::
   :label: macro_cuachy

   \sigma_{ij} n_i da \stackrel{\text{def}}{=} \int_{da} \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)} }

.. math::
   :label: macro_force

   \rho f_{j} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} f_{j}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)} }

.. math::
   :label: macro_accel

   \rho a_{j} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} a_{j}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)} }

which are the definitions for the macro-scale Cauchy stress :math:`\mathbf{\sigma}`,
the macro-scale body force :math:`\mathbf{f}`, and the macro-scale acceleration :math:`\mathbf{a}`.
This means we can write

.. math::

   \int_{\partial B} \sigma_{ij} n_i \,{da} + \int_B \rho \left(f_i - a_i \right) \,{dv} &= 0

   \int_B \left[\sigma_{ij,i} + \rho \left(f_i - a_i \right) \right] \,{dv} &= 0

Localizing the integral leads to the familiar form.

.. math::
   :label: balance_of_linear_momentum

   \sigma_{ij,i} + \rho \left(f_i - a_i \right) = 0

===========================
Balance of Angular Momentum
===========================

Letting :math:`\phi^{\left(\alpha\right)} = \varepsilon_{mkj} x_{k}^{\left(\alpha\right)}`
and recalling that :math:`x_{k}^{\left(\alpha\right)} = x_k + \xi_k`
we can solve for the balance of angular momentum and find

.. math::

   \int_{\partial B} \left\{ \int_{da}\varepsilon_{mkj}  x_{k}^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}
      n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\} &+ \int_{B} \left\{ \int_{dv} \left[-\varepsilon_{mkj}
      \delta_{ki} \sigma_{ij}^{\left(\alpha\right)} + \varepsilon_{mkj} x_{k}^{\left(\alpha\right)} \rho^{\left(\alpha\right)}
      \left(f_{j}^{\left(\alpha\right)} - a_{j}^{\left(\alpha\right)}\right)\right] \,{dv^{\left(\alpha\right)}}\right\} = 0

   \int_{\partial B} \left\{ \int_{da} \varepsilon_{mkj}  x_{k}^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}
      n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\} &+ \int_{B} \left\{ \int_{dv} \left[-\varepsilon_{mkj}
      \sigma_{kj}^{\left(\alpha\right)} + \varepsilon_{mkj} x_{k}^{\left(\alpha\right)} \rho^{\left(\alpha\right)}
      \left(f_{j}^{\left(\alpha\right)} - a_{j}^{\left(\alpha\right)}\right)\right] \,{dv^{\left(\alpha\right)}}\right\} = 0

We note that

.. math::

   a_{j}^{\left(\alpha\right)} &= \frac{D^2}{Dt^2} \left(x_j + \xi_j\right)

   &= \ddot{x_j} + \ddot{\xi_j}

   \ddot{\xi_j} = \frac{D}{Dt} \dot{\xi_j} &= \frac{D}{Dt} \left(v_{jk} \xi_k\right)

   &= \dot{v_{jk}} \xi_k + v_{jk} \dot{\xi_k}

   &= \dot{v_{jk}} \xi_k + v_{jk} v_{kn} \xi_n

   &= \left(\dot{v_{jk}} + v_{jn} v_{nk} \right) \xi_k

We now investigate the first term in the balance of angular momentum which becomes

.. math::

   \int_{\partial B} \left\{ \int_{da} \varepsilon_{mkj} x_{k}^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}
      n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}} \right\} &= \int_{\partial B} \left\{ \int_{da} \varepsilon_{mkj}
      \left(x_k + \xi_k \right) \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)}\,{da^{\left(\alpha\right)}} \right\}

   &= \int_{\partial B} \varepsilon_{mkj} x_k \left\{ \int_{da} \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)}
      \,{da^{\left(\alpha\right)}} \right\} + \int_{\partial B} \varepsilon_{mkj} \left\{ \int_{da} \xi_k
      \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}} \right\}

We now make the definition

.. math::
   :label: high_order_def

   m_{ijk} n_i da \stackrel{\text{def}}{=} \int_{da} \sigma_{ij}^{\left(\alpha\right)} \xi_k
      n_{i}^{\left(\alpha\right)}\,{da^{\left(\alpha\right)}}

where :math:`m_{ijk}` is a higher-order stress.
We can understand this quantity as a measure of the traction induced by the micro-scale Cauchy stress
on the surface of the differential element acting on the lever-arm of the micro-position vector :math:`\mathbf{\xi}`.
This stress is not only a measure of the induced moment but also includes stretching and shearing actions
scaled by the micro-position vector.
This definition allows us to write

.. math::

   \int_{\partial B} \left\{\int_{da} \varepsilon_{mkj} x_{k}^{\left(\alpha\right)}
      \sigma_{ij}^{\left(\alpha\right)} n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\}
      = \varepsilon_{mkj} \int_{\partial B} \left[x_k \sigma_{ij} n_i + m_{ijk} n_i \right] \,{da}

This allows us to then put the first term into the form

.. math::
   :label: angbal_first

   \int_{\partial B} \left\{\int_{da} \varepsilon_{mkj} x_{k}^{\left(\alpha\right)} sigma_{ij}^{\left(\alpha\right)}
      n_{i}^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\} = \varepsilon_{mkj} \int_{B} \left[\sigma_{kj}
      - x_k \sigma_{ij,i} + m_{ijk,i} \right] \,{dv}.

We now study the second term in the balance of angular momentum to find

.. math::

   \int_B \int_{dv} \varepsilon_{mkj} \sigma_{kj}^{\left(\alpha\right)} \,{dv{\left(\alpha\right)}} = \int_B \varepsilon_{mkj}
      \int_{dv} \sigma_{kj}{\left(\alpha\right)} \,{dv{\left(\alpha\right)}}.

We now make the definition

.. math::
   :label: sym_stress_def

   s_{ij} dv \stackrel{\text{def}}{=} \int_{dv} \sigma_{ij}^{\left(\alpha\right)} \,{dv{\left(\alpha\right)}},

where :math:`s_{ij}` is the symmetric micro stress.
We now substitute back into the balance of angular momentum to write

.. math::

   \varepsilon_{ijk} \int_B \left\{\left( \sigma_{kj} + m_{ijk,i} - s_{kj} \right) \,{dv}
      + \int_{dv} \left[\rho^{\left(\alpha\right)} \xi_k f_{j}^{\left(\alpha\right)} - \rho^{\left(\alpha\right)}
      \xi_k a_{j}^{\left(\alpha\right)} - x_k \left( \sigma_{ij,i}^{\left(\alpha\right)} + \rho^{\left(\alpha\right)}
      \left(f_j^{\left(\alpha\right)} - a_j^{\left(\alpha\right)}\right)\right)\right] \,{dv^{\left(\alpha\right)}}\right\} = 0.

Incorporating the balance of linear momentum and defining the body-force couple as

.. math::
   :label: body_force_couple

   \rho l_{jk} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \xi_k f_j^{\left(\alpha\right)}
      \,{dv^{\left(\alpha\right)}}

yields

.. math::

	\varepsilon_{ijk} \int_B \left\{ \left( \sigma_{kj} + m_{ijk,i} - s_{kj} + \rho l_{jk}\right) \,{dv} - \int_{dv} \rho^{\left(\alpha\right)} \xi_k a_j^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}} \right\} = 0.

Turning our attention to

.. math::

   \int_{dv} \rho^{\left(\alpha\right)} \xi_k a_j^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}} &=
      \int_{dv} \rho^{\left(\alpha\right)} \left[\xi_k \ddot{x_j} + \xi_k \ddot{\xi_j}\right]
      \,{dv^{\left(\alpha\right)}}

   &= \ddot{x_j} \int_{dv} \rho^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}}
      + \int_{dv} \rho^{\left(\alpha\right)} \xi_k \ddot{\xi_j} \,{dv^{\left(\alpha\right)}},

because :math:`\mathbf{\xi}` is defined relative to the centroid of :math:`dv`,

.. math::

   \ddot{x_j} \int_{dv} \rho^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}} = 0,

and defining the micro-spin inertia as

.. math::
   :label: micro_spin_inertia

   \rho \omega_{jk} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} \xi_k \ddot{\xi_j},
      \,{dv^{\left(\alpha\right)}},

the balance of angular momentum may be written as

.. math::
   :label: int_bal_ang_mom

   \varepsilon_{ijk} \int_{B} \left\{ \sigma_{kj} + m_{ijk,i} - s_{kj} + \rho \left( l_{jk} - \omega_{jk} \right) \right\} = 0.

Recognizing that :math:`\varepsilon_{ijk} A_{jk} = skew\left(A_{jk}\right)` we may write after localization

.. math::
   :label: loc_bal_ang_mom

   \sigma_{\left[ kj \right]} + m_{i\left[jk \right],i} + \rho \left( l_{\left[jk\right]} - \omega_{\left[jk\right]}  \right) = 0,

where :math:`\mathbf{A}_{\left[ ... \right] }` indicates the skew/anti-symmetric part of :math:`\mathbf{A}`.
Note that the balance of angular momentum indicates that the Cauchy stress is no longer guaranteed to be symmetric
indicating that point-wise couples are now allowable unlike in standard continuum theory.
This provides an additional three equations which must be solved.

=======================================
Balance of the First Moment of Momentum
=======================================

Because :math:`\mathbf{\chi}` has nine terms, we require an additional six equations beyond what the balance
of angular momentum can provide. We therefore will insist that the balance of the first moment of momentum
is satisfied. Note that moment in this case is used in the sense of a mathematical moment. We can do this by
setting  :math:`\phi' = x'_m` and writing

.. math::

   \int_{\partial B} \left\{\int_{da} x_m^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}
      n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\} &+ \int_B \left\{\int_{dv}
      \left[-\delta_{mi} \sigma_{ij}^{\left(\alpha\right)} + x_m^{\left(\alpha\right)} \rho^{\left(\alpha\right)}
      \left(f_j^{\left(\alpha\right)} - a_j^{\left(\alpha\right)}\right)\right] \,{dv^{\left(\alpha\right)} }\right\} = 0

   \int_{\partial B} \left\{\int_{da} x_m^{\left(\alpha\right)} \sigma_{ij}^{\left(\alpha\right)}
      n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}\right\} &+ \int_B \left\{\int_{dv}
      \left[-\sigma_{mj}^{\left(\alpha\right)} + x_m^{\left(\alpha\right)} \rho^{\left(\alpha\right)}
      \left(f_j^{\left(\alpha\right)} - a_j^{\left(\alpha\right)}\right)\right] \,{dv^{\left(\alpha\right)}} \right\} = 0

which leads us to a similar result as the balance of angular momentum except that we are no longer restricted to the skew/anti-symmetric part but rather

.. math::
   :label: int_bal_mom_mom

   \int_B \left\{ \sigma_{kj} + m_{ijk,i} - s_{kj} + \rho \left( l_{jk} - \omega_{jk}\right)\right\} = 0

After localization we find

.. math::
   :label: loc_bal_mom_mom

   \sigma_{kj} + m_{ijk,i} - s_{kj} + \rho \left( l_{jk} - \omega_{jk}\right) = 0,

which shows that the balance of first moment of momentum contains the balance of angular momentum as a special case.
Note that in the absence of higher order effects, (in order words :math:`\mathbf{m}`, :math:`\mathbf{l}`, and
:math:`\mathbf{\omega}` are all zero), :math:`\mathbf{\sigma} = \mathbf{s}`, i.e. the Cauchy stress is symmetric.
In this way the balance equations of micromorphic continuum theory reduce to the standard continuum theory in the limit
of no higher-order effects. Conceptually, this occurs when the deformation and/or loading is homogeneous across
:math:`dv` which is achieved when the length of the micro-position vector goes to zero among other cases.

*****************
Balance of Energy
*****************

Assuming that the classical balance of energy holds in the micro-element, we can write the energy of :math:`dv` as

.. math::

   \int_{dv} \rho^{\left(\alpha\right)} \dot{e}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}} = \int_{dv}
      \left(\sigma_{ij}^{\left(\alpha\right)} v_{j,i}^{\left(\alpha\right)} - q_{i,i}^{\left(\alpha\right)}
      + \rho^{\left(\alpha\right)} r^{\left(\alpha\right)}\right) \,{dv^{\left(\alpha\right)}},


where :math:`\dot{e}^{\left(\alpha\right)}` is the time rate of change of the micro-internal energy per unit mass,
:math:`v_{i,j}^{\left(\alpha\right)}` is the micro-velocity gradient, :math:`q_i^{\left(\alpha\right)}`
is the micro-heat flux, and :math:`r^{\left(\alpha\right)}` is the micro-heat source density per unit mass.
We can compute the balance of energy for the body :math:`B` via

.. math::

   \int_B \left\{ \int_{dv} \rho^{\left(\alpha\right)} \dot{e}^{\left(\alpha\right)}\,{dv^{\left(\alpha\right)}}\right\}
      = \int_B \left\{ \int_{dv} \left(\sigma_{ij}^{\left(\alpha\right)} v_{j,i}^{\left(\alpha\right)}
      - q_{i,i}^{\left(\alpha\right)} + \rho^{\left(\alpha\right)} r^{\left(\alpha\right)}\right)\right\}

The first term may be written as

.. math::

   \int_{dv} \rho^{\left(\alpha\right)} \dot{e}^{\left(\alpha\right)}\,{dv^{\left(\alpha\right)}}
      = \int_{dV} \rho_0^{\left(\alpha\right)} \dot{e}^{\left(\alpha\right)}\,{dV^{\left(\alpha\right)}}
      = \frac{D}{Dt} \int_{dV} \rho_0^{\left(\alpha\right)} e^{\left(\alpha\right)}\,{dV^{\left(\alpha\right)}},

which allows us to define

.. math::
   :label: micro_macro_energy

   \rho_0 e dV = \rho e dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} e^{\left(\alpha\right)}
      \,{dv^{\left(\alpha\right)}} = \int_{dV} \rho_0^{\left(\alpha\right)} e^{\left(\alpha\right)}
      \,{dV^{\left(\alpha\right)}}.

where :math:`e` is the macro-energy per unit mass.
Using :math:`x_i^{\left(\alpha\right)} = x_i + \xi_i, v_i^{\left(\alpha\right)} = v_i + \dot{\xi_i} = v_i + v_{ij}\xi_k`,
the second term becomes

.. math::

   \int_{dv} \sigma_{ij}^{\left(\alpha\right)} v_{j,i}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      &= \int_{dv} \left[\left(\sigma_{ij}^{\left(\alpha\right)} v_j^{\left(\alpha\right)} \right)_{,i}
      - \sigma_{ij,i}^{\left(\alpha\right)} v_i^{\left(\alpha\right)} \right] \,{dv^{\left(\alpha\right)}}

   &= \int_{da} \sigma_{ij}^{\left(\alpha\right)} v_j^{\left(\alpha\right)} n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}
      - \int_{dv} \sigma_{ij,i}^{\left(\alpha\right)} v_j^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}

   &= \int_{da} \sigma_{ij}^{\left(\alpha\right)} \left(v_j + v_{jk}\xi_k\right) n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}
      - \int_{dv} \sigma_{ij,i}^{\left(\alpha\right)} \left(v_j + v_{kj}\xi_k\right) \,{dv^{\left(\alpha\right)}}

   &= v_j \int_{da} \sigma_{ij}^{\left(\alpha\right)} n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}
      + v_{jk} \int_{dv} \sigma_{ij}^{\left(\alpha\right)} \xi_k n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}
      - v_j \int_{dv} \sigma_{ij,i}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      - v_{jk} \int_{dv} \sigma_{ij,i}^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}}.

Using the definitions from before and the micro-element balance of linear momentum

.. math::

   \int_{dv} \sigma_{ij}^{\left(\alpha\right)} v_{j,i}^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      &= v_j\sigma_{ij}n_{i}da + v_{jk}m_{ijk}n_{i}da

   &+ v_j \int_{dv} \rho^{\left(\alpha\right)} f_j^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      - v_j \int_{dv} \rho^{\left(\alpha\right)} a_j^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      + v_{jk} \int_{dv} \rho^{\left(\alpha\right)} f_j^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}}
      - v_{jk} \int_{dv} \rho^{\left(\alpha\right)} a_j^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}}

   &= \left(v_j\sigma_{ij} + v_{jk}m_{ijk}\right) n_{i}da + v_j\rho\left(f_j - a_j\right)dv + v_{jk}\rho l_{jk}dv
      - v_{jk} \int_{dv} \rho^{\left(\alpha\right)} \left(a_j + \ddot{\xi_j}\right)\xi_k dv^{\left(\alpha\right)}

   &= \left(v_j\sigma_{ij} + v_{jk}m_{ijk}\right) n_{i}da + v_j\rho\left(f_j - a_j\right)dv + v_{jk}\rho \left(l_{jk}
      - \omega_{jk}\right)dv - v_{jk}a_j \int_{dv} \rho^{\left(\alpha\right)} \xi_k dv^{\left(\alpha\right)}.

Recall that

.. math::

   \int_{dv} \rho^{\left(\alpha\right)} \xi_k \,{dv^{\left(\alpha\right)}} = 0.

We also make the definitions

.. math::
   :label: micro_macro_flux

   q_{i}n_{i}da \stackrel{\text{def}}{=} \int_{da} q_i^{\left(\alpha\right)} n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}

.. math::
   :label: micro_macro_heat_source

   \rho rdv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)} r^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}},

where :math:`\mathbf{q}` is the macro-scale heat flux and :math:`r` is the macro-heat source per unit mass.
The balance of energy becomes

.. math::

   \int_B \rho \dot{e} \,{dv} = \int_{\partial B} \left\{v_j\sigma_{ij}n_i + v_{jk}m_{ijk}n_i\right\} \,{da}
      - \int_{\partial B} q_{i}n_{i} \,{da} + \int_B \rho r \,{dv} + \int_B \left\{v_{j}\rho \left(f_j
      - a_j\right) + v_{jk}\rho\left(l_{jk} - \omega_{jk}\right)\right\} \,{dv},

which can also be written as

.. math::

   \int_B \rho\dot{e} \,{dv} &= \int_B \left\{\left(v_{j}\sigma_{ij}\right)_{,i} + \left(v_{jk}m_{ijk}\right)_{,i}
      + v_j\rho\left(f_j - a_j\right) + v_{jk}\rho\left(l_{jk} - \omega_{jk}\right)\right\} \,{dv}
      - \int_{\partial B} q_{i}n_{i}\,{da} + \int_B \rho r \,{dv}

   &= \int_B \left\{v_{j,i}\sigma_{ij} + v_{j}\sigma_{ij,i} + v_{jk,i}m_{ijk} + v_{jk}m_{ijk,i} + v_{j}\rho \left(f_j
      - a_j\right) + v_{jk} \rho \left(l_{jk} - \omega_{jk}\right)\right\} \,{dv}

   &- \int_{\partial B} q_{i}n_{i} \,{da} + \int_B \rho r \,{dv}.

Collecting terms,

.. math::

   \int_B \rho\dot{e} \,{dv} = \int_B \left\{v_{j,i}\sigma_{ij} + v_{jk,i}m_{ijk} + v_j\left[\sigma_{ij,i}
      + \rho\left(f_j - a_j\right)\right] + v_{jk} \left[m_{ijk,i} + \rho\left(l_{jk}
      - \omega_{jk}\right)\right]\right\} \,{dv}.

By substituting in the balance of linear momentum, and using the balance of the first moment of momentum to replace
:math:`m_{ijk,i} + \rho\left(l_{jk} - \omega_{jk}\right)` with :math:`s_{kj} - \sigma_{kj}`, we can write

.. math::
   :label: int_bal_eng

   \int_B \rho\dot{e} \,{dv} = \int_B \left\{ v_{j,i}\sigma_{ij} + v_{jk,i}m_{ijk} + v_{jk} \left[s_{kj}
      - \sigma_{kj}\right] - q_{i,i} + \rho r \right\} \,{dv}.

Upon localization,

.. math::
   :label: loc_bal_eng

   \rho\dot{e} = v_{j,i}\sigma_{ij} + v_{jk,i}m_{ijk} + v_{jk} \left[s_{kj} - \sigma_{kj}\right] - q_{i,i} + \rho r.

In the absence of higher-order effects, :math:`\mathbf{s} = \mathbf{\sigma}` and :math:`\mathbf{m} = 0`,
which means the standard continuum energy balance is recovered.

****************************
Second Law of Thermodynamics
****************************

The second law is assumed to hold in the micro-element as it does in the typical continuum such that

.. math::

   \frac{D}{Dt} \int_{dv} \rho^{\left(\alpha\right)} \eta^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}
      + \int_{da} \frac{1}{\theta^{\left(\alpha\right)}} q_i^{\left(\alpha\right)} n_i^{\left(\alpha\right)}
      \,{da^{\left(\alpha\right)}} - \int_{dv} \frac{\rho^{\left(\alpha\right)}
      r^{\left(\alpha\right)}}{\theta^{\left(\alpha\right)}} \,{dv^{\left(\alpha\right)}} \geq 0,

where :math:`\eta^{\left(\alpha\right)}` is the micro-entropy per unit mass and
:math:`\theta^{\left(\alpha\right)}` is the micro-temperature.
We now define

.. math::
   :label: micro_macro_temp

   \theta dv \stackrel{\text{def}}{=} \int_{dv} \theta^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}}

.. math::
   :label: micro_macro_entropy

   \rho \dot{\eta} dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)}
      \frac{D\\eta^{\left(\alpha\right)}}{Dt} \,{dv^{\left(\alpha\right)}}

.. math::
   :label: micro_macro_flux2

   \frac{1}{\theta}q_{i}n_{i} da \stackrel{\text{def}}{=} \int_{da} \frac{1}{\theta^{\left(\alpha\right)}}
      q_i^{\left(\alpha\right)} n_i^{\left(\alpha\right)} \,{da^{\left(\alpha\right)}}

.. math::
   :label: micro_macro_heatgen

   \frac{\rho r}{\theta} dv \stackrel{\text{def}}{=} \int_{dv} \frac{\rho^{\left(\alpha\right)} r^{\left(\alpha\right)}}
      {\theta^{\left(\alpha\right)}} \,{dv^{\left(\alpha\right)}},

where :math:`\theta` is the macro-scale temperature and :math:`\dot{\eta}` is the macro total time rate of change of
the macro-scale entropy per unit mass.
Integrate over the body

.. math::

   \int_B \left\{\rho\dot{\eta} + \left(\frac{1}{\theta} q_i\right)_{,i} - \frac{\rho r}{\theta} \right\} \,{dv} \geq 0

and expand to find

.. math::
   :label: int_2nd_law

   \int_B \left\{\rho\dot{\eta} + \frac{1}{\theta}q_{i,i} - \frac{1}{\theta^{2}} q_{i}\theta_{,i}
      - \frac{\rho r}{\theta}\right\} \,{dv} \geq 0.

Localizing the integral yields

.. math::
   :label: loc_2nd_law

   \rho\theta\dot{\eta} + q_{i,i} - \frac{1}{\theta} q_{i}\theta_{,i} - \rho r \geq 0.

Note that there is no significant variation of the second law from the classical form and
the introduction of a micro temperature has no effect on the homogenized response.

*************************
Clausius-Duhem Inequality
*************************

A macro-scale Helmholtz free energy per unit mass may introduced as

.. math::
   :label: helm

   \psi = e - \theta\eta.

Alternatively, the free energy could be expressed as

.. math::
   :label: micro_macro_helm

   \psi dv \stackrel{\text{def}}{=} \int_{dv} \psi^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}},

or, possibly more appropriately

.. math::
   :label: micro_macro_helmholtz

   \rho \psi dv \stackrel{\text{def}}{=} \int_{dv} \rho^{\left(\alpha\right)}\psi^{\left(\alpha\right)} \,{dv^{\left(\alpha\right)}},

since it is consistent with the definitions of :math:`e` and :math:`\eta`.
These definitions have the advantage that the macroscale Helmholtz free energy is directly tied to the microscale.
This definition provides an explicit connection between the behavior of a macroscale
micromorphic material and the micro-scale material it is \underline{derived} from.
In this way, a micromorphic interpretation of a classical constitutive model is possible as a consequence of defining
the fundamental length scale.

This allows us to write

.. math::

   \dot{\psi} &= \dot{e} - \dot{\theta}\eta - \theta\dot{\eta}

   \Rightarrow \rho\theta\dot{\eta} &= \rho\dot{e} - \rho\dot{\theta}\eta - \rho\dot{\psi},

which means the second law can be rewritten as

.. math::

   \rho\dot{e} - \rho\dot{\theta}\eta - \rho\dot{\psi} + q_{i,i} - \frac{1}{\theta} q_{i}\theta_i - \rho r \geq 0.

The first law of thermodynamics may be substituted into the second law and simplified as

.. math::
   :label: CDI

   v_{j,i}\sigma_{ij} + v_{jk,i}m_{ijk} + v_{jk}\left[s_{kj} - \sigma_{kj}\right] - \rho \left(\dot{\theta}\eta
   + \dot{\psi}\right) - \frac{1}{\theta}q_{i}\theta_{,i} \geq 0.

This form will reduce to the classical form so long as higher-order terms are absent.