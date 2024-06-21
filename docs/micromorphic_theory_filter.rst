.. _micromorphic_theory_filter:

*******************
Micromorphic Filter
*******************

.. note::

   The purpose of this section is to describe the theory and implementation of the Micromorphic
   Filter, however, this effort is a work-in-progress and otherwise incomplete.

As described in sections :ref:`workflow_theory_filter` and :ref:`workflow_homogenization`,
the Micromorphic Filter calculates stress and deformation measures (see equations
:math:numref:`{number} <homogenized_quantities>` and :math:numref:`{number} <deformation_measures>`)
from DNS data over one or more macroscopic filter domains (:math:`\mathcal{G}`)
which contain one or more micro-averaging domains (:math:`\mathcal{G}^\beta`),
where :math:`\beta` is the index corresponding to the number and ID of the micro-averaging
domains.
The following figure depicts a single filter domain with a single micro-averaging domain
highlighted for a collection of micro- (or DNS) points, where :math:`\alpha` is an index
corresponding to the number and ID of the micro-points.

.. figure:: averaging_domains.png
   :name: averaging-domains
   :align: center
   :width: 35%

   Micromorphic Filter :math:`\mathcal{G}` showing the macroscale volume associated with the :math:`\beta\text{th}` Gauss domain (:math:`\mathcal{G}^\beta`) and micro-scale point :math:`\alpha`.

The Micromorphic Filter first determines the macro-scale kinematic quantities for displacements
:math:`\mathbf{u}` and micro-deformation :math:`\mathbf{\phi}` followed by the volume
and surface integral quantities defined in Eq. :math:numref:`{number} <homogenized_quantities>`.
The Micromorphic Filter does not actually calculate the macroscale deformation gradient
:math:`\mathbf{F}`, the micro-deformation tensor :math:`\mathbf{\chi}` or the other deformation measures
in equation :math:numref:`{number} <deformation_measures>`, however, it does provide the terms
necessary to calculate these quantities including
:math:`u_i`, :math:`\frac{\partial u_i}{\partial X_I}`, :math:`\phi_{iI}`, and
:math:`\frac{\partial \phi_{iI}}{\partial X_J}`.
Refer to :py:mod:`model_package.xdmf_reader_tools.compute_deformations` to see how quantities
in equation :math:numref:`{number} <deformation_measures>` are calculated from the
Micromorphic Filter output.

It should be noted that the micro-deformation tensor may be related to a micro-displacement tensor
:math:`\mathbf{\phi}` as :math:`\mathbf{\chi} = \mathbf{1} + \mathbf{\phi}`.
However, as pointed out by Isbuga 2011 :cite:`Isbuga2011`, this micro-displacement tensor
only exists as a numerical means to interpolate and solve for nine additional
dofs that are used to calculate :math:`\mathbf{\chi}`. Although this is only a small detail
within the context of this theory, it is an important for user's of the Micromorphic Filter.

Note that more than 1 micro-averaging domain is needed to properly
determine the relevant deformations.