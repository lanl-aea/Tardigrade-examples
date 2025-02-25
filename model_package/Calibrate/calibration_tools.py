#!python
import yaml

import matplotlib.pyplot
import numpy

try:
    import micromorphic
except:
    print('micromorphic library not imported')

try:
    import linear_elastic_parameter_constraint_equations as constraints
except:
    print('linear_elastic_parameter_constraint_equations not imported')


def average_quantities(quantities, type, elem):
    '''Average tensor quantites over 8 quadrature points

    :param dict quantities: A 2nd or 3rd order tensor dictionary with keys for quadrature points and values storing an array where indices correspond to time, element number, and tensor components
    :param str type: A string specifying the type of tensor to average. Use "3" for a vector. Use "3x3" for a regular second order tensor. Use "9" for a flattened second order tensor. Use "3x3x3" for a third order tensor.
    :param int elem: The macro (filter) element to calibrate

    :returns: ``output`` dict with same indices as ``quantities`` and a single key
    '''

    output = {}
    shapes = numpy.shape(quantities[0])

    if type == '9':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for k in range(9):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,elem,k])
            means = numpy.mean(mean_field, axis=0)
            output[0][:,elem,k] = means
    elif type == '3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3]))
        for i in range(3):
            for j in range(3):
                mean_field = []
                for qp in quantities.keys():
                    mean_field.append(quantities[qp][:,elem,i,j])
                means = numpy.mean(mean_field, axis=0)
                output[0][:,0,i,j] = means
    elif type == '3x3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3], shapes[4]))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    mean_field = []
                    for qp in quantities.keys():
                        mean_field.append(quantities[qp][:,elem,i,j,k])
                    means = numpy.mean(mean_field, axis=0)
                    output[0][:,0,i,j,k] = means
    elif type == '3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for i in range(3):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,elem,i])
            means = numpy.mean(mean_field, axis=0)
            output[0][:,0,i] = means

    return(output)


def isolate_element(quantities, type, elem):
    '''Isolate the homogenized quantities for a specified element

    :param dict quantities: A 2nd or 3rd order tensor dictionary with keys for quadrature points and values storing an array where indices correspond to time, element number, and tensor components
    :param str type: A string specifying the type of tensor to average. Use "3" for a vector. Use "3x3" for a regular second order tensor. Use "9" for a flattened second order tensor. Use "3x3x3" for a third order tensor.
    :param int elem: The macro (filter) element to calibrate

    :returns: ``output`` dict with same indices as ``quantities`` and a single key
    '''

    output = {}
    shapes = numpy.shape(quantities[0])

    for qp in range(8):
        if type == '9':
            output[qp] = numpy.zeros((shapes[0], 1, shapes[2]))
            for k in range(9):
                output[qp][:,0,k] = quantities[qp][:,elem,k]
        elif type == '3x3':
            output[qp] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3]))
            for i in range(3):
                for j in range(3):
                    output[qp][:,0,i,j] = quantities[qp][:,elem,i,j]
        elif type == '3x3x3':
            output[qp] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3], shapes[4]))
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        output[qp][:,0,i,j,k] = quantities[qp][:,elem,i,j,k]
        elif type == '3':
            output[qp] = numpy.zeros((shapes[0], 1, shapes[2]))
            for i in range(3):
                output[qp][:,0,i] = quantities[qp][:,elem,i]

    return(output)


def isolate_element_and_qp(quantities, type, elem, qp):
    '''solate the homogenized quantities for a specified element and quadrature point

    :param dict quantities: A 2nd or 3rd order tensor dictionary with keys for quadrature points and values storing an array where indices correspond to time, element number, and tensor components
    :param str type: A string specifying the type of tensor to average. Use "3" for a vector. Use "3x3" for a regular second order tensor. Use "9" for a flattened second order tensor. Use "3x3x3" for a third order tensor.
    :param int elem: The macro (filter) element to calibrate
    :param int qp: The quadrature point of the macro (filter) element to calibrate

    :returns: ``output`` dict with same indices as ``quantities`` and a single key
    '''

    output = {}
    shapes = numpy.shape(quantities[0])

    if type == '9':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for k in range(9):
            output[0][:,0,k] = quantities[qp][:,elem,k]
    elif type == '3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3]))
        for i in range(3):
            for j in range(3):
                output[0][:,0,i,j] = quantities[qp][:,elem,i,j]
    elif type == '3x3x3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2], shapes[3], shapes[4]))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    output[0][:,0,i,j,k] = quantities[qp][:,elem,i,j,k]
    elif type == '3':
        output[0] = numpy.zeros((shapes[0], 1, shapes[2]))
        for i in range(3):
            output[0][:,0,i] = quantities[qp][:,elem,i]

    return(output)


def Isbuga_micrormorphic_elasticity_parameters(Emod, nu, Lc, case_1_override=False):
    '''Calculate initial estimate of 18 parameter micromorphic linear elasticity model parameters using method defined in https://doi.org/10.1016/j.ijengsci.2011.04.006

    :param float Emod: An estimate of homogenized elastic modulus
    :param float nu: An estimate of the homogenized Poisson ratio
    :param float Lc: An estimate of the length scale parameter

    :returns: array of estimated micromorphic linear elasticity parameters
    '''

    # calculate "classic" lame parameters
    lame_lambda = Emod*nu/((1.+nu)*(1.-2*nu))
    lame_mu     = Emod/(2*(1.+nu)) #shear modulus, K

    # estimate micromorphic parameters
    lamb = 0.7435*lame_lambda
    mu = 0.583*lame_mu
    eta = 1.53*lame_lambda
    tau = 0.256*lame_lambda
    kappa = 0.833*lame_mu
    nu_new = 0.667*lame_mu
    sigma = 0.4167*lame_mu

    tau_1 = 0.111*(lame_lambda*Lc*Lc)
    tau_2 = 0.185*(lame_lambda*Lc*Lc)
    tau_3 = 0.185*(lame_lambda*Lc*Lc)
    tau_4 = 0.204*(lame_lambda*Lc*Lc)
    tau_5 = 0.1*(lame_lambda*Lc*Lc)
    tau_6 = 0.256*(lame_lambda*Lc*Lc)
    tau_7 = 0.670*(lame_mu*Lc*Lc)
    tau_8 = 0.495*(lame_mu*Lc*Lc)
    tau_9 = 0.495*(lame_mu*Lc*Lc)
    tau_10 = 0.408*(lame_mu*Lc*Lc)
    tau_11 = 0.495*(lame_mu*Lc*Lc)

    if case_1_override == True:
        lamb = lame_lambda
        mu = lame_mu

    # collect
    parameters = numpy.array([lamb, mu, eta, tau, kappa, nu_new, sigma,
                              tau_1, tau_2, tau_3, tau_4, tau_5, tau_6,
                              tau_7, tau_8, tau_9, tau_10, tau_11])

    return(parameters)


def plot_stresses(strain, stress, stress_sim, output_name, element, nqp, x_label_base, y_label_base,
                  increment=None, find_bounds=False):
    '''Plot comparison of stress vs strain between homogenized DNS results against calibrated model predictions

    :param dict strain: The quantities dict storing a strain measure
    :param dict stress: The quantities dict storing a homogenized DNS stress
    :param dict stress_sim: The quantities dict storing a calibrated stress
    :param str output_name: The output plot name
    :param int element: The macro (filter) element considered for calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param str x_label_base: A string to include in the plot x-label
    :param str y_label_base: A string to include in the plot y-label
    :param list increment: An optional list of one or more increments to plot restults
    :param bool find_bounds: Boolean specifying whether or not to identify common y-axis bounds for all subplots

    :returns: ``output_name``
    '''

    name = output_name.replace('.PNG','')
    fig1 = matplotlib.pyplot.figure(name, figsize=(11,10))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [0., 1.]

    if increment:
        inc = [int(i) for i in increment]
    else:
        inc = [i for i in range(0, numpy.shape(strain[0][:,0,0,0])[0])]

    colors = matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    e = 0
    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            x_label = r"$" + str(x_label_base) + "_{" + str(i+1) + str(j+1) + "}$"
            y_label = r"$" + str(y_label_base) + "_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if nqp == 1:
                ax1.plot(strain[0][inc,e,i,j], stress[0][inc,e,i,j], 'o', label='Filter')
                ax1.plot(strain[0][inc,e,i,j], stress_sim[0][inc,e,i,j], '-', label='Fit')
                ybounds[0] = numpy.min([ybounds[0], numpy.min([stress[0][inc,e,i,j], stress_sim[0][inc,e,i,j]])])
                ybounds[1] = numpy.max([ybounds[1], numpy.max([stress[0][inc,e,i,j], stress_sim[0][inc,e,i,j]])])
            else:
                for qp in range(nqp):
                    ax1.plot(strain[qp][inc,e,i,j], stress[qp][inc,e,i,j], 'o', color=colors[qp], label=f'Filter, qp #{qp+1}')
                    ax1.plot(strain[qp][inc,e,i,j], stress_sim[qp][inc,e,i,j], '-', color=colors[qp], label=f'Fit, qp #{qp+1}')
                    ybounds[0] = numpy.min([ybounds[0], numpy.min([stress[qp][inc,e,i,j], stress_sim[qp][inc,e,i,j]])])
                    ybounds[1] = numpy.max([ybounds[1], numpy.max([stress[qp][inc,e,i,j], stress_sim[qp][inc,e,i,j]])])
            ax1.set_xlabel(x_label, fontsize=14)
            ax1.set_ylabel(y_label, fontsize=14)
            ax1.set_xticks(ax1.get_xticks(), ax1.get_xticklabels(), rotation=45)
            ax1.xaxis.set_major_formatter(matplotlib.pyplot.ScalarFormatter())
            ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    handles, labels = ax1.get_legend_handles_labels()
    if nqp != 1:
        labels = ['Filter', 'Fit']
    if find_bounds == True:
        for i in range(3):
            for j in range(3):
                mult=1.02
                ybounds = [mult*ybounds[0], mult*ybounds[1]]
                axes1[i][j].set_ybound(ybounds)
    fig1.legend(handles[0:2], labels[0:2], loc='lower center', bbox_to_anchor=(0.52, 0.00), ncols=2, fontsize=14)
    fig1.tight_layout()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.subplots_adjust(bottom=0.12)
    fig1.savefig(f'{output_name}')

    return 0


def plot_higher_order_stresses(Gamma, M, M_sim, output_name, element, nqp, increment=None, find_bounds=False):
    '''Plot comparison of higher order stress and micro-deformation gradient between homogenized DNS results against calibrated model predictions

    :param dict Gamma: The quantities dict storing the micro-deformation gradient
    :param dict M: The quantities dict storing higher order stress from the homogenized DNS results
    :param dict M_sim: The quantities dict storing higher order stress from the calibrated model predictions
    :param str output_name: The output plot name
    :param int element: The macro (filter) element considered for calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param list increment: An optional list of one or more increments to plot restults
    :param bool find_bounds: Boolean specifying whether or not to identify common y-axis bounds for all subplots

    :returns: ``output_name``
    '''

    name = output_name.replace('.PNG','')
    fig1 = matplotlib.pyplot.figure(name, figsize=(12,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [-0.01, 0.01]

    if increment:
        inc = [int(i) for i in increment]
    else:
        inc = [i for i in range(0, numpy.shape(M[0][:,0,0,0,0])[0])]

    colors = matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    e = 0
    for i in range(3):
        for j in range(3):
            #for k in range(3)
            ax1 = axes1[i][j]
            plot_label = r"$M_{" + str(i+1) + str(j+1) + "K}$ (MPa)"

            if nqp == 1:
                ax1.plot(Gamma[0][inc,e,i,j,0], M[0][inc,e,i,j,0], 'o', color=colors[0], label='Filter, K=1')
                ax1.plot(Gamma[0][inc,e,i,j,0], M_sim[0][inc,e,i,j,0], '-', color=colors[0], label='Fit, K=1')
                ax1.plot(Gamma[0][inc,e,i,j,1], M[0][inc,e,i,j,1], 'o', color=colors[1], label='Filter, K=2')
                ax1.plot(Gamma[0][inc,e,i,j,1], M_sim[0][inc,e,i,j,1], '-', color=colors[1], label='Fit, K=2')
                ax1.plot(Gamma[0][inc,e,i,j,2], M[0][inc,e,i,j,2], 'o', color=colors[2], label='Filter, K=3')
                ax1.plot(Gamma[0][inc,e,i,j,2], M_sim[0][inc,e,i,j,2], '-', color=colors[2], label='Fit, K=3')
                lowers = numpy.min([M[0][inc,e,i,j,:], M_sim[0][inc,e,i,j,:]])
                uppers = numpy.max([M[0][inc,e,i,j,:], M_sim[0][inc,e,i,j,:]])
                ybounds[0] = numpy.min([ybounds[0], lowers])
                ybounds[1] = numpy.max([ybounds[1], uppers])
            else:
                for qp in range(nqp):
                    ax1.plot(Gamma[qp][inc,e,i,j,0], M[qp][inc,e,i,j,0], 'o', color=colors[qp], label=f'Filter, K=1')
                    ax1.plot(Gamma[qp][inc,e,i,j,0], M_sim[qp][inc,e,i,j,0], '-', color=colors[qp], label=f'Fit, K=1')
                    ax1.plot(Gamma[qp][inc,e,i,j,1], M[qp][inc,e,i,j,1], '^', color=colors[qp], label=f'Filter, K=2')
                    ax1.plot(Gamma[qp][inc,e,i,j,1], M_sim[qp][inc,e,i,j,1], ':', color=colors[qp], label=f'Fit, K=2')
                    ax1.plot(Gamma[qp][inc,e,i,j,2], M[qp][inc,e,i,j,2], 'v', color=colors[qp], label=f'Filter, K=3')
                    ax1.plot(Gamma[qp][inc,e,i,j,2], M_sim[qp][inc,e,i,j,2], '-.', color=colors[qp], label=f'Fit, K=3')
                    lowers = numpy.min([M[qp][inc,e,i,j,:], M_sim[qp][inc,e,i,j,:]])
                    uppers = numpy.max([M[qp][inc,e,i,j,:], M_sim[qp][inc,e,i,j,:]])
                    ybounds[0] = numpy.min([ybounds[0], lowers])
                    ybounds[1] = numpy.max([ybounds[1], uppers])
            ax1.set_xlabel(r"$\Gamma_{" + str(i+1) + str(j+1) + "K}$", fontsize=14)
            ax1.set_ylabel(plot_label, fontsize=14)
            ax1.set_xticks(ax1.get_xticks(), ax1.get_xticklabels(), rotation=45)
            ax1.xaxis.set_major_formatter(matplotlib.pyplot.ScalarFormatter())
            ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    handles, labels = ax1.get_legend_handles_labels()
    if find_bounds == True:
        for i in range(3):
            for j in range(3):
                mult=1.02
                ybounds = [mult*ybounds[0], mult*ybounds[1]]
                axes1[i][j].set_ybound(ybounds)
    fig1.legend(handles[0:6], labels[0:6], loc='lower center', bbox_to_anchor=(0.52, 0.), ncols=6, fontsize=12)
    fig1.tight_layout()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.subplots_adjust(bottom=0.12)
    fig1.savefig(f'{output_name}')

    return 0


def plot_stress_norm_calibration_comparison(PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim, E, Ecal, Gamma,
                                            output_name, nqp, increment=None):
    '''Plot the infinity norms of deviatoric Cauchy, symmetric micro, and higher stresses for both homogenized DNS and calibration

    :param dict PK2: The quantities dict storing homogenized DNS second Piola-Kirchhoff stress
    :param dict PK2_sim: The quantities dict storing calibrated second Piola-Kirchhoff stress
    :param dict SIGMA: The quantities dict storing homogenized DNS symmetric micro stress
    :param dict SIGMA_sim: The quantities dict storing calibrated symmetric micro stress
    :param dict M: The quantities dict storing homogenized DNS higher order stress
    :param dict M_sim: The quantities dict storing calibrated higher order stress
    :param dict E: The quantities dict storing homogenized DNS Green-Lagrange strain
    :param dict Ecal: The quantities dict storing homogenized DNS micro strain
    :param dict Gamma: The quantities dict storing homogenized DNS micro-deformation gradient
    :param str output_name: Output filename
    :param int nqp: The number of quadrature points
    :param list increment: An optional list of one or more increments to plot restults

    :returns: ``output_name`` plot
    '''

    fig, axes = matplotlib.pyplot.subplots(1, 3)

    if increment:
        inc = [int(i) for i in increment]
    else:
        inc = [i for i in range(0, numpy.shape(E[0][:,0,0,0])[0])]

    e = 0

    colors = matplotlib.pyplot.rcParams['axes.prop_cycle'].by_key()['color']
    for qp in range(nqp):
        # take norms
        PK2_norm, SIGMA_norm, PK2_norm_sim, SIGMA_norm_sim = [], [], [], []
        M1_norm, M2_norm, M3_norm = [], [], []
        M1_norm_sim, M2_norm_sim, M3_norm_sim = [], [], []
        E_norm, Ecal_norm, Gamma_norm = [], [], []
        for t in inc:
            PK2_norm = numpy.hstack([PK2_norm, deviatoric_norm(PK2[qp][t,e,:,:])])
            PK2_norm_sim = numpy.hstack([PK2_norm_sim, deviatoric_norm(PK2_sim[qp][t,e,:,:])])
            SIGMA_norm = numpy.hstack([SIGMA_norm, deviatoric_norm(SIGMA[qp][t,e,:,:])])
            SIGMA_norm_sim = numpy.hstack([SIGMA_norm_sim, deviatoric_norm(SIGMA_sim[qp][t,e,:,:])])
            M_norms = deviatoric_norm(M[qp][t,e,:,:,:], third_order=True)
            M_norms_sim = deviatoric_norm(M_sim[qp][t,e,:,:,:], third_order=True)
            M1_norm = numpy.hstack([M1_norm, M_norms[0]])
            M1_norm_sim = numpy.hstack([M1_norm_sim, M_norms_sim[0]])
            M2_norm = numpy.hstack([M2_norm, M_norms[1]])
            M2_norm_sim = numpy.hstack([M2_norm_sim, M_norms_sim[1]])
            M3_norm = numpy.hstack([M3_norm, M_norms[2]])
            M3_norm_sim = numpy.hstack([M3_norm_sim, M_norms_sim[2]])
            E_norm = numpy.hstack([E_norm, deviatoric_norm(E[qp][t,e,:,:])])
            Ecal_norm = numpy.hstack([Ecal_norm, deviatoric_norm(Ecal[qp][t,e,:,:])])
            Gamma_norm = numpy.hstack([Gamma_norm, deviatoric_norm(Gamma[qp][t,e,:,:,:], third_order=True, full_norm=True)])

        axes[0].plot(E_norm, PK2_norm, 'o', color=colors[qp])
        axes[0].plot(E_norm, PK2_norm_sim, '-', color=colors[qp])
        axes[1].plot(Ecal_norm, SIGMA_norm, 'o', color=colors[qp])
        axes[1].plot(Ecal_norm, SIGMA_norm_sim, '-', color=colors[qp])
        axes[2].plot(Gamma_norm, M1_norm, 'o', label=f"Filter, K=1", color=colors[qp])
        axes[2].plot(Gamma_norm, M1_norm_sim, '-', label=f"Fit, K=1", color=colors[qp])
        axes[2].plot(Gamma_norm, M2_norm, '^', label=f"Filter, K=2", color=colors[qp])
        axes[2].plot(Gamma_norm, M2_norm_sim, ':', label=f"Fit, K=2", color=colors[qp])
        axes[2].plot(Gamma_norm, M3_norm, 'v', label=f"Filter, K=3", color=colors[qp])
        axes[2].plot(Gamma_norm, M3_norm_sim, '-.', label=f"Fit, K=3", color=colors[qp])

    axes[0].set_ylabel(r'$||dev\left(S_{IJ}\right)|| \left( MPa \right)$', fontsize=14)
    axes[0].set_xlabel(r'$||dev\left(E_{IJ}\right)||$', fontsize=14)

    axes[1].set_ylabel(r'$||dev\left(\Sigma_{IJ}\right)|| \left( MPa \right)$', fontsize=14)
    axes[1].set_xlabel(r'$||dev\left(\mathcal{E}_{IJ}\right)|| $', fontsize=14)

    axes[2].set_ylabel(r'$||dev\left(M_{IJK}\right)|| \left( MPa \cdot mm^2 \right)$', fontsize=14)
    axes[2].set_xlabel(r'$||dev\left(\Gamma_{IJK}\right)||$', fontsize=14)

    fig.set_figheight(5)
    fig.set_figwidth(12)
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(handles[0:6], labels[0:6], loc='lower center', bbox_to_anchor=(0.52, 0.), ncols=6, fontsize=12)
    fig.tight_layout()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.subplots_adjust(bottom=0.2)
    fig.savefig(f'{output_name}')

    return 0


def deviatoric_norm(stress, third_order=False, full_norm=False):
    '''Calculate the norm(s) of the deviatoric part of a stress quantity

    :param array-like stress: A second or third order stress tenosr
    :param bool third_order: A boolean specifying whether the stress tensor is third order
    :param bool full_norm: A boolean specifying whether a third order tensor norm should be across all indices

    :returns: list of norm of deviatoric stresses (1 item for second order stress, 3 for third order)
    '''

    # Third order tenosr norms
    if third_order == True:
        # norm of each "K index"
        if full_norm ==False:
            norm = []
            for k in range(3):
                dev = stress[:,:,k] - (1/3)*numpy.eye(3)*numpy.trace(stress[:,:,k])
                norm.append(numpy.linalg.norm(dev, ord='fro'))
        # Full norm
        else:
            dev = 0
            for k in range(3):
                dev = dev + stress[:,:,k] - (1/3)*numpy.eye(3)*numpy.trace(stress[:,:,k])
            norm = [numpy.linalg.norm(dev, ord='fro')]
    # Norm of second order tensor
    else:
        dev = stress - (1/3)*numpy.eye(3)*numpy.trace(stress)
        norm = [numpy.linalg.norm(dev, ord='fro')]

    return norm


def collect_deviatoric_norm_errors(nqp, t, e, PK2, PK2_sim, SIGMA, SIGMA_sim, M, M_sim):
    '''Calculate the errors between filtered and simulated deviatoric norms for second and third order stresses

    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param int t: The current time increment
    :param int e: The current macro element being considered for calibration
    :param dict PK2: The quantities dict storing a homogenized DNS second Piola-Kirchhoff stress
    :param dict PK2_sim: The quantities dict storing a simulated second Piola-Kirchhoff stress
    :param dict SIGMA: The quantities dict storing a homogenized DNS symmetric micro stress
    :param dict SIGMA_sim: The quantities dict storing a simulated symmetric micro stress
    :param dict M: The quantities dict storing a homogenized DNS higher order stress
    :param dict M_sim: The quantities dict storing a simulated higher order stress

    :returns: error between homogenized and simulations deviatoric second Piola-Kirchhoff, symmetric micro, and higher order stresses
    '''

    PK2_dev_norms = numpy.array([deviatoric_norm(PK2[q][t,e,:,:]) for q in range(nqp)]).flatten()
    PK2_sim_dev_norms = numpy.array([deviatoric_norm(PK2_sim[q][t,e,:,:]) for q in range(nqp)]).flatten()
    SIGMA_dev_norms = numpy.array([deviatoric_norm(SIGMA[q][t,e,:,:]) for q in range(nqp)]).flatten()
    SIGMA_sim_dev_norms = numpy.array([deviatoric_norm(SIGMA_sim[q][t,e,:,:]) for q in range(nqp)]).flatten()
    M_dev_norms = numpy.array([deviatoric_norm(M[q][t,e,:,:,:], third_order=True) for q in range(nqp)]).flatten()
    M_sim_dev_norms = numpy.array([deviatoric_norm(M_sim[q][t,e,:,:,:], third_order=True) for q in range(nqp)]).flatten()

    return PK2_dev_norms-PK2_sim_dev_norms, SIGMA_dev_norms-SIGMA_sim_dev_norms, M_dev_norms-M_sim_dev_norms


def evaluate_constraints(parameters, svals=None):
    '''Evaluate Smith conditions by calling tardigrade_micromorphic_linear_elasticity/src/python/linear_elastic_parameter_constraint_equations

    :param array-like parameters: an array of 18 micromorphic linear elasticity parameters
    :param array-like svals: TODO figure out what this is for

    :returns: a dictionary of constants from evaluating the Smith conditions
    '''
 
    elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                                  'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                                  'tau8', 'tau9', 'tau10', 'tau11']
 
    parameter_dictionary = dict(zip(elastic_parameter_ordering, parameters[:18]))
 
    consts = [constraints.evaluate_g1, 
              constraints.evaluate_g2, 
              constraints.evaluate_g3,
              constraints.evaluate_g4,
              constraints.evaluate_g5,
              constraints.evaluate_g6, 
              constraints.evaluate_g7, 
              constraints.evaluate_g8,
              constraints.evaluate_g9,
              constraints.evaluate_g10,
              constraints.evaluate_g11,
              constraints.evaluate_g12,
              constraints.evaluate_g13]

    if (svals is None):
        svals = dict([(f's{i+1}', 0) for i in range(len(consts))])
 
    parameter_dictionary.update(svals)
 
    return [const(**parameter_dictionary) for const in consts]


def parse_input_parameters(input_parameters):
    ''' Parse material parameters from a YAML file

    :param str input_parameters: YAML file containing calibration results

    :returns: array of elastic parameters and plastic parameters (if there are any)
    '''

    stream = open(input_parameters, 'r')
    UI = yaml.load(stream, Loader=yaml.FullLoader)
    stream.close()
    # elastic parameters case
    if len(UI.keys()) <= 5:
        e_params = numpy.hstack([[2], [float(i) for i in UI['line 1'].split(' ')[1:]],
                                 [5], [float(i) for i in UI['line 2'].split(' ')[1:]],
                                 [11],[float(i) for i in UI['line 3'].split(' ')[1:]],
                                 [2], [float(i) for i in UI['line 4'].split(' ')[1:]]])
        return e_params, None
    else:
        e_params = numpy.hstack([[2], [float(i) for i in UI['line 10'].split(' ')[1:]],
                                 [5], [float(i) for i in UI['line 11'].split(' ')[1:]],
                                 [11],[float(i) for i in UI['line 12'].split(' ')[1:]],
                                 [2], [float(i) for i in UI['line 13'].split(' ')[1:]]])
        p_params = numpy.hstack([[float(i) for i in UI['line 01'].split(' ')[1:]],
                                 [float(i) for i in UI['line 02'].split(' ')[1:]],
                                 [float(i) for i in UI['line 03'].split(' ')[1:]]])
        return e_params, p_params


def parse_fparams_file(parameter_file, material_type='elastic'):
    '''Parse material parameters from a YAML file into an array with parameter names listed

    :param str input_parameters: YAML file containing calibration results
    :param str material_type: The material type: 'elastic', 'plastic', or 'full_plastic'

    :returns: array of parameters and list of parameter names
    '''

    stream = open(parameter_file, 'r')
    UI = yaml.load(stream, Loader=yaml.FullLoader)
    stream.close()

    elastic_parameter_ordering = ['lambda', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                                  'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                                  'tau8', 'tau9', 'tau10', 'tau11']
    plastic_parameter_ordering = ['cu0', 'Hu', 'cchi0', 'Hchi', 'cnablachi0', 'Hnablachi']

    if material_type == 'elastic':
        parameter_ordering = elastic_parameter_ordering + ['obj_func_value']
        params = numpy.hstack([[float(i) for i in UI['line 1'].split(' ')[1:]],
                               [float(i) for i in UI['line 2'].split(' ')[1:]],
                               [float(i) for i in UI['line 3'].split(' ')[1:]],
                               float(UI['obj_func_value'])])
    elif material_type == 'plastic':
        parameter_ordering = plastic_parameter_ordering + elastic_parameter_ordering + ['obj_func_value']
        params = numpy.hstack([[float(i) for i in UI['line 01'].split(' ')[1:]],
                                [float(i) for i in UI['line 02'].split(' ')[1:]],
                                [float(i) for i in UI['line 03'].split(' ')[1:]],
                                [float(i) for i in UI['line 10'].split(' ')[1:]],
                                [float(i) for i in UI['line 11'].split(' ')[1:]],
                                [float(i) for i in UI['line 12'].split(' ')[1:]],
                                float(UI['obj_func_value'])])
    elif material_type == 'full_plastic':
        raise NotImplementedError("'full_plastic' option has not been implemented yet!")
    else:
        raise NameError("Specify a valid material_type!")

    return params, parameter_ordering


def evaluate_model(inputs, parameters, model_name, parameters_to_fparams, nsdvs, element, nqp, maxinc=None, dim=3, maxsubiter=5):
    """Evaluate the model given the parameters. Adapted from overlap_coupling/src/python/read_xdmf_output.py.
   
    :param list inputs: A list storing DNS quantities for Green-Lagrange strain (dict), displacements (dict), displacement gradient (dict), micro-deformation (dict), micro-deformation gradient (dict), and time increments (list)
    :param numpy.ndarray parameters: The array of parameters
    :param str model_name: The name of the model
    :param func parameters_to_fparams: A function that converts the parameters vector to the fparams vector required
        for the function
    :param int nsdvs: The number of solution dependant state variables
    :param int element: The macro (filter) element to calibration
    :param int nqp: The number of quadrature points (1 if filter data is averaged, 8 otherwise)
    :param int maxinc: The maximum increment to evaluate
    :param int dim: The spatial dimension of the problem, default=3
    :param int maxsubiter: The maximum number of sub iterations, default=5

    :returns: evaluated micromorphic simulation quantities for PK2, SIGMA, M, and SDVS
    """

    E, displacement, grad_u, phi, grad_phi, time = inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5]
    ninc = E[0].shape[0]

    nel = 1

    if maxinc is None:
        maxinc = ninc-1

    PK2_sim   = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    SIGMA_sim = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    M_sim     = dict([(qp,numpy.zeros((maxinc+1,nel,dim * dim * dim))) for qp in range(nqp)])
    SDVS_sim  = dict([(qp,numpy.zeros((maxinc+1,nel,nsdvs))) for qp in range(nqp)])
 
    keys = ['errorCode', 'PK2', 'SIGMA', 'M', 'SDVS',\
            'DPK2Dgrad_u', 'DPK2Dphi', 'DPK2Dgrad_phi',\
            'DSIGMADgrad_u', 'DSIGMADphi', 'DSIGMADgrad_phi',\
            'DMDgrad_u', 'DMDphi', 'DMDgrad_phi',\
            'ADD_TERMS', 'ADD_JACOBIANS', 'output_message']

    tp = 0

    nsubiter = 0

    elem = element
    e = 0
    for qp in range(nqp):
        for i in range(maxinc+1):
            #print("increment: ", i)
            # Map the parameters vector to the function parameters
            fparams = parameters_to_fparams(parameters)

            sp = 0
            ds = 1.

            if (i == 0):
                previous_SDVS_s = numpy.zeros(nsdvs)
            else:
                previous_SDVS_s = numpy.copy(SDVS_sim[qp][i-1,e,:])
       
            while (sp < 1.0):

                s = sp + ds

                time_1     = time[i]
                grad_u_1   = grad_u[qp][i, e, :, :]
                phi_1      = phi[qp][i, e, :, :]
                grad_phi_1 = grad_phi[qp][i, e, :, :, :]

                if (i == 0):
                    time_0     = 0
                    grad_u_0   = numpy.zeros((3,3))
                    phi_0      = numpy.zeros((3,3))
                    grad_phi_0 = numpy.zeros((3,3,3))

                else:
                    time_0     = time[i-1]
                    grad_u_0   = grad_u[qp][i-1, e, :, :]
                    phi_0      = phi[qp][i-1, e, :, :]
                    grad_phi_0 = grad_phi[qp][i-1, e, :, :]

                t                = (time_1 - time_0) * s + time_0
                current_grad_u   = (grad_u_1 - grad_u_0) * s + grad_u_0
                current_phi      = (phi_1 - phi_0) * s + phi_0
                current_grad_phi = (grad_phi_1 - grad_phi_0) * s + grad_phi_0

                tp                = (time_1 - time_0) * sp + time_0
                previous_grad_u   = (grad_u_1 - grad_u_0) * sp + grad_u_0
                previous_phi      = (phi_1 - phi_0) * sp + phi_0
                previous_grad_phi = (grad_phi_1 - grad_phi_0) * sp + grad_phi_0

                current_phi = current_phi.flatten()
                previous_phi = previous_phi.flatten()
               
                current_grad_phi = current_grad_phi.reshape((dim * dim, dim))
                previous_grad_phi = previous_grad_phi.reshape((dim * dim, dim))

                #TODO: add dof and add grad dof not currently used
                current_ADD_DOF = numpy.zeros((1))
                current_ADD_grad_DOF = numpy.zeros((1,3))

                previous_ADD_DOF = numpy.zeros((1))
                previous_ADD_grad_DOF = numpy.zeros((1,3))

                # Evaluate the model
                values = micromorphic.evaluate_model(model_name, numpy.array([t, t - tp]), fparams,
                                                     current_grad_u, current_phi, current_grad_phi,
                                                     previous_grad_u, previous_phi, previous_grad_phi,
                                                     previous_SDVS_s,
                                                     current_ADD_DOF, current_ADD_grad_DOF,
                                                     previous_ADD_DOF, previous_ADD_grad_DOF)

                results = dict(zip(keys, values))

                if (results['errorCode'] == 1):
                    #print("error")
                    ds = 0.5 * ds
                    nsubiter += 1

                    if (nsubiter > maxsubiter):
                        break

                elif (results['errorCode'] == 2):
                    errormessage = f"evaluate_model return error code {results['errorCode']}\n\n"
                    errormessage += results['output_message'].decode("utf-8")
                    raise IOError(errormessage)

                else:
                    sp += ds
                    nsubiter = 0

                    if numpy.isclose(sp, 1):
                        ds = 1
                    else:
                        ds = 1 - sp

                    previous_SDVS_s = numpy.copy(results['SDVS'])

            if (results['errorCode'] != 0):
                errormessage = f"evaluate_model returned error code {results['errorCode']}\n\n"
                errormessage += results['output_message'].decode('utf-8')
                print(parameters, 'fail')

                return numpy.nan

            PK2_sim[qp][i,e,:]   = results['PK2']
            SIGMA_sim[qp][i,e,:] = results['SIGMA']
            M_sim[qp][i,e,:]     = results['M']
            SDVS                 = results['SDVS']
            SDVS_sim[qp][i,e,:]  = results['SDVS']

    return PK2_sim, SIGMA_sim, M_sim, SDVS_sim
