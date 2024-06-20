# ------------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------------

import sys
import os
import argparse
import inspect
import datetime
import yaml

import numpy as np
import h5py
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import scipy
from itertools import compress

import micromorphic
import xdmf_reader_tools as XRT
import linear_elastic_parameter_constraint_equations as constraints

elastic_parameter_ordering = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                              'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                              'tau8', 'tau9', 'tau10', 'tau11']

def average_quantities(quantities, type):
    output = {}
    output[0] = np.zeros_like(quantities[0])
    
    if type == '9':
        for k in range(9):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,0,k])
            means = np.mean(mean_field, axis=0)
            output[0][:,0,k] = means
    elif type == '3x3':
        for i in range(3):
            for j in range(3):
                mean_field = []
                for qp in quantities.keys():
                    mean_field.append(quantities[qp][:,0,i,j])
                means = np.mean(mean_field, axis=0)
                output[0][:,0,i,j] = means
    elif type == '3':
        for i in range(3):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,i,0])
            means = np.mean(mean_field, axis=0)
            output[0][:,i,0] = means
            
    return(output)


# Functions
def plot_stresses(e, stress, stress_sim, options, output_name):
    name = output_name.replace('.PNG','')
    fig1 = plt.figure(name, figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [-1, 1]
 
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    k = 0
    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
           
            for qp in stress.keys():
                if 'cauchy' in output_name:
                    plot_label = r"$\sigma_{" + str(i+1) + str(j+1) + "}$ (MPa)"
                    #stress = cauchy[qp][:,0,i,j]
                if 'symm' in output_name:
                    plot_label = r"$s_{" + str(i+1) + str(j+1) + "}$ (MPa)"
                if 'PK2' in output_name:
                    plot_label = r"$PK2_{" + str(i+1) + str(j+1) + "}$ (MPa)"
                if 'SIGMA' in output_name:
                    plot_label = r"$SIGMA_{" + str(i+1) + str(j+1) + "}$ (MPa)"
                    #stress = symm[qp][:,0,i,j]
                if 'average' in options:
                    ax1.plot(e[qp][:,0,i,j], stress[qp][:,0,i,j], '-')
                    ax1.plot(e[qp][:,0,i,j], stress_sim[qp][:,0,i,j], 'o')
                else:
                    ax1.plot(e[qp][:,0,i,j], stress[qp][:,0,i,j], '-', color=colors[qp], label=f'qp #{qp+1}')
                    ax1.plot(e[qp][:,0,i,j], stress_sim[qp][:,0,i,j], 'o', color=colors[qp], label=f'qp #{qp+1}')
                #ybounds[0] = min(PK2[qp][:,0,i,j].min() - 10, ybounds[0])
                #ybounds[1] = min(PK2[qp][:,0,i,j].min() - 10, ybounds[1])
 
            ax1.set_xlabel(r"$e_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.set_ylabel(plot_label, fontsize=14)
            #plt.ticklabel_format(style='sci', axis='x', scilimits=(-10,10))
            #plt.ticklabel_format(style='sci', axis='y', scilimits=(-10,10))
            plt.ticklabel_format(style='sci', axis='x')
            plt.ticklabel_format(style='sci', axis='y')
            if (i == 2) and (j ==2):
                plt.xticks(rotation=45)
            k = k + 1
 
    ax1 = axes1[1][2]
    if 'average' not in options:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.figure(name)
    plt.tight_layout()
    #plt.show()
    plt.savefig(f'{output_name}')
 
    return 0
   
def initial_estimate(Emod, nu, L):
 
    print(f"E, nu, L = {Emod}, {nu}, {L}")
 
    # calculate "classic" lame parameters
    lame_lambda = Emod*nu/((1.+nu)*(1.-2*nu))
    lame_mu     = Emod/(2*(1.+nu)) #shear modulus, K
   
    # estimate characteristic length
    Lc = np.sqrt(3*(L**2))
   
    # estimate micromorphic parameters
    lamb = lame_lambda
    #lamb = 59.0
    mu = lame_mu
    #mu = 70.0
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
   
    # collect
    parameters = np.array([lamb, mu, eta, tau, kappa, nu_new, sigma,
                           tau_1, tau_2, tau_3, tau_4, tau_5, tau_6,
                           tau_7, tau_8, tau_9, tau_10, tau_11])
                          
    return(parameters)
    
def evaluate_constraints(parameters, svals=None):
 
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
   
def evaluate_model(inputs, parameters, model_name, parameters_to_fparams, nsdvs, maxinc=None, dim=3, maxsubiter=5):
    """
    Evaluate the model given the parameters. Copied from overlap_coupling/src/python/read_xdmf_output.py.
   
    :param dict data: The data dictionary
    :param np.ndarray parameters: The array of parameters
    :param str model_name: The name of the model
    :param func parameters_to_fparams: A function that converts the parameters vector to the fparams vector required
        for the function
    :param int nsdvs: The number of solution dependant state variables
    :param int maxinc: The maximum increment to evaluate
    :param int dim: The spatial dimension of the problem
    :param int maxsubiter: The maximum number of sub iterations
    """
    #position, grad_u, phi, grad_phi = XRT.construct_degrees_of_freedom(data, dim = dim)
    E, position, grad_u, phi, grad_phi, time = inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5]
    ninc = E[0].shape[0]
 
    nel = E[0].shape[1]
 
    nqp = len([key for key in E.keys()])
  
    if maxinc is None:
        maxinc = ninc-1
 
    PK2_sim   = dict([(qp,np.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    SIGMA_sim = dict([(qp,np.zeros((maxinc+1,nel,dim * dim))) for qp in range(nqp)])
    M_sim     = dict([(qp,np.zeros((maxinc+1,nel,dim * dim * dim))) for qp in range(nqp)])
    SDVS_sim  = dict([(qp,np.zeros((maxinc+1,nel,nsdvs))) for qp in range(nqp)])
 
    keys = ['errorCode', 'PK2', 'SIGMA', 'M', 'SDVS',\
            'DPK2Dgrad_u', 'DPK2Dphi', 'DPK2Dgrad_phi',\
            'DSIGMADgrad_u', 'DSIGMADphi', 'DSIGMADgrad_phi',\
            'DMDgrad_u', 'DMDphi', 'DMDgrad_phi',\
            'ADD_TERMS', 'ADD_JACOBIANS', 'output_message']
   
    #time = data['time'][1:]
 
    tp = 0
 
    nsubiter = 0
       
    for e in range(nel):
        for qp in range(nqp):
            for i in range(maxinc+1):
                #print("increment: ", i)
                # Map the parameters vector to the function parameters
                fparams = parameters_to_fparams(parameters)
           
                sp = 0
                ds = 1.
               
                if (i == 0):
                    previous_SDVS_s = np.zeros(nsdvs)
                else:
                    previous_SDVS_s = np.copy(SDVS_sim[qp][i-1,e,:])
           
                while (sp < 1.0):
               
                    s = sp + ds
               
                    time_1     = time[i]
                    grad_u_1   = grad_u[qp][i, e, :, :]
                    phi_1      = phi[qp][i, e, :, :]
                    grad_phi_1 = grad_phi[qp][i, e, :, :, :]
                   
                    if (i == 0):
                        time_0     = 0
                        grad_u_0   = np.zeros((3,3))
                        phi_0      = np.zeros((3,3))
                        grad_phi_0 = np.zeros((3,3,3))
                        
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
                    current_ADD_DOF = np.zeros((1))
                    current_ADD_grad_DOF = np.zeros((1,3))
 
                    previous_ADD_DOF = np.zeros((1))
                    previous_ADD_grad_DOF = np.zeros((1,3))
                   
                    # Evaluate the model
                    values = micromorphic.evaluate_model(model_name, np.array([t, t - tp]), fparams,\
                                                         current_grad_u, current_phi, current_grad_phi,
                                                         previous_grad_u, previous_phi, previous_grad_phi,
                                                         previous_SDVS_s,\
                                                         current_ADD_DOF, current_ADD_grad_DOF,
                                                         previous_ADD_DOF, previous_ADD_grad_DOF)
                   
                    results = dict(zip(keys, values))
 
#                    print(results['output_message'].decode("utf-8"))
 
                    #print("    nsubiter, ds: ", nsubiter, ds)
                   
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
                       
                        if np.isclose(sp, 1):
                            ds = 1
                        else:
                            ds = 1 - sp
                           
                        previous_SDVS_s = np.copy(results['SDVS'])
                        
                if (results['errorCode'] != 0):
                    errormessage = f"evaluate_model returned error code {results['errorCode']}\n\n"
                    errormessage += results['output_message'].decode('utf-8')
                    print(parameters, 'fail')
                   
                    return np.nan
               
                PK2_sim[qp][i,e,:]   = results['PK2']
                SIGMA_sim[qp][i,e,:] = results['SIGMA']
                M_sim[qp][i,e,:]     = results['M']
                SDVS                 = results['SDVS']
                SDVS_sim[qp][i,e,:]  = results['SDVS']
               
    return PK2_sim, SIGMA_sim, M_sim, SDVS_sim
   
def parameters_to_fparams(parameters):
    """
    Map the elastic parameters to the fparams vector for use in the micromorphic material models
 
    :param np.ndarray parameters: The parameters vector
        lambda, mu, eta, tau, kappa, nu, sigma, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11
    """
   
    fparams = np.hstack([[2], parameters[:2], [5], parameters[2:7], [11], parameters[7:18], 2, [parameters[3], parameters[6]]])
 
    return fparams
 
def objective(x0, Y, inputs, cal_norm, nu_targ, options, return_error=False, stresses_to_include=['S','SIGMA','M']):
 
    Xstore = []
    Lstore = []
    Fstore = []
   
    parameter_names = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',
                   'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',
                   'tau8', 'tau9', 'tau10', 'tau11']
 
    model_name=r'LinearElasticity'
    XX = x0
    #consts = evaluate_constraints(XX, parameter_names)
    consts = evaluate_constraints(XX)
    cvals = np.array([c[0] for c in consts])
   
    # enforce positivity condition
    if (cvals.min() < 0):
        if return_error:
            return -np.inf, np.inf
        else:
            return np.inf
            #return 1.e16
            
    # enforce Poisson ratio constraint within 2%
    if '1' in options:
        E = inputs[0]
        poisson = XX[0]/(2.*(XX[0] + XX[1]))
        #poisson = np.average([E[0][-1,0,0,0],E[0][-1,0,1,1]])
        if (poisson >= 1.01*nu_targ) or (poisson <= .99*nu_targ):
            return np.inf
           
    # Evaluate stresses from DNS strain inputs
    PK2_sim, SIGMA_sim, M_sim, SDVS_sim = evaluate_model(inputs, XX, model_name, parameters_to_fparams, 0)
    position, grad_u, phi, grad_phi = inputs[1], inputs[2], inputs[3], inputs[4]
    
    # convert flattened PK2 and SIGMA to 3x3 structure, ignore m and M for now!
    #PK2_sim = XRT.map_sim(PK2_sim)
    #SIGMA_sim = XRT.map_sim(SIGMA_sim)
   
    # Parse out stresses from DNS stress data Y
    PK2, SIGMA, M = Y[0], Y[1], Y[2]
   
    # map all stresses to current configuration, ignore m and M for now!
    #cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, position, grad_u, phi, grad_phi)
    #cauchy, symm = XRT.get_current_configuration_stresses(PK2, SIGMA, position, grad_u, phi, grad_phi)
   
    # Number of time steps and elements
    steps = PK2[0].shape[0]
    num_elem = PK2[0].shape[1]
   
    # Initialize errors and objective
    PK2_error   = []
    SIGMA_error = []
    M_error     = []
    obj = 0
   
   # define time steps to calibrate against
    if 'final' in options:
        time_steps = [steps-1]
    else:
        time_steps = range(steps)
   
    # Loop over quadrature points, time, and elements
    if '1' in options:
        for qp in PK2.keys():
            for i in range(1, PK2[qp].shape[0]):
                for e in range(num_elem):
                    PK2_error = np.hstack([PK2_error, PK2[qp][i,e,2,2] - PK2_sim[qp][i,e,-1]])
                    SIGMA_error = np.hstack([PK2_error, PK2[qp][i,e,2,2] - PK2_sim[qp][i,e,-1]])
                    M_error = []
    else:
        for qp in PK2.keys():
            for i in range(1, PK2[qp].shape[0]):
                for e in range(num_elem):
                    PK2_error = np.hstack([PK2_error, PK2[qp][i,e,:,:].flatten() - PK2_sim[qp][i,e,:]])
                    SIGMA_error = np.hstack([SIGMA_error, SIGMA[qp][i,e,:,:].flatten() - SIGMA_sim[qp][i,e,:]])
                    M_error = np.hstack([M_error, M[qp][i,e,:,:].flatten() - M_sim[qp][i,e,:]])
               
    # collect errors
    errors    = {'S':PK2_error, 'SIGMA':SIGMA_error, 'M':M_error}
   
    sparsity_control = 0.1
    # calculate residual, L1 norm
    for stress in stresses_to_include:
        if cal_norm == 'L1':
            obj += np.abs(errors[stress]).sum()
        elif cal_norm == 'L2':
            obj += np.dot(errors[stress], errors[stress])
        elif cal_norm == 'L1_sparse':
            obj += np.abs(errors[stress]).sum() + sparsity_control*SOMETHING
            print('NOT READY YET!')
        elif cal_norm == 'L2_sparse':
            obj += np.dot(errors[stress], errors[stress]) + sparsity_control*SOMETHING
            print('NOT READY YET!')
        elif cal_norm == 'L1-L2':
            obj += 0.5*(np.abs(errors[stress]).sum()) + 0.5*(np.dot(errors[stress], errors[stress]))
        else:
            print('Specify valid objective!')
       
    Xstore.append(np.copy(x0))
    Lstore.append(obj)
   
    print(f'obj = {obj}')
    return obj


def opti_options_1(X,Y, inputs, cal_norm, nu_targ, options, calibrate=True):
    X1 = X[:2]
    others = [0e-3, 0e-3, 0e-3, 0e-3, 0e-3,
              0e-3, 0e-3, 0e-3, 0e-3, 0, 0, 1e-3, 0e-3, 0, 0e-3, 0e-3]
    XX = np.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, options, stresses_to_include=['S','SIGMA']))
    else:
        return(np.hstack([X1, others]))


def opti_options_2(X,Y, inputs, cal_norm, nu_targ, options, calibrate=True):
    X1 = X[:7]
    others = [0e-3, 0e-3, 0e-3, 0e-3, 0, 0, 1e-3, 0e-3, 0, 0e-3, 0e-3]
    XX = np.hstack([X1, others])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, options, stresses_to_include=['S', 'SIGMA']))
    else:
        return(np.hstack([X1, others]))


def opti_options_3(X,Y, inputs, cal_norm, nu_targ, options, calibrate=True):
    X1, X2 = X[:7], X[-1]
    tau1to6 = [0e-3, 0e-3, 0e-3, 0e-3, 0, 0]
    tau8to11 = [0e-3, 0, 0e-3, 0e-3]
    XX = np.hstack([X1, tau1to6, X2, tau8to11])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, options, stresses_to_include=['S', 'SIGMA']))
    else:
        return(np.hstack([X1, tau1to6, X2, tau8to11]))


def opti_options_4(X,Y, inputs, cal_norm, nu_targ, options, calibrate=True):
    XX = X
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, options, stresses_to_include=['S', 'SIGMA', 'M']))
    else:
        return(XX)


def opti_options_5(X,Y, inputs, cal_norm, options, calibrate=True):
    X1, X2 = X[:5], X[-1]
    tau1to6 = [0e-3, 0e-3, 0e-3, 0e-3, 0, 0]
    tau8to11 = [0e-3, 0, 0e-3, 0e-3]
    XX = np.hstack([[59.25348, 70.39516], X1, tau1to6, X2, tau8to11])
    if calibrate:
        return(objective(XX, Y, inputs, cal_norm, nu_targ, options, stresses_to_include=['S', 'SIGMA']))
    else:
        return(np.hstack([[59.25348, 70.39516], X1, tau1to6, X2, tau8to11]))


def calibrate(input_file, output_file, options, Emod, nu, L):
    # Main
    PK2_sdstore = []
    SIGMA_sdstore = []
    M_sdstore = []
    Xstore = []
    Lstore = []
    Fstore = []
 
    # Read in the data
    filename = input_file
    data, geometry, topology = XRT.parse_xdmf_output(filename)
    nqp = 8
    ninc = np.shape(data['time'])[0]
    nel = np.shape(data['cauchy_stress_0_0'])[0]
 
    # Read in the position information
    position, gradu, phi, gradphi = XRT.construct_degrees_of_freedom(data, nqp, nel)

    # Read in the stress information
    cauchy, symm, m = XRT.collect_stresses(data, nqp, nel)
    PK2, SIGMA, M = XRT.get_reference_configuration_stresses(data, nqp, nel)

    # Read in the strain information
    E, Ecal, Gamma, F, chi, grad_chi, e, h = XRT.compute_deformations(data, nqp, nel)
    
    # Get times
    times = np.unique(data['time'])
    ninc = len(times)
    
    # average fields if selected as an option
    if 'average' in options:
        cauchy      = average_quantities(cauchy, '3x3')
        symm        = average_quantities(symm, '3x3')
        PK2         = average_quantities(PK2, '3x3')
        SIGMA       = average_quantities(SIGMA, '3x3')
        #F           = average_quantities(F, '3x3')
        E           = average_quantities(E, '3x3')
        #position    = average_quantities(position, '3')
        gradu       = average_quantities(gradu, '3x3')
        phi         = average_quantities(phi, '3x3')
        e           = average_quantities(e, '3x3')
        h           = average_quantities(h, '3x3')
        gradphi     = average_quantities(gradphi, '3x3')
 
    # store data for calibration
    Y = [PK2, SIGMA, M]
    inputs = [E, position, gradu, phi, gradphi, times]
    
    # get target nu from E
    #nu_targ = -1*np.average([E[-1][-1,0,0,0],E[-1][-1,0,1,1]])/E[-1][-1,0,2,2]
    nu_targ = np.average([-1*np.average([E[key][-1,0,0,0],E[key][-1,0,1,1]])/E[key][-1,0,2,2] for key in E.keys()])
    print(nu_targ)
 
    # Estimate initial parameters
    param_est = initial_estimate(Emod, nu, L)
 
    # Define the elastic bounds
    parameter_bounds = [[0.0, 1.e5]] + [[-1.e5, 1.e5] for _ in range(6)] + [[-1.e5, 1.e5] for _ in range(11)]
 
    # define the distributions
    #distributions = [scipy.stats.uniform(loc=lb, scale=(ub-lb) for lb, ub in parameter_bounds]
 
    # compute a test parameter
    parameter_names = ['l', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',
                       'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',
                       'tau8', 'tau9', 'tau10', 'tau11']
    if len(parameter_names) != len(parameter_bounds):
        raise ValueError(f"The parameter bounds and the parameter names do not have the same length {len(parmaeter_bounds)} vs. {len(parameter_names)}")

    # calibrate!
    cal_norm = 'L1'
    maxit = 2000
    # calibrate just lambda and mu
    if '1' in options:
        param_mask = [True, True, False, False, False, False, False,
                      False, False, False, False, False, False,
                      False, False, False, False, False]
        param_est = list(compress(param_est, param_mask))
        print('initial parameter estimation:')
        print(param_est)
        #param_est = [59.25, 70.395]
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_1, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, options))
        print(f"res = {res}")
        print(f"fit params = {list(res.x)}")
        params = opti_options_1(list(res.x), Y, inputs, cal_norm, nu_targ, options, calibrate=False)
    # calibrate first 7 parameters
    elif '2' in options:
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      False, False, False, False, False] 
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_2, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, options))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_2(res.x, Y, inputs, cal_norm, nu_targ, options, calibrate=False)
        #fparams = np.array([
    # calibrate first 7 parameters and tau 7
    elif '3' in options:
        param_mask = [True, True, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False] 
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_3, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, options))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_3(res.x, Y, inputs, cal_norm, nu_targ, options, calibrate=False)
    # calibrate all parameters
    elif '4' in options:
        param_mask = [True, True, True, True, True, True, True,
                      True, True, True, True, True, True, True, True, True, True, True]
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_4, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, options))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_4(res.x, Y, inputs, cal_norm, nu_targ, options, calibrate=False)
    # calibrate first seven parameters and tau 7 except force specifc values for lamba and mu
    elif '5' in options:
        param_mask = [False, False, True, True, True, True, True,
                      False, False, False, False, False, False,
                      True, False, False, False, False] 
        param_est = list(compress(param_est, param_mask))
        parameter_bounds = list(compress(parameter_bounds,param_mask))
        res = scipy.optimize.differential_evolution(func=opti_options_5, bounds=parameter_bounds, maxiter=maxit, x0=param_est, args=(Y, inputs, cal_norm, nu_targ, options))
        print(f"res = {res}")
        print(f"fit params = {res.x}")
        params = opti_options_5(res.x, Y, inputs, cal_norm, nu_targ, options, calibrate=False)
    else:
        print('Select valid option!')
   
    #print(f"res = {res}")
    #print(f"params = {res.x}")
   
    # plot resulting calibration
    if 'plot' in options:
        print('plotting...')
        model_name=r'LinearElasticity'
        PK2_sim, SIGMA_sim, M_sim, SDVS_sim = evaluate_model(inputs, params, model_name, parameters_to_fparams, 0)
        PK2_sim = XRT.map_sim(PK2_sim, ninc)
        SIGMA_sim = XRT.map_sim(SIGMA_sim, ninc)
        cauchy_sim, symm_sim = XRT.get_current_configuration_stresses(PK2_sim, SIGMA_sim, inputs[1], inputs[2], inputs[3], inputs[4])
        #print(PK2_sim)
        #print(np.shape(PK2_sim[0]))

        # average if in options
        if 'average' in options:
            PK2_sim     = average_quantities(PK2_sim, '3x3')
            SIGMA_sim   = average_quantities(SIGMA_sim, '3x3')
            cauchy_sim  = average_quantities(cauchy_sim, '3x3')
            symm_sim    = average_quantities(symm_sim, '3x3')

        plot_stresses(e, cauchy, cauchy_sim, options, 'cauchy_fit.PNG')
        plot_stresses(e, symm, symm_sim, options, 'symm_fit.PNG')

    # output parameters!
    output_filename = output_file
    output_dict = {}
    #output_dict['res'] = str(res)
    #output_dict['fit_params'] = [str(item) for item in res.x]
    p = params
    output_dict['line 1'] = f"2 {p[0]} {p[1]}"
    output_dict['line 2'] = f"5 {p[2]} {p[3]} {p[4]} {p[5]} {p[6]}"
    output_dict['line 3'] = f"11 {p[7]} {p[8]} {p[9]} {p[10]} {p[11]} {p[12]} {p[13]} {p[14]} {p[15]} {p[16]} {p[17]}"
    output_dict['line 4'] = f"2 {p[3]} {p[6]}"
    with open(output_filename, 'w') as f:
        yaml.dump(output_dict, f)
    return


def get_parser():
 
    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Calibrate "
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str,
        help="The homogenized XDMF file output by the Micromorphic Filter")
    parser.add_argument('-o', '--output-file', type=str,
        help="The resulting list of parameters stored in a ____ file")
    parser.add_argument('--options', type=str, required=False, default='',
        help="Specify the calibration options")
    parser.add_argument('--Emod', type=float,
        help="DNS elastic modulus, used for initial parameter estimation.")
    parser.add_argument('--nu', type=float,
        help="DNS Poisson's ratio, used for initial parameter estimation.")
    parser.add_argument('--L', type=float,
        help="DNS max dimension (width, height, depth, etc.), used for initial parameter estimation.")
       
    return parser


if __name__ == '__main__':
    parser = get_parser()
   
    args, unknown = parser.parse_known_args()
    sys.exit(calibrate(
                    input_file=args.input_file,
                    output_file=args.output_file,
                    options=args.options,
                    Emod=args.Emod,
                    nu=args.nu,
                    L=args.L))