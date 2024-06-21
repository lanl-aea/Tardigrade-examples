import sys
import inspect
import os
import argparse

import numpy
import scipy.stats
import scipy.optimize
import h5py
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import h5py
import pandas
from scipy.linalg import norm
from scipy.linalg import polar
import argparse

import xdmf_reader_tools as XRT


def str2bool(v):
    '''Function for converting string to Boolean. Borrowed from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse

    :param str/bool v: A string or boolean indicating a True or False value

    :returns: True or False
    '''

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def plot_reference_stresses(E, stress, average, output_name):
    ''' Plot stress vs strain in the reference configuration

    :param dict E: The quantities dict storing Green-Lagrance strain
    :param dict stress: The quantities dict storing either second Piola Kirchhoff or Symmetric micro stress
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    name = output_name.replace('.PNG','')
    fig1 = plt.figure(name, figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [-1, 1]

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            if 'PK2' in output_name:
                plot_label = r"$PK2_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if 'SIGMA' in output_name:
                plot_label = r"$\Sigma_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if average:
                ax1.plot(E[0][:,0,i,j], stress[0][:,0,i,j], '-o')
            else:
                for qp in stress.keys():
                    ax1.plot(E[qp][:,0,i,j], stress[qp][:,0,i,j], '-o', color=colors[qp], label=f'qp #{qp+1}')
            if (i == 2) and (j == 2):
                plt.xticks(rotation=45)
                        
            ax1.set_xlabel(r"$E_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.set_ylabel(plot_label, fontsize=14)
            plt.ticklabel_format(style='sci', axis='x')
            plt.ticklabel_format(style='sci', axis='y')
            ax1.tick_params('x', labelrotation=45)

    ax1 = axes1[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.figure(name)
    plt.tight_layout()
    plt.savefig(f'{output_name}')

    return 0


def plot_current_stresses(e, stress, average, output_name):
    ''' Plot stress vs strain in the current configuration

    :param dict e: The quantities dict storing Euler-Almansi strain
    :param dict stress: The quantities dict storing either Cauchy or symmetric micro stress
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    name = output_name.replace('.PNG','')
    fig1 = plt.figure(name, figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    ybounds = [-1, 1]

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            if 'cauchy' in output_name:
                plot_label = r"$\sigma_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if 'symm' in output_name:
                plot_label = r"$s_{" + str(i+1) + str(j+1) + "}$ (MPa)"
            if average:
                ax1.plot(e[0][:,0,i,j], stress[0][:,0,i,j], '-o')
            else:
                for qp in stress.keys():
                    ax1.plot(e[qp][:,0,i,j], stress[qp][:,0,i,j], '-o', color=colors[qp], label=f'qp #{qp+1}')
            if (i == 2) and (j ==2):
                plt.xticks(rotation=45)

            ax1.set_xlabel(r"$e_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.set_ylabel(plot_label, fontsize=14)
            plt.ticklabel_format(style='sci', axis='x')
            plt.ticklabel_format(style='sci', axis='y')
            ax1.tick_params('x', labelrotation=45)

    ax1 = axes1[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.figure(name)
    plt.tight_layout()

    plt.savefig(f'{output_name}')

    return 0


def plot_axial_vector_mag(cauchy, times, average, output_name):
    ''' Plot the magnitude of the Cauchy couple

    :param dict cauchy: The quantities dict storing Cauchy stress
    :param array-like times: The time increments
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    fig1 = plt.figure('axial', figsize=(3,5),dpi=300)
    
    for qp in cauchy.keys():
        axial = []
        for i,t in enumerate(times):
            # skew symmetric 
            sig = cauchy[qp][i,0,:,:]
            sig_t = numpy.transpose(sig)
            skew = 0.5*(sig - sig_t)

            # off diagonals
            ax = numpy.sqrt((skew[0][1]**2)+(skew[0][2]**2)+(skew[1][2]**2))

            # Append results
            axial.append(ax)

        # Plot magnitude of axial vector
        if average:
            plt.plot(times, axial, '-o')
        else:
            plt.plot(times, axial, '-o', color=colors[qp], label=f'qp #{qp+1}')
        
    plt.figure('axial')
    plt.xlabel('Simulation time (s)', fontsize=14)
    plt.ylabel('|$\omega^{\sigma}$| (MPa)', fontsize=14)
    if average == False:
        plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_name}')

    return 0 


def plot_stress_diffs(cauchy, symm, times, average, output_name):
    ''' Plot differences between Cauchy and micro symmetric stresses

    :param dict cauchy: The quantities dict storing Cauchy stress
    :param dict symm: The quantities dict storing symmetric micro stress
    :param array-like times: The time increments
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    fig1 = plt.figure('stress_diff', figsize=(11,9),dpi=300)
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            
            for qp in cauchy.keys():
                sig = cauchy[qp][:,0,i,j]
                sym = symm[qp][:,0,i,j]
                diff = sig-sym
                if average:
                    ax1.plot(times, diff, '-o')
                else:
                    ax1.plot(times, diff, '-o', color=colors[qp], label=f'qp #{qp+1}')

            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax1.set_ylabel(r"($\sigma_{" + str(i+1) + str(j+1) + "}$ - $s_{" + str(i+1) + str(j+1) + "}$)", fontsize=14)
    ax1 = axes1[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.figure('stress_diff')
    plt.tight_layout()
    plt.savefig(f'{output_name}')

    return 0


def get_R_and_U(PK2, F, chi, times):
    ''' Calculate macro and micro stretches and deformations from deformation gradients

    :param dict PK2: The quantities dict storing second Piola Kirchhoff stress
    :param dict F: The quantities dict storing the macro deformation gradient
    :param dict chi: The quantities dict storing the micro deformation
    :param array-like times: The time increments

    :returns: macro rotation (R), macro stretch (U), micro rotation (Rchi), and micro stretch (Uchi)
    '''

    R, U = {}, {}
    Rchi, Uchi = {}, {}
    for qp in PK2.keys():
        Rs, Us = [], []
        Rcs, Ucs = [], []
        for i,t in enumerate(times):
            #if i > 0:
            Rr, Uu = polar(F[qp][i,0,:,:],side='left')
            Rcr, Ucu = polar(chi[qp][i,0,:,:], side='left')

            #check orthogonality
            tol = 1.e-9
            if (norm((numpy.dot(Rr,numpy.transpose(Rr)) - numpy.eye(3)),ord=2)) > tol:
                print('Error!!! R is not orthogonal!')
            if (norm((numpy.dot(Rcr,numpy.transpose(Rcr)) - numpy.eye(3)),ord=2)) > tol:
                print('Error!!! Rc is not orthogonal!')

            Rs.append(Rr)
            Us.append(Uu)
            Rcs.append(Rcr)
            Ucs.append(Ucu)
        R[qp] = Rs
        U[qp] = Us
        Rchi[qp] = Rcs
        Uchi[qp] = Ucs

    return(R,U,Rchi,Uchi)

    
def plot_rot_diffs(PK2, times, R, Rchi, average, output_name):
    ''' Plot differences between macro and micro rotations

    :param dict PK2: The quantities dict storing second Piola Kirchhoff stress
    :param array-like times: The time increments
    :param dict R: The quantities dict storing macro rotations
    :param dict Rchi: The quantities dict storing micro rotations
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    fig1 = plt.figure('Rot_diff', figsize=(12,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            
            for qp in PK2.keys():
                diffs = []
                #for k in range(len(times)-1):
                for k,t in enumerate(times):
                    diff = R[qp][k][i][j]-Rchi[qp][k][i][j]
                    diffs.append(diff)
                if average:
                    ax1.plot(times[:], diffs, '-o')
                else:
                    ax1.plot(times[:],diffs, '-o', color=colors[qp], label=f'qp #{qp+1}')
                
            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax1.set_ylabel(r"$R_{" + str(i+1) + str(j+1) + "}$ - $R^{\chi}_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.ticklabel_format(axis='y',style='sci')
    ax1 = axes1[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.ticklabel_format(axis='y',style='sci')
    plt.figure('Rot_diff')
    plt.tight_layout()
    plt.savefig(f'{output_name}')
    
    return 0 


def plot_stretch_diffs(PK2, times, U, Uchi, average, output_name):
    ''' Plot differences between macro and micro stretches

    :param dict PK2: The quantities dict storing second Piola Kirchhoff stress
    :param array-like times: The time increments
    :param dict U: The quantities dict storing macro stretches
    :param dict Uchi: The quantities dict storing micro stretches
    :param bool average: Average over quadrature points if True
    :param str output_name: The output plot name

    :returns: ``output_name``
    '''

    fig1 = plt.figure('U', figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            
            for qp in PK2.keys():
                diffs = []
                for k,t in enumerate(times):
                    diff = U[qp][k][i][j]-Uchi[qp][k][i][j]
                    diffs.append(diff)
                if average:
                    ax1.plot(times[:], diffs, '-o')
                else:
                    ax1.plot(times[:],diffs, '-o', color=colors[qp], label=f'qp #{qp+1}')

                
            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax1.set_ylabel(r"$U_{" + str(i+1) + str(j+1) + "}$ - $U^{\chi}_{" + str(i+1) + str(j+1) + "}$", fontsize=14)
            ax1.ticklabel_format(axis='y',style='sci')
    ax1 = axes1[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2, 0.9))
    plt.figure('U')
    plt.tight_layout()
    plt.savefig(f'{output_name}')
    
    return 0 


def plot_first_moment_of_momentum_measures(coup, spins, times, average, plot_body_couples=None, plot_spin_inertias=None, plot_spin_diff=None):
    '''Plot the first moment of momentum measures

    :param dict coup: The quantities dict storing body couples
    :param dict spins: THe quantities dict stroring micro spin inertias
    :param array-like times: The time increments
    :param bool average: Average over quadrature points if True
    :param str plot_body_couples: Optional filename to plot body couples vs. simulation time
    :param str plot_spin_inertias: Optional filename to plot micro spin inertias vs. simulation time
    :param str plot_spin_diff: Optional filenmae to plot difference between body couples and micro spin inertias

    :returns: ``plot_body_couples``, ``plot_spin_inertias``, and ``plot_spin_diff``
    '''
    
    name1 = 'body_couples'
    name2 = 'micro_spin_inertias'
    name3 = 'couple_spin_diff'
    fig1 = plt.figure(name1, figsize=(11,9))
    axes1 = [[fig1.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    fig2 = plt.figure(name2, figsize=(11,9))
    axes2 = [[fig2.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    fig3 = plt.figure(name3, figsize=(11,9))
    axes3 = [[fig3.add_subplot(3,3,3 * i + j + 1) for j in range(3)] for i in range(3)]
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    for i in range(3):
        for j in range(3):
            ax1 = axes1[i][j]
            ax2 = axes2[i][j]
            ax3 = axes3[i][j]
            
            if average:
                ax1.plot(times[:], coup[0][:,0,i,j], '-o')
                ax2.plot(times[:], spins[0][:,0,i,j], '-o')
                ax3.plot(times[:], coup[0][:,0,i,j] - spins[0][:,0,i,j], '-o')
            else:
                for qp in coup.keys():
                    ax1.plot(times[:], coup[qp][:,0,i,j], '-o', color=colors[qp], label=f'qp #{qp+1}')
                    ax2.plot(times[:], spins[qp][:,0,i,j], '-o', color=colors[qp], label=f'qp #{qp+1}')
                    ax3.plot(times[:], coup[qp][:,0,i,j] - spins[qp][:,0,i,j], '-o', color=colors[qp], label=f'qp #{qp+1}')
            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax1.set_ylabel(r'$l_{' + str(i+1) + str(j+1) + '}$', fontsize=14)
            ax2.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax2.set_ylabel(r'$\omega_{' + str(i+1) + str(j+1) + '}$', fontsize=14)
            ax3.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax3.set_ylabel(r'$l_{' + str(i+1) + str(j+1) + '}$ - $\omega_{' + str(i+1) + str(j+1) + '}$', fontsize=14)
            
    ax1, ax2, ax3 = axes1[1][2], axes2[1][2], axes3[1][2]
    if average == False:
        ax1.legend(bbox_to_anchor=(1.2,0.9))
        ax2.legend(bbox_to_anchor=(1.2,0.9))
        ax3.legend(bbox_to_anchor=(1.2,0.9))
        
    plt.figure(name1)
    plt.tight_layout()
    plt.figure(name2)
    plt.tight_layout()
    plt.figure(name3)
    plt.tight_layout()

    if plot_body_couples:
        fig1.savefig(plot_body_couples)
    if plot_spin_inertias:
        fig2.savefig(plot_spin_inertias)
    if plot_spin_diff:
        fig3.savefig(plot_spin_diff)
    
    return 0


def deviatoric(stress):
    '''Calculate the deviatoric component of a stress quantity

    :param array-like stress: A second order tensor or 9-component slice of a third order tensor

    :returns: deviatoric of ``stress``
    '''

    return(stress - (1/3)*numpy.trace(stress))


def plot_stress_norm(cauchy, symm, m, nqp, nel, ninc, times, output_name):
    '''Plot the infinity norms of deviatoric Cauchy, symmetric micro, and higher stresses

    :param dict cauchy: The quantities dict storing Cauchy stress
    :param dict symm: The quantities dict storing symmetric micro stress
    :param dict m: THe quantities dict storing higher order stress
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param array-like times: The time increments
    :param str output_name: Output filename

    :returns: ``output_name`` plot
    '''

    fig = plt.figure('stress_norms', figsize=(9,9), dpi=300)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)

    for el in range(nel):
        for qp in range(nqp):
            # take norms
            cauchy_norm, symm_norm, diff_norm = [], [], []
            m1_norm, m2_norm, m3_norm = [], [], []
            for t in range(ninc):
                cauchy_norm.append(numpy.linalg.norm(deviatoric(cauchy[qp][t,el,:,:]), ord='fro'))
                symm_norm.append(numpy.linalg.norm(deviatoric(symm[qp][t,el,:,:]), ord='fro'))
                diff_norm.append(numpy.linalg.norm(deviatoric(cauchy[qp][t,el,:,:]-symm[qp][t,el,:,:]), ord='fro'))
                m1_norm.append(numpy.linalg.norm(deviatoric(m[qp][t,el,:,:,0]), ord='fro'))
                m2_norm.append(numpy.linalg.norm(deviatoric(m[qp][t,el,:,:,1]), ord='fro'))
                m3_norm.append(numpy.linalg.norm(deviatoric(m[qp][t,el,:,:,2]), ord='fro'))
            label = f"qp  #{(qp+1)+(nel*8)}"
            ax1.plot(times, cauchy_norm, label=label)
            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax1.set_ylabel(r'$||dev\left(\sigma_{ij}\right)||$', fontsize=14)
            ax2.plot(times, symm_norm, label=label)
            ax2.set_ylabel(r'$||dev\left(s_{ij}\right)||$', fontsize=14)
            ax2.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax3.plot(times, diff_norm, label=label)
            ax3.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax3.set_ylabel(r'$||dev\left(\sigma_{ij} - s_{ij}\right)||$', fontsize=14)
            ax4.plot(times, m1_norm, '-o', label=f"index k = 1, {label}")
            ax4.plot(times, m2_norm, '-^', label=f"index k = 2, {label}")
            ax4.plot(times, m3_norm, '-v', label=f"index k = 3, {label}")
            ax4.set_xlabel(r'Simulation time (s)', fontsize=14)
            ax4.set_ylabel(r'$||dev\left(m_{ijk}\right)||$', fontsize=14)
            # ax4.legend([(label='index 1', color='b'),
                        # (label='index 2', color='b'),
                        # (label='index 3', color='b')])
            ax4.legend(['index k = 1', 'index k = 2', 'index k = 3'])
            ax1.set_xlabel(r'Simulation time (s)', fontsize=14)

    fig.savefig(output_name)

    return 0


def csv_machine(quantity, output_file, nqp, nel, time, three_point=False):
    '''Write a text dump summarizing a variety of statistics

    :param dict quantities: A 2nd or third order tensor dictionary with keys for quadrature points and values storing a 4d array where indices correspond to time, element number, component i, and component j
    :param str output_file: Output filename
    :param int nqp: The dictionary key of the quantityt to output
    :param int nel: The desired element of the quantity to output
    :param int time: The desired time increment to output the quantity
    :param bool three_point: Parameter indicating if tensor is 2nd or 3rd order. "False" by default

    :returns: ``output_file``
    '''

    comps, means, mins, maxs, devs, covs = [], [], [], [], [], []
    if three_point == False:
        for i in range(3):
            for j in range(3):
                comps.append(f"{i+1}{j+1}")
                everything = numpy.array([quantity[qp][time,el,i,j] for qp in range(nqp) for el in range(nel)]).flatten()
                means.append(numpy.mean(everything))
                mins.append(numpy.min(everything))
                maxs.append(numpy.max(everything))
                devs.append(numpy.std(everything))
                covs.append(numpy.std(everything)/numpy.mean(everything))
    else:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    comps.append(f"{i+1}{j+1}{k+1}")
                    everything = numpy.array([quantity[qp][time,el,i,j,k] for qp in range(nqp) for el in range(nel)]).flatten()
                    means.append(numpy.mean(everything))
                    mins.append(numpy.min(everything))
                    maxs.append(numpy.max(everything))
                    devs.append(numpy.std(everything))
                    covs.append(numpy.std(everything)/numpy.mean(everything))

    df = pandas.DataFrame({'component': comps,
                           'mean': means,
                           'min': mins,
                           'max': maxs,
                           'dev': devs,
                           'cov': covs})
    print(f'building csv file for {output_file}')
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0

def dump_all(quantity, output_file, icomp, jcomp, nqp, nel, time):
    '''Dump all values for particular components of a 2nd order tensor for a given time value and element number

    :param dict quantities: A 2nd order tensor dictionary with keys for quadrature points and values storing a 4d array where indices correspond to time, element number, component i, and component j
    :param str output_file: Output filename
    :param int icomp: the desired i component of the quantity to output
    :param int jcomp: the desired j component of the quantity to output
    :param int nqp: The number of quadrature points
    :param int nel: The desired element of the quantity to output
    :param int time: The desired time increment to output the quantity

    :returns: ``output_file``
    '''

    everything = numpy.array([quantity[qp][time,el,icomp,jcomp] for qp in range(nqp) for el in range(nel)]).flatten()
    df = pandas.DataFrame({f'quantity': everything})
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def average_quantities(quantities):
    '''Average 2nd order tensor quantites over 8 quadrature points

    :param dict quantities: A 2nd order tensor dictionary with keys for quadrature points and values storing a 4d array where indices correspond to time, element number, component i, and component j

    :returns: ``output`` dict with same indices as ``quantities`` and a single key
    '''

    output = {}
    output[0] = numpy.zeros_like(quantities[0])
    for i in range(3):
        for j in range(3):
            mean_field = []
            for qp in quantities.keys():
                mean_field.append(quantities[qp][:,0,i,j])
            means = numpy.mean(mean_field, axis=0)
            output[0][:,0,i,j] = means

    return(output)


def csv_all_quantities(output_file, times, nqp, quantities, quantity_names):
    '''Create csv file of a list of quantities for a single domain

    :param str output_file: Output filename
    :param array-like times: Array of times
    :param int nqp: The number of quadrature points
    :param list quantities: A list containing Micromorphic Filter output quantities
    :param list quantity_names: A list containing the names of Micromorphic Filter output quantities

    :returns: Write ``output_file``
    '''

    # Collect all results
    output_dict = {}
    output_dict['time'] = numpy.array([[time for qp in range(nqp) for time in times]]).flatten()
    output_dict['qp'] = numpy.array([[qp for qp in range(nqp) for time in times]]).flatten()
    for quantity, name in zip(quantities, quantity_names):
        # volume fractions
        if len(numpy.shape(quantity[0])) == 0:
            key = f'{name}'
            qp_data = []
            for qp in range(nqp):
                for t in range(len(times)):
                    qp_data.append(quantity[qp])
            output_dict[key] = numpy.array(qp_data).flatten()
        # vectors
        elif len(numpy.shape(quantity[0])) == 3:
            for i in range(3):
                key = f'{name}_{i+1}'
                qp_data = []
                for qp in range(nqp):
                    qp_data.append(quantity[qp][:,0,i])
                output_dict[key] = numpy.array(qp_data).flatten()
        # second order tensors
        elif len(numpy.shape(quantity[0])) == 4:
            for i in range(3):
                for j in range(3):
                    key = f'{name}_{i+1}{j+1}'
                    qp_data = []
                    for qp in range(nqp):
                        qp_data.append(quantity[qp][:,0,i,j])
                    output_dict[key] = numpy.array(qp_data).flatten()
        # third order tensors
        elif len(numpy.shape(quantity[0])) == 5:
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        key = f'{name}_{i+1}{j+1}{k+1}'
                        qp_data = []
                        for qp in range(nqp):
                            qp_data.append(quantity[qp][:,0,i,j,k])
                        output_dict[key] = numpy.array(qp_data).flatten()
        else:
            print(f'Data shape not supported! quantity = {quantity}')

    # Output CSV
    df = pandas.DataFrame(output_dict)
    print(f'building csv file for {output_file}')
    df.to_csv(output_file, header=True, sep=',', index=False)

    return 0


def calculate_initial_grain_volume_fraction(density, rho_binder, rho_grain, nqp):
    '''Calculate the initial grain volume fraction given the homogenized density and two reference densities

    :param dict density: Dictionary containing homogenized densities
    :param str rho_binder: The density of the binder material, required if 'csv-all-quantities-single-domain' is specified
    :param str rho_grain: The density of the grain material, required if 'csv-all-quantities-single-domain' is specified
    :param int nqp: The number of quadrature points

    :returns: Dictionary of grain volume fractions for each quadrature point
    '''

    fgs = {}
    for qp in range(nqp):
        fgs[qp] = (density[qp][0]-rho_binder)/(rho_grain-rho_binder)

    return fgs


def visualize_results(input_file, average, num_domains,
                      plot_cauchy_couple=None,
                      plot_cauchy_stress=None,
                      plot_PK2_stress=None,
                      plot_symm_stress=None,
                      plot_SIGMA_stress=None,
                      plot_stress_diff=None,
                      plot_body_couples=None,
                      plot_spin_inertias=None,
                      plot_spin_diff=None,
                      plot_rotation_diff=None,
                      plot_stretch_diff=None,
                      plot_stress_norms=None,
                      csv_cauchy=None,
                      csv_PK2=None,
                      csv_GLstrain=None,
                      csv_estrain=None,
                      csv_ref_mod=None,
                      csv_cur_mod=None,
                      csv_symm=None,
                      csv_stress_diff=None,
                      csv_m=None,
                      csv_M=None,
                      csv_stress33_all=None,
                      csv_all_quantities_single_domain=None,
                      rho_binder=None,
                      rho_grain=None):
    ''' Post-process Micromorphic Filter output

    :param str input_file: The XDMF Micromorphic Filter results file
    :param bool average: Average over quadrature points if True
    :param int num_domains: The number of filter domains
    :param str plot_cauchy_couple: Optional filename to plot Cauchy couple vs. simulation time
    :param str plot_cauchy_stress: Optional filename to plot Cauchy stress vs. Eulerian strain
    :param str plot_PK2_stress: Optional filename to plot PK2 stress vs. Green-Lagrange strain
    :param str plot_symm_stress: Optional filename to plot symmetric micro stress vs. Eulerian strain
    :param str plot_SIGMA_stress: Optional filename to plot Symmetric micro stress vs. Green-Lagrange strain
    :param str plot_stress_diff: Optional filename to plot difference between Cauchy and symmetric micro stresses vs. simulation time
    :param str plot_body_couples: Optional filename to plot body couples vs. simulation time
    :param str plot_spin_inertias: Optional filename to plot micro spin inertias vs. simulation time
    :param str plot_spin_diff: Optional filename to plot difference between body couples and micro spin inertias vs. simulation time
    :param str plot_rotation_diff: Optional filename to plot difference between macro and micro rotations vs. simulation time
    :param str plot_stress_diff: Optional filename to plot differences between macro and micro stretches vs. simulation time
    :param str plot_stress_norms: Optional filename to plot norms of cauchy stress, symmetric micro stress, difference between Cauchy and symmetric micro stresses, and higher order stress
    :param str csv_cauchy: Optional filename for csv output of Cauchy stress summary statistics
    :param str csv_PK2: Optional filename for csv output of PK2 stress summary statistics
    :param str csv_GLstrain: Optional filename for csv output of Green-Lagrange strain summary statistics
    :param str csv_ref_mod: Optional filename for csv output of 'moduli' calculation (S_{ij} / E_{ij}) in reference configuration summary statistics
    :param str csv_cur_mod: Optional filename for csv output of 'moduli' calculation (\sigma_{ij} / e_{ij}) in the current configuration summary statistics
    :param str csv_strain: Optional filename for csv output of Eulerian strain summary statistics
    :param str csv_symm: Optional filename for csv output of symmetric micro stress summary statistics
    :param str csv_stress_diff: Optional filename for csv output of difference between Cauchy and symmetric micro stresses summary statistics
    :param str csv_m: Optional filename for csv output of couple stress (current configuration) summary statistics
    :param str csv_M: Optional filename for csv output of couple stress (reference configuration) summary statistics
    :param str csv_stress33_all: Optional filename for csv output of all Cauchy 33 values
    :param str csv_all_quantities_single_domain: Optional filename for csv output of all quantities for a single domain
    :param str rho_binder: The density of the binder material, required if 'csv-all-quantities-single-domain' is specified
    :param str rho_grain: The density of the grain material, required if 'csv-all-quantities-single-domain' is specified
    '''

    # Read in data
    data, geometry, topology = XRT.parse_xdmf_output(input_file)
    nqp = 8
    ninc = numpy.shape(data['time'])[0]
    nel = num_domains

    # Read in the DOF information
    displacement, gradu, phi, gradphi = XRT.construct_degrees_of_freedom(data, nqp, nel)

    # Read in the stress information
    cauchy, symm, m = XRT.collect_stresses(data, nqp, nel)
    PK2, SIGMA, M = XRT.get_reference_configuration_stresses(data, nqp, nel)

    # Read in the strain information
    E, Ecal, Gamma, F, chi, grad_chi, e, h = XRT.compute_deformations(data, nqp, nel)

    # Read in first moment of momentum measures
    coup, spins = XRT.collect_first_moment_of_momentum_measures(data, nqp, nel)

    # Get times
    times = numpy.unique(data['time'])
    
    # average fields if selected as an option
    if average:
        cauchy = average_quantities(cauchy)
        symm   = average_quantities(symm)
        E      = average_quantities(E)
        PK2    = average_quantities(PK2)
        SIGMA  = average_quantities(SIGMA)
        F      = average_quantities(F)
        chi    = average_quantities(chi)
        e      = average_quantities(e)
        h      = average_quantities(h)
        coup   = average_quantities(coup)
        spins  = average_quantities(spins)

    # Plot stresses
    ## plot cauchy and PK2 stress
    if plot_cauchy_stress:
        plot_current_stresses(e, cauchy, average, output_name=plot_cauchy_stress)
    if plot_PK2_stress:
        plot_reference_stresses(E, PK2, average, output_name=plot_PK2_stress)
    ## plot symmetric micro-stress, s and SIGMA
    if plot_symm_stress:
        plot_current_stresses(e, symm, average, output_name=plot_symm_stress)
    if plot_SIGMA_stress:
        plot_reference_stresses(E, SIGMA, average, output_name=plot_SIGMA_stress)

    ## Check #1: Symmetry of Cauchy Stress using Axial Vector
    if plot_cauchy_couple:
        plot_axial_vector_mag(cauchy, times, average, output_name=plot_cauchy_couple)

    ## Check #2: Difference between Cauchy and Symmetric Micro-Stresses
    if plot_stress_diff:
        plot_stress_diffs(cauchy, symm, times, average, output_name=plot_stress_diff)

    ## Check #3: Micromorphic Effects - Difference between macro- and micro-deformation gradients
    # Get R and U
    if plot_rotation_diff or plot_stretch_diff:
        R, U, Rchi, Uchi = get_R_and_U(PK2, F, chi, times)
    # plot difference between rotations
    if plot_rotation_diff:
        plot_rot_diffs(PK2, times, R, Rchi, average, output_name=plot_rotation_diff)
    # plot difference between stretches
    if plot_stretch_diff:
        plot_stretch_diffs(PK2, times, R, Rchi, average, output_name=plot_stress_diff)
    
    # Check #4: Different between body_couple and micro_spin_inertias, plot them too
    if plot_body_couples or plot_spin_inertias or plot_spin_diff:
        plot_first_moment_of_momentum_measures(coup, spins, times, average, plot_body_couples, plot_spin_inertias, plot_spin_diff)

    # Plot stress norms
    if plot_stress_norms:
        plot_stress_norm(cauchy, symm, m, nqp, nel, ninc, times, plot_stress_norms)

    # Output csvs
    #TODO: add an argument for specifying what time to gather statistics on
    if csv_cauchy:
        csv_machine(cauchy, csv_cauchy, nqp, nel, ninc-1)
    if csv_PK2:
        csv_machine(PK2, csv_PK2, nqp, nel, ninc-1)
    if csv_GLstrain:
        csv_machine(E, csv_GLstrain, nqp, nel, ninc-1)
    if csv_estrain:
        csv_machine(e, csv_estrain, nqp, nel, ninc-1)
    if csv_ref_mod:
        mod = {}
        for qp in range(nqp):
            values = numpy.zeros((ninc, nel, 3, 3))
            for el in range(nel):
                for t in range(ninc):
                    values[t,el,:,:] = PK2[qp][t,el,:,:] / E[qp][t,el,:,:]
                mod.update({qp:numpy.copy(values)})
        csv_machine(mod, csv_ref_mod, nqp, nel, ninc-1)
    if csv_cur_mod:
        mod = {}
        for qp in range(nqp):
            values = numpy.zeros((ninc, nel, 3, 3))
            for el in range(nel):
                for t in range(ninc):
                    values[t,el,:,:] = cauchy[qp][t,el,:,:] / e[qp][t,el,:,:]
                mod.update({qp:numpy.copy(values)})
        csv_machine(mod, csv_cur_mod, nqp, nel, ninc-1)
    if csv_symm:
        csv_machine(symm, csv_symm, nqp, nel, ninc-1)
    # TODO: make this better lol
    if csv_stress_diff:
        diff = {}
        for qp in range(nqp):
            values = numpy.zeros((ninc, nel, 3, 3))
            for el in range(nel):
                for t in range(ninc):
                    values[t,el,:,:] = cauchy[qp][t,el,:,:] - symm[qp][t,el,:,:]
                diff.update({qp:numpy.copy(values)})
        csv_machine(diff, csv_stress_diff, nqp, nel, ninc-1)
    # Higher order stresses
    if csv_m:
        csv_machine(m, csv_m, nqp, nel, ninc-1, three_point=True)
    if csv_M:
        csv_machine(M, csv_M, nqp, nel, ninc-1, three_point=True)

    # Dump all cauchy
    if csv_stress33_all:
        dump_all(cauchy, csv_stress33_all, 2, 2, nqp, nel, ninc-1)

    # Summary csv for quantities of a single single domains
    if csv_all_quantities_single_domain:
        # get initial volume fractions
        assert rho_binder
        assert rho_grain
        fgs = calculate_initial_grain_volume_fraction(data['density_0'], rho_binder, rho_grain, nqp)
        print(f'fgs = {fgs}')
        quantities = [fgs, displacement, gradu, phi, gradphi, F, chi, E, Ecal, Gamma, e, h, cauchy, symm, m]
        quantities_names = ['f_g', 'u', 'gradu', 'phi', 'gradphi', 'F', 'chi', 'E', 'Ecal', 'Gamma', 'e', 'h', 'cauchy', 'symm', 'm']
        csv_all_quantities(csv_all_quantities_single_domain, times, nqp, quantities, quantities_names)

    return 0


def get_parser():

    filename = inspect.getfile(lambda: None)
    basename = os.path.basename(filename)
    basename_without_extension, extension = os.path.splitext(basename)
    cli_description = "Post-process Micromorphic Filter Output"
    parser = argparse.ArgumentParser(description=cli_description,
                                     prog=os.path.basename(filename))
    parser.add_argument('-i', '--input-file', type=str, required=True,
        help="The XDMF Micromorphic Filter results file")
    parser.add_argument('--average', type=str, required=False, default=False,
        help='Boolean whether or not homogenized DNS results will be averaged')
    parser.add_argument('--num-domains', type=int, required=False, default=1,
        help="Specify the number of filter domains")
    parser.add_argument('--plot-cauchy-couple', type=str, required=False, default=None,
        help="Optional filename to plot Cauchy couple vs. simulation time")
    parser.add_argument('--plot-cauchy-stress', type=str, required=False, default=None,
        help="Optional filename to plot Cauchy stress vs. Eulerian strain")
    parser.add_argument('--plot-PK2-stress', type=str, required=False, default=None,
        help="Optional filename to plot PK2 stress vs. Green-Lagrange strain")
    parser.add_argument('--plot-symm-stress', type=str, required=False, default=None,
        help="Optional filename to plot symmetric micro stress vs. Eulerian strain")
    parser.add_argument('--plot-SIGMA-stress', type=str, required=False, default=None,
        help="Optional filename to plot Symmetric micro stress vs. Green-Lagrange strain")
    parser.add_argument('--plot-stress-diff', type=str, required=False, default=None,
        help="Optional filename to plot difference between Cauchy and symmetric micro\
              stresses vs. simulation time")
    parser.add_argument('--plot-body-couples', type=str, required=False, default=None,
        help="Optional filename to plot body couples vs. simulation time")
    parser.add_argument('--plot-spin-inertias', type=str, required=False, default=None,
        help="Optional filename to plot micro spin inertias vs. simulation time")
    parser.add_argument('--plot-spin-diff', type=str, required=False, default=None,
        help="Optional filename to plot difference between body couples and micro spin\
              inertias vs. simulation time")
    parser.add_argument('--plot-rotation-diff', type=str, required=False, default=None,
        help="Optional filename to plot difference between macro and micro rotations\
              vs. simulation time")
    parser.add_argument('--plot-stretch-diff', type=str, required=False, default=None,
        help="Optional filename to plot differences between macro and micro stretches\
              vs. simulation time")
    parser.add_argument('--plot-stress-norms', type=str, required=False, default=None,
        help="Optional filename to plot norms of cauchy stress, symmetric micro stress,\
              difference between Cauchy and symmetric micro stresses, and higher order\
              stress.")
    parser.add_argument('--csv-cauchy', type=str, required=False, default=None,
        help="Optional filename for csv output of Cauchy stress summary statistics")
    parser.add_argument('--csv-PK2', type=str, required=False, default=None,
        help="Optional filename for csv output of PK2 stress summary statistics")
    parser.add_argument('--csv-GLstrain', type=str, required=False, default=None,
        help="Optional filename for csv output of Green-Lagrange strain summary statistics")
    parser.add_argument('--csv-ref-mod', type=str, required=False, default=None,
        help="Optional filename for csv output of 'moduli' calculation (S_{ij} / E_{ij})\
              in reference configuration summary statistics")
    parser.add_argument('--csv-cur-mod', type=str, required=False, default=None,
        help="Optional filename for csv output of 'moduli' calculation (\sigma_{ij} / e_{ij})\
              in the current configuration summary statistics")
    parser.add_argument('--csv-estrain', type=str, required=False, default=None,
        help="Optional filename for csv output of Eulerian strain summary statistics")
    parser.add_argument('--csv-symm', type=str, required=False, default=None,
        help="Optional filename for csv output of symmetric micro stress summary statistics")
    parser.add_argument('--csv-stress-diff', type=str, required=False, default=None,
        help="Optional filename for csv output of difference between Cauchy and symmetric\
              micro stresses summary statistics")
    parser.add_argument('--csv-m', type=str, required=False, default=None,
        help="Optional filename for csv output of couple stress (current configuration)\
              summary statistics")
    parser.add_argument('--csv-M', type=str, required=False, default=None,
        help="Optional filename for csv output of couple stress (reference configuration)\
              summary statistics")
    parser.add_argument('--csv-stress33-all', type=str, required=False, default=None,
        help="Optional filename for csv output of all Cauchy 33 values")
    parser.add_argument('--csv-all-quantities-single-domain', type=str, required=False, default=None,
        help="Optional filename for csv output of all quantities for a single domain")
    parser.add_argument('--rho-binder', type=float, required=False, default=None,
        help="The density of the binder material, required if '--csv-all-quantities-single-domain' is specified")
    parser.add_argument('--rho-grain', type=float, required=False, default=None,
        help="The density of the grain material, required if '--csv-all-quantities-single-domain' is specified")

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(visualize_results(
                input_file=args.input_file,
                average=str2bool(args.average),
                num_domains=args.num_domains,
                plot_cauchy_couple=args.plot_cauchy_couple,
                plot_cauchy_stress=args.plot_cauchy_stress,
                plot_PK2_stress=args.plot_PK2_stress,
                plot_symm_stress=args.plot_symm_stress,
                plot_SIGMA_stress=args.plot_SIGMA_stress,
                plot_stress_diff=args.plot_stress_diff,
                plot_body_couples=args.plot_body_couples,
                plot_spin_inertias=args.plot_spin_inertias,
                plot_spin_diff=args.plot_spin_diff,
                plot_rotation_diff=args.plot_rotation_diff,
                plot_stretch_diff=args.plot_stretch_diff,
                plot_stress_norms=args.plot_stress_norms,
                csv_cauchy=args.csv_cauchy,
                csv_PK2=args.csv_PK2,
                csv_GLstrain=args.csv_GLstrain,
                csv_estrain=args.csv_estrain,
                csv_ref_mod=args.csv_ref_mod,
                csv_cur_mod=args.csv_cur_mod,
                csv_symm=args.csv_symm,
                csv_stress_diff=args.csv_stress_diff,
                csv_m=args.csv_m,
                csv_M=args.csv_M,
                csv_stress33_all=args.csv_stress33_all,
                csv_all_quantities_single_domain=args.csv_all_quantities_single_domain,
                rho_binder=args.rho_binder,
                rho_grain=args.rho_grain,
                ))