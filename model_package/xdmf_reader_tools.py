import sys
import inspect
import os
import argparse
import numpy as np
import h5py
import xml.etree.ElementTree as ET

from scipy.linalg import norm
from scipy.linalg import polar
import argparse


def get_attribute_data(xml_node, path='.'):
    '''Collect attribute data from xml formatted file

    :param node xml_node: The xml child node containing path to attribute data
    :param str path: a path for locating data within an HDF5 file, default='.'

    :returns: value and keyname of attribute data
    '''

    if ("Attribute" != xml_node.tag):
        raise ValueError("XML node is not of type 'Attribute'")
    value = None
    for child in xml_node:
        if "DataItem" == child.tag:
            value, keyname = read_data(child, path=path)
            break

    return value, keyname


def get_set(xml_node, path='.'):
    '''Collect set data from xml formatted file

    :param node xml_node: The xml child node containing path to attribute data
    :param str path: a path for locating data within an HDF5 file, default='.'

    :returns: value of set data
    '''

    if ("Set" != xml_node.tag):
        raise ValueError("XML node is not of type 'Set'")
    for child in xml_node:
        if "DataItem" == child.tag:
            value = read_data(child, path=path)
            break

    return value


def get_geometry(xml_node, path='.'):
    '''Collect geometry data from xml formatted file

    :param node xml_node: The xml child node containing path to geometry data
    :param str path: a path for locating data within an HDF5 file, default='.'

    :returns: value of geometry data
    '''

    if ("Geometry" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Geometry'")
    for child in xml_node:
        if "DataItem" == child.tag:
            value = read_data(child, path=path)
            break

    return value


def get_topology(xml_node, path='.'):
    '''Collect topology data from xml formatted file

    :param node xml_node: The xml child node containing path to topology data
    :param str path: a path for locating data within an HDF5 file, default='.'

    :returns: value of topology data
    '''

    if ("Topology" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Topology'")
    for child in xml_node:
        if "DataItem" == child.tag:
            value = read_data(child, path=path)
            break

    return value


def read_data(xml_node, path='.'):
    '''Collect data from HDF5 file using path specified in xml formatted XDMF file

    :param node xml_node: The xml child node containing path to data
    :param str path: a path for locating data within an HDF5 file, default='.'

    :returns: value of topology data
    '''

    if ("DataItem" != xml_node.tag):
        raise ValueError("XML node is not of type 'DataItem'")
    shape = [int(v) for v in xml_node.attrib["Dimensions"].split()]

    if xml_node.attrib["Format"] == "XML":
        value = np.hstack([np.array(line.split()).astype(dtype) for line in xml_node.text.split("\n")]).reshape(shape)
    else:
        filename, keyname = xml_node.text.split(':')
        with h5py.File(os.path.join(path, filename), 'r') as h5file:
            value = np.array(h5file[keyname]).reshape(shape)

    return value, keyname


def parse_xdmf_output(input_file):
    '''Parse XDMF and HDF5 file contents into attributes, geometry, topology, and time

    :param str input_file: The XDMF filename

    :returns: dictionaries for simulation data, geometry, and topology
    '''

    tree = ET.parse(input_file)
    root = tree.getroot()

    data, geometry, topology = {}, {}, {}

    for domain in root:
        for collection in domain:
            if collection.tag != "Grid":
                continue

            for grid in collection:
                if grid.tag != "Grid":
                    continue

                for child in grid:
                    # Attribute
                    if ("Attribute" == child.tag):
                        name = child.attrib["Name"]
                        value, keyname = get_attribute_data(child)

                        if keyname in data.keys():
                            data[keyname].append(value)
                        else:
                            data.update({keyname:[value]})

                    # Geometry
                    if ("Geometry" == child.tag):
                        name = "geometry"
                        value = get_geometry(child)
                        if name in geometry.keys():
                            geometry[name].append(value)
                        else:
                            geometry.update({name:[value]})

                    #Topology
                    if ("Topology" == child.tag):
                        name = "topology"
                        value = get_topology(child)
                        if name in topology.keys():
                            topology[name].append(value)
                        else:
                            topology.update({name:[value]})

                    # Time
                    if ("Time" == child.tag):
                        name = "time"
                        value = float(child.attrib['Value'])
                        if name in data.keys():
                            data[name].append(value)
                        else:
                            data.update({name:[value]})

    return(dict([(name, np.vstack(data[name])) for name in data.keys()]), geometry, topology)


def construct_degrees_of_freedom(data, nqp, nel, dim = 3):
    '''Collect quadrature point data for displacement, displacement gradient, micro-displacement, and micro-displacement gradient from Micromorphic Filter output

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for displacement, displacement gradient, micro-displacement, and micro-displacement gradient
    '''

    ninc = np.shape(data['time'])[0]

    displacement, grad_u, phi, grad_phi = {}, {}, {}, {}
    for qp in range(nqp):
        #for el in range(nel):
            # collect the displacement
            values = np.zeros((ninc, nel, dim))
            for t in range(ninc):
                root_string = f'displacement_qpt_{qp}_{t}'
                values[t,:,:] = data[root_string]
            displacement.update({qp:np.copy(values)})

            # collect the gradient of the macro deformation
            values = np.zeros((ninc, nel, dim, dim))
            for t in range(ninc):
                root_string = f'current_displacement_gradient_qpt_{qp}_{t}'
                values[t,:,:,:] = data[root_string].reshape((nel, dim,dim))
            grad_u.update({qp:np.copy(values)})

            # collect the micro deformation
            values = np.zeros((ninc, nel, dim, dim))
            for t in range(ninc):
                root_string = f'micro_displacement_qpt_{qp}_{t}'
                values[t,:,:,:] = data[root_string].reshape((nel, dim,dim))
            phi.update({qp:np.copy(values)})

            # collect the gradient of the micro deformation
            values = np.zeros((ninc, nel, dim, dim, dim))
            for t in range(ninc):
                root_string = f'current_micro_displacement_gradient_qpt_{qp}_{t}'
                values[t,:,:,:,:] = data[root_string].reshape((nel,dim,dim,dim))
            grad_phi.update({qp:np.copy(values)})

    return(displacement, grad_u, phi, grad_phi)


def collect_stresses(data, nqp, nel, dim=3):
    '''Collect quadrature point data for Cauchy, symmetric micro-, and higher order stresses from Micromorphic Filter output (all in current configuration)

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for Cauchy, symmetric micro-, and higher order stresses
    '''

    ninc = np.shape(data['time'])[0]

    cauchy_stress, symmetric_micro_stress, higher_order_stress = {}, {}, {}

    for qp in range(nqp):

        # Collect the Cauchy Stress
        values = np.zeros((ninc, nel, dim, dim))
        for t in range(ninc):
            root_string = f'cauchy_stress_{qp}_{t}'
            values[t,:,:,:] = data[root_string].reshape((nel,dim,dim))
        cauchy_stress.update({qp:np.copy(values)})

        # collect the symmetric micro stress
        values = np.zeros((ninc, nel, dim, dim))
        for t in range(ninc):
            root_string = f'symmetric_micro_stress_{qp}_{t}'
            values[t,:,:,:] = data[root_string].reshape((nel,dim,dim))
        symmetric_micro_stress.update({qp:np.copy(values)})

        # collect the higher order stress
        values = np.zeros((ninc, nel, dim, dim, dim))
        for t in range(ninc):
            root_string = f'higher_order_stress_{qp}_{t}'
            values[t,:,:,:,:] = data[root_string].reshape((nel,dim,dim,dim))
        higher_order_stress.update({qp:np.copy(values)})

    return(cauchy_stress, symmetric_micro_stress, higher_order_stress)


def get_reference_configuration_stresses(data, nqp, nel, dim=3):
    '''Map Cauchy, symmetric micro-, and higher order stresses to the reference configuraiton

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for Second Piola Kirchhoff, Symmetric micro-, and Higher order stresses (all in reference configuration)
    '''

    position, grad_u, phi, grad_phi = construct_degrees_of_freedom(data, nqp, nel, dim=dim)
    cauchy_stress, symmetric_micro_stress, higher_order_stress = collect_stresses(data, nqp, nel, dim=dim)
    PK2, SIGMA, M = {}, {}, {}

    ninc = np.shape(data['time'])[0]

    for qp in range(nqp):
        PK2.update({qp:np.zeros((ninc, nel, dim, dim))})
        SIGMA.update({qp:np.zeros((ninc, nel, dim, dim))})
        M.update({qp:np.zeros((ninc, nel, dim, dim, dim))})

        for inc in range(ninc):
            for el in range(nel):
                # construct the deformation gradient
                F = grad_u[qp][inc,el,:,:] + np.eye(dim)
                Finv = np.linalg.inv(F)
                J = np.linalg.det(F)

                # construct the micro-displacement
                chi = phi[qp][inc,el,:,:] + np.eye(dim)
                chiinv = np.linalg.inv(chi)

                # pull back (right??)
                PK2[qp][inc,el,:,:] = J*np.einsum("Kk,kl,Ll->KL", Finv, cauchy_stress[qp][inc,el,:,:], Finv)
                SIGMA[qp][inc,el,:,:] = J*np.einsum("Kk,kl,Ll->KL", Finv, symmetric_micro_stress[qp][inc,el,:,:], Finv)
                M[qp][inc,el,:,:,:] = J*np.einsum("Kk,Ll,Mm,klm->KLM", Finv, Finv, chiinv, higher_order_stress[qp][inc,el,:,:,:])

    return(PK2, SIGMA, M)


def map_sim(stress, ninc, dim=3):
    '''Map a flattened 2nd order stress tensor to index component notation. This function is used for converting output from ``micromorphic.evaluate_model`` to a convenient form for post-processing against Micromorphic Filter output data.

    :param dict stress: The dictionary of flattened 2nd order stress tensor
    :param int ninc: The number of time increments
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionary with reshaped stress data
    '''

    # ignore m and M for now!!!
    nel = 1
    nqp = len([key for key in stress.keys()])
    new_stress = {}

    for qp in range(nqp):
        new_stress.update({qp:np.zeros((ninc, nel, dim, dim))})
        k = 0
        for i in range(dim):
            for j in range(dim):
                new_stress[qp][:,:,i,j] = stress[qp][:,:,k]
                k = k + 1

    return(new_stress)


def get_current_configuration_stresses(PK2, SIGMA, grad_u, phi, dim=3):
    '''Convert Second Piola Kirchhoff and Symmetric micro- stresses to the current configuration

    :param dict PK2: A dictionary containing Second Piola Kirchhoff stress data
    :param dict SIGMA: A dictionary containing Symmetric micro-stress data
    :param dict grad_u: A dictionary containing displacement gradient data
    :param dict phi: A dicionary containing micro displacement data
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for Cauchy and symmetric micro-stresses (all in current configuration)
    '''

    # ignore m and M for now!!!
    ninc = PK2[0].shape[0]
    nel = PK2[0].shape[1]
    nqp = len([key for key in PK2.keys()])
    cauchy, symm, m = {}, {}, {}

    for qp in range(nqp):
        cauchy.update({qp:np.zeros((ninc, nel, dim, dim))})
        symm.update({qp:np.zeros((ninc, nel, dim, dim))})
        m.update({qp:np.zeros((ninc, nel, dim, dim, dim))})

        for inc in range(ninc):
            for el in range(nel):
                # construct the deformation gradient
                F = grad_u[qp][inc,el,:,:] + np.eye(dim)
                Ft = np.transpose(F)
                J = np.linalg.det(F)

                # construct the micro-displacement
                chi = phi[qp][inc,el,:,:] + np.eye(dim)

                # push forward, for m use eq. 3.33 in Miller Thesis
                cauchy[qp][inc,el,:,:] = (1./J)*np.einsum("kK,KL,Ll->kl", F, PK2[qp][inc,el,:,:], Ft)
                symm[qp][inc,el,:,:] = (1./J)*np.einsum("kK,KL,Ll->kl", F, SIGMA[qp][inc,el,:,:], Ft)
                #m[qp][inc,el,:,:,:] = (1./J)*np.einsum("kK,lL,mM,KLM->klm", F, Ft, chi, M[qp][inc,el,:,:,:])

    return(cauchy, symm)


def compute_deformations(data, nqp, nel, dim=3):
    '''Compute quadrature point data for a variety of deformation measures

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for Green-Lagrange strain, Micro-Green-Lagrange strain, micro-deformation gradient, deformation gradient, micro-deformation tensor, gradient of micro-deformation tensor, Euler-Almansi strain, and Hencky strain
    '''

    position, grad_u, phi, grad_phi = construct_degrees_of_freedom(data, nqp, nel, dim=dim)
    GreenLagrangeStrain, MicroGreenLagrangeStrain, Gamma = {}, {}, {}
    all_F, all_chi, all_grad_chi = {}, {}, {}
    EulerAlmansiStrain, HenckyStrain = {}, {}

    ninc = np.shape(data['time'])[0]

    for qp in range(nqp):
        GreenLagrangeStrain.update({qp:np.zeros((ninc,nel,dim,dim))})
        MicroGreenLagrangeStrain.update({qp:np.zeros((ninc,nel,dim,dim))})
        Gamma.update({qp:np.zeros((ninc, nel, dim, dim, dim))})
        all_F.update({qp:np.zeros((ninc,nel,dim,dim))})
        all_chi.update({qp:np.zeros((ninc,nel,dim,dim))})
        all_grad_chi.update({qp:np.zeros((ninc, nel, dim, dim, dim))})
        EulerAlmansiStrain.update({qp:np.zeros((ninc,nel,dim,dim))})
        HenckyStrain.update({qp:np.zeros((ninc,nel,dim,dim))})

        for inc in range(ninc):
            for el in range(nel):
                # calculate some intermediate terms
                F = grad_u[qp][inc,el,:,:] + np.eye(dim)
                Finv = np.linalg.inv(F)
                FT = np.transpose(F)
                FTinv = np.transpose(Finv)
                b = np.dot(F, FT)
                binv = np.linalg.inv(b)
                chi = phi[qp][inc,el,:,:] + np.eye(dim)
                grad_chi = grad_phi[qp][inc,el,:,:,:]
                # Get Hencky strains through eigenvalue decomposition
                bw, bv = np.linalg.eig(b)
                lam = np.log(np.sqrt(bw))
                Hencky_p = np.diag(lam) # principal components of Hencky
                Hencky = np.dot(bv, np.dot(Hencky_p, np.transpose(bv))) # rotate
                # Store intermediate terms
                all_F[qp][inc,el,:,:] = F
                all_chi[qp][inc,el,:,:] = chi
                all_grad_chi[qp][inc,el,:,:,:] = grad_chi
                # Store deformation measures
                GreenLagrangeStrain[qp][inc,el,:,:] = 0.5*(np.dot(FT, F) - np.eye(dim))
                MicroGreenLagrangeStrain[qp][inc,el,:,:] = np.dot(FT, chi) - np.eye(dim)
                Gamma[qp][inc,el,:,:,:] = np.einsum('iI,iJK->IJK', F, grad_chi)
                EulerAlmansiStrain[qp][inc,el,:,:] =  0.5*(np.eye(dim) - binv)
                HenckyStrain[qp][inc,el,:,:] = Hencky
    return(GreenLagrangeStrain, MicroGreenLagrangeStrain, Gamma, all_F, all_chi, all_grad_chi, EulerAlmansiStrain, HenckyStrain)


def collect_first_moment_of_momentum_measures(data, nqp, nel, dim=3):
    '''Collect body couples and micro-spin inertias 

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: dictionaries for body couples and micro-spin inertias
    '''

    body_couples = {}
    micro_spin_inertias = {}

    ninc = np.shape(data['time'])[0]

    for qp in range(nqp):

        # Collect body couples
        values = np.zeros((ninc, nel, dim, dim))
        for t in range(ninc):
            root_string = f'body_force_couple_{t}'
            values[t,:,:,:] = data[root_string][qp].reshape((3,3))
        body_couples.update({qp:np.copy(values)})

        # Collect micro_spin_inertias
        values = np.zeros((ninc, nel, dim, dim))
        for t in range(ninc):
            root_string = f'micro_spin_inertia_{t}'
            values[t,:,:,:] = data[root_string][qp].reshape((3,3))
        micro_spin_inertias.update({qp:np.copy(values)})

    return(body_couples, micro_spin_inertias)


def get_R_and_U(data, F, chi, nqp, nel, dim=3):
    '''Calculate stretch and rotation tensors for macro deformation and micro deformation tensors using polar decomposition

    :param dict data: The data dictionary containing output from Micromorphic Filter
    :param dict F: A dictionary containing macro deformation gradient information
    :param dict chi: A dictionary containing micro deformation tensor informaiton
    :param int nqp: The number of quadrature points
    :param int nel: The number of elements
    :param int dim: The number of spatial dimensions, default=3

    :returns: R, U, Rchi, and Uchi
    '''

    ninc = np.shape(data['time'])[0]

    R, U, Rchi, Uchi = {}, {}, {}, {}

    for qp in range(nqp):

        R_temp, U_temp = [], []
        R_chi_temp, U_chi_temp = [], []

        # calculate polar decompositions
        for t in range(ninc):

            Rr, Ur = polar(F[qp][t,:,:,:], side='left')
            Rchir, Uchir = polar(chi[qp][t,:,:,:], side='left')

            # check orthogonality
            tol = 1.e-9
            if (norm((np.dot(Rr, np.transpose(Rr)) - np.eye(3)),ord=2)) > tol:
                print('Error!!! R is not orthogonal!')
            if (norm((np.dot(Rchir, np.transpose(Rchir)) - np.eye(3)),ord=2)) > tol:
                print('Error!!! Rchi is not orthogonal!')

            # append results
            R_temp.append(Rr)
            U_temp.append(Ur)
            R_chi_temp.append(Rchir)
            U_chi_temp.append(Uchir)

        # Collect R
        values = np.zeros((ninc, nel, dim, dim))
        for i in range(dim):
            for j in range(dim):
                for t in range(ninc):
                    values[t,:,i,j] = R_temp[t][i][j]
        R.update({qp:np.copy(values)})

        # Collect U
        values = np.zeros((ninc, nel, dim, dim))
        for i in range(dim):
            for j in range(dim):
                for t in range(ninc):
                    values[t,:,i,j] = U_temp[t][i][j]
        U.update({qp:np.copy(values)})

        # Collect Rchi
        values = np.zeros((ninc, nel, dim, dim))
        for i in range(dim):
            for j in range(dim):
                for t in range(ninc):
                    values[t,:,i,j] = R_chi_temp[t][i][j]
        Rchi.update({qp:np.copy(values)})

        # Collect Uchi
        values = np.zeros((ninc, nel, dim, dim))
        for i in range(dim):
            for j in range(dim):
                for t in range(ninc):
                    values[t,:,i,j] = U_chi_temp[t][i][j]
        Uchi.update({qp:np.copy(values)})

    return(R, U, Rchi, Uchi)