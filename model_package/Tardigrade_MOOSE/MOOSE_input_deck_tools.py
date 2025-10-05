#!python
import argparse
import os
import pathlib
import sys
import yaml


def write_variables(file):
    '''Write the variables block for displacements and micro-displacements

    :params file file: The file to write to
    '''

    file.write('[Variables]\n')
    file.write('  [./disp_x]\n')
    file.write('  [../]\n')
    file.write('  [./disp_y]\n')
    file.write('  [../]\n')
    file.write('  [./disp_z]\n')
    file.write('  [../]\n')
    file.write('  [./phi_xx]\n')
    file.write('  [../]\n')
    file.write('  [./phi_yy]\n')
    file.write('  [../]\n')
    file.write('  [./phi_zz]\n')
    file.write('  [../]\n')
    file.write('  [./phi_yz]\n')
    file.write('  [../]\n')
    file.write('  [./phi_xz]\n')
    file.write('  [../]\n')
    file.write('  [./phi_xy]\n')
    file.write('  [../]\n')
    file.write('  [./phi_zy]\n')
    file.write('  [../]\n')
    file.write('  [./phi_zx]\n')
    file.write('  [../]\n')
    file.write('  [./phi_yx]\n')
    file.write('  [../]\n')

    return 0


def write_kernels(file, phis='on', internal_force='InternalForce', internal_couple='InternalCouple'):
    '''Write the kernels block for coupled kinematic variables

    :params file file: The file to write to
    :param str phis: Either "on" to activate phi coupling kernels, or "off" to deactivate
    '''

    file.write('[Kernels]\n')
    file.write('  #Define the internal force balance equations\n')
    file.write('  [./force_1]\n')
    file.write(f'    type = {internal_force}\n')
    file.write('    component = 0\n')
    file.write('    dof_num   = 0\n')
    file.write('    variable  = disp_x\n')
    file.write('    save_in = force_x\n')
    file.write('\n')
    file.write('    #Coupled variables\n')
    file.write('    u1     = disp_x\n')
    file.write('    u2     = disp_y\n')
    file.write('    u3     = disp_z\n')
    file.write('    phi_11 = phi_xx\n')
    file.write('    phi_22 = phi_yy\n')
    file.write('    phi_33 = phi_zz\n')
    file.write('    phi_23 = phi_yz\n')
    file.write('    phi_13 = phi_xz\n')
    file.write('    phi_12 = phi_xy\n')
    file.write('    phi_32 = phi_zy\n')
    file.write('    phi_31 = phi_zx\n')
    file.write('    phi_21 = phi_yx\n')
    file.write('  [../]\n')
    file.write('  [./force_2]\n')
    file.write(f'    type = {internal_force}\n')
    file.write('    component = 1\n')
    file.write('    dof_num   = 1\n')
    file.write('    variable  = disp_y\n')
    file.write('    save_in = force_y\n')
    file.write('\n')
    file.write('    #Coupled variables\n')
    file.write('    u1     = disp_x\n')
    file.write('    u2     = disp_y\n')
    file.write('    u3     = disp_z\n')
    file.write('    phi_11 = phi_xx\n')
    file.write('    phi_22 = phi_yy\n')
    file.write('    phi_33 = phi_zz\n')
    file.write('    phi_23 = phi_yz\n')
    file.write('    phi_13 = phi_xz\n')
    file.write('    phi_12 = phi_xy\n')
    file.write('    phi_32 = phi_zy\n')
    file.write('    phi_31 = phi_zx\n')
    file.write('    phi_21 = phi_yx\n')
    file.write('  [../]\n')
    file.write('  [./force_3]\n')
    file.write(f'    type = {internal_force}\n')
    file.write('    component = 2\n')
    file.write('    dof_num   = 2\n')
    file.write('    variable  = disp_z\n')
    file.write('    save_in = force_z\n')
    file.write('\n')
    file.write('    #Coupled variables\n')
    file.write('    u1     = disp_x\n')
    file.write('    u2     = disp_y\n')
    file.write('    u3     = disp_z\n')
    file.write('    phi_11 = phi_xx\n')
    file.write('    phi_22 = phi_yy\n')
    file.write('    phi_33 = phi_zz\n')
    file.write('    phi_23 = phi_yz\n')
    file.write('    phi_13 = phi_xz\n')
    file.write('    phi_12 = phi_xy\n')
    file.write('    phi_32 = phi_zy\n')
    file.write('    phi_31 = phi_zx\n')
    file.write('    phi_21 = phi_yx\n')
    file.write('  [../]\n')
    if phis == 'on':
        file.write('  #Define the internal couple balance equations\n')
        file.write('  [./couple_11]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 0\n')
        file.write('    component_j = 0\n')
        file.write('    dof_num     = 3\n')
        file.write('    variable    = phi_xx\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_12]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 0\n')
        file.write('    component_j = 1\n')
        file.write('    dof_num     = 4\n')
        file.write('    variable    = phi_xy\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_13]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 0\n')
        file.write('    component_j = 2\n')
        file.write('    dof_num     = 5\n')
        file.write('    variable    = phi_xz\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_21]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 1\n')
        file.write('    component_j = 0\n')
        file.write('    dof_num     = 6\n')
        file.write('    variable    = phi_yx\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_22]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 1\n')
        file.write('    component_j = 1\n')
        file.write('    dof_num     = 7\n')
        file.write('    variable    = phi_yy\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_23]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 1\n')
        file.write('    component_j = 2\n')
        file.write('    dof_num     = 8\n')
        file.write('    variable    = phi_yz\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_31]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 2\n')
        file.write('    component_j = 0\n')
        file.write('    dof_num     = 9\n')
        file.write('    variable    = phi_zx\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_32]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 2\n')
        file.write('    component_j = 1\n')
        file.write('    dof_num     = 10\n')
        file.write('    variable    = phi_zy\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
        file.write('  [./couple_33]\n')
        file.write(f'    type = {internal_couple}\n')
        file.write('    component_i = 2\n')
        file.write('    component_j = 2\n')
        file.write('    dof_num     = 11\n')
        file.write('    variable    = phi_zz\n')
        file.write('\n')
        file.write('    #Coupled variables\n')
        file.write('    u1     = disp_x\n')
        file.write('    u2     = disp_y\n')
        file.write('    u3     = disp_z\n')
        file.write('    phi_11 = phi_xx\n')
        file.write('    phi_22 = phi_yy\n')
        file.write('    phi_33 = phi_zz\n')
        file.write('    phi_23 = phi_yz\n')
        file.write('    phi_13 = phi_xz\n')
        file.write('    phi_12 = phi_xy\n')
        file.write('    phi_32 = phi_zy\n')
        file.write('    phi_31 = phi_zx\n')
        file.write('    phi_21 = phi_yx\n')
        file.write('  [../]\n')
    else:
        file.write('  [./couple_11]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_xx\n')
        file.write('  []\n')
        file.write('  [./couple_12]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_xy\n')
        file.write('  []\n')
        file.write('  [./couple_13]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_xz\n')
        file.write('  []\n')
        file.write('  [./couple_21]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_yx\n')
        file.write('  []\n')
        file.write('  [./couple_22]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_yy\n')
        file.write('  []\n')
        file.write('  [./couple_23]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_yz\n')
        file.write('  []\n')
        file.write('  [./couple_31]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_zx\n')
        file.write('  []\n')
        file.write('  [./couple_32]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_zy\n')
        file.write('  []\n')
        file.write('  [./couple_33]\n')
        file.write('    type = NullKernel\n')
        file.write('    variable = phi_zz\n')
        file.write('  []\n')
    

    return 0


def write_dynamic_kernels(file, ref_density):

    file.write('  [./inertia_kernel_1]\n')
    file.write('    type = MicromorphicInertialForce\n')
    file.write('    variable = disp_x\n')
    file.write('    component = 0\n')
    file.write('    dof_num = 0\n')
    file.write(f'    reference_density = {ref_density}\n')
    file.write('    u1 = disp_x\n')
    file.write('    u2 = disp_y\n')
    file.write('    u3 = disp_z\n')
    file.write('  []\n')
    file.write('  [./inertia_kernel_2]\n')
    file.write('    type = MicromorphicInertialForce\n')
    file.write('    variable = disp_y\n')
    file.write('    component = 1\n')
    file.write('    dof_num = 1\n')
    file.write(f'    reference_density = {ref_density}\n')
    file.write('    u1 = disp_x\n')
    file.write('    u2 = disp_y\n')
    file.write('    u3 = disp_z\n')
    file.write('  []\n')
    file.write('  [./inertia_kernel_3]\n')
    file.write('    type = MicromorphicInertialForce\n')
    file.write('    variable = disp_z\n')
    file.write('    component = 2\n')
    file.write('    dof_num = 2\n')
    file.write(f'    reference_density = {ref_density}\n')
    file.write('    u1 = disp_x\n')
    file.write('    u2 = disp_y\n')
    file.write('    u3 = disp_z\n')
    file.write('  []\n')

    return 0


def write_default_auxvariables(file):
    '''Write the default aux variables for forces and normal components of PK2 and Sigma stresses

    :params file file: The file to write to
    '''

    file.write('  [force_x][]\n')
    file.write('  [force_y][]\n')
    file.write('  [force_z][]\n')
    file.write('  [./pk2_11]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_22]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_33]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_11]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_22]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_33]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')

    return 0


def write_extra_second_order_auxvariables(file):
    '''Write the aux variables for off-diagonal components of PK2 and Sigma stresses

    :params file file: The file to write to
    '''

    file.write(' [./pk2_12]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_13]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_21]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_23]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_31]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./pk2_32]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_12]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_13]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./sigma_23]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')

    return 0


def write_higher_order_stress_auxvariables(file, dim=3):
    '''Write the aux variables for higher order stress components

    :params file file: The file to write to
    :params into dim: The problem dimension
    '''

    for i in range(0,dim):
        for j in range(0,dim):
            for k in range(0,dim):
                file.write(f'  [./M_{i+1}{j+1}{k+1}]\n')
                file.write('    order = CONSTANT\n')
                file.write('    family = MONOMIAL\n')
                file.write('  [../]\n')

    return 0


def write_plastic_auxvariables(file):
    '''Write the aux variables for plastic multipliers and internal state variables

    :params file file: The file to write to
    '''

    file.write('  [./macro_gamma]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gamma]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_gamma_1]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_gamma_2]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_gamma_3]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./macro_isv]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_isv]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_isv_1]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_isv_2]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')
    file.write('  [./micro_gradient_isv_3]\n')
    file.write('    order = CONSTANT\n')
    file.write('    family = MONOMIAL\n')
    file.write('  [../]\n')

    return 0


def write_default_auxkernels(file):
    '''Write the default aux kernels for normal components of PK2 and Sigma stresses

    :params file file: The file to write to
    '''

    file.write('[AuxKernels]\n')
    file.write('  [./pk2_11]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 0\n')
    file.write('    variable = pk2_11\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_22]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 4\n')
    file.write('    variable = pk2_22\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_33]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 8\n')
    file.write('    variable = pk2_33\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_11]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 0\n')
    file.write('    variable = sigma_11\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_22]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 4\n')
    file.write('    variable = sigma_22\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_33]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 8\n')
    file.write('    variable = sigma_33\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')

    return 0


def write_extra_second_order_stress_auxkernels(file):
    '''Write the off-diagonal aux kernels for normal components of PK2 and Sigma stresses

    :params file file: The file to write to
    '''

    file.write('[AuxKernels]\n')
    file.write('  [./pk2_12]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 1\n')
    file.write('    variable = pk2_12\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_13]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 2\n')
    file.write('    variable = pk2_13\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_21]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 3\n')
    file.write('    variable = pk2_21\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_23]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 5\n')
    file.write('    variable = pk2_23\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_31]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 6\n')
    file.write('    variable = pk2_31\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./pk2_32]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = PK2\n')
    file.write('    index = 7\n')
    file.write('    variable = pk2_32\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_12]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 1\n')
    file.write('    variable = sigma_12\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_13]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 3\n')
    file.write('    variable = sigma_13\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./sigma_23]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SIGMA\n')
    file.write('    index = 5\n')
    file.write('    variable = sigma_23\n')
    file.write('  [../]\n')
    file.write('[]\n')

    return 0


def write_plastic_auxkernels(file):
    '''Write the aux kernels for plastic multipliers and internal state variables

    :params file file: The file to write to
    '''

    file.write('## plastic Aux kernels\n')
    file.write('[AuxKernels]\n')
    file.write('  [./macro_gamma]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 45\n')
    file.write('    variable = macro_gamma\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gamma]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 46\n')
    file.write('    variable = micro_gamma\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_gamma_1]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 47\n')
    file.write('    variable = micro_gradient_gamma_1\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_gamma_2]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 48\n')
    file.write('    variable = micro_gradient_gamma_2\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_gamma_3]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 49\n')
    file.write('    variable = micro_gradient_gamma_3\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./macro_isv]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 50\n')
    file.write('    variable = macro_isv\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_isv]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 51\n')
    file.write('    variable = micro_isv\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_isv_1]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 52\n')
    file.write('    variable = micro_gradient_isv_1\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_isv_2]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 53\n')
    file.write('    variable = micro_gradient_isv_2\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')
    file.write('[AuxKernels]\n')
    file.write('  [./micro_gradient_isv_3]\n')
    file.write('    type = MaterialStdVectorAux\n')
    file.write('    property = SDVS\n')
    file.write('    index = 54\n')
    file.write('    variable = micro_gradient_isv_3\n')
    file.write('  [../]\n')
    file.write('[]\n')

    return 0


def write_higher_order_stress_auxkernels(file, dim=3):
    '''Write the aux kernels for higher order stress components

    :params file file: The file to write to
    :params into dim: The problem dimension
    '''

    m = 0
    for i in range(0,dim):
        for j in range(0,dim):
            for k in range(0,dim):
                file.write('[AuxKernels]\n')
                file.write(f'  [./M_{i+1}{j+1}{k+1}]\n')
                file.write('    type = MaterialStdVectorAux\n')
                file.write('    property = M\n')
                file.write(f'    index = {m}\n')
                file.write(f'    variable = M_{i+1}{j+1}{k+1}\n')
                file.write('  [../]\n')
                file.write('[]\n')
                m += 1

    return 0


def write_phi_BCs(file, phi_BC):
    '''Write the boundary conditions to fix micro-displacements

    :params file file: The file to write to
    :param str phi_BC: Nodeset to force micro deformation components to be zero
    '''


    phis = ['phi_xx', 'phi_xy', 'phi_xz',
            'phi_yx', 'phi_yy', 'phi_yz',
            'phi_zx', 'phi_zy', 'phi_zz']

    for phi in phis:
        file.write(f'  [fix_{phi}]\n')
        file.write('    type = DirichletBC\n')
        file.write(f'    variable = {phi}\n')
        file.write(f'    boundary = "{phi_BC}"\n')
        file.write('    value = 0\n')
        file.write('  [../]\n')

    return 0


def write_elastic_material_card(file, yaml_data, elem_num=None):
    '''Write an elastic material block

    :params file file: The file to write to
    :param dict yaml_data: Dictionary containing material information unpacked from a yaml file
    :param int elem_num: Optional element number to make the material card unique to a specific element block
    '''

    mat_line_1 = yaml_data['line 1']
    mat_line_2 = yaml_data['line 2']
    mat_line_3 = yaml_data['line 3']
    mat_line_4 = yaml_data['line 4']
    # Write in material info
    if elem_num is not None:
        file.write(f'  [./linear_elastic_{elem_num}]\n')
    else:
        file.write(f'  [./linear_elastic]\n')
    file.write('    type = MicromorphicMaterial\n')
    file.write(f'    material_fparameters = "{mat_line_1}\n')
    file.write(f'                            {mat_line_2}\n')
    file.write(f'                            {mat_line_3}\n')
    file.write(f'                            {mat_line_4}"\n')
    file.write(f'    model_name = "LinearElasticity"\n')
    file.write('\n')
    file.write('    #Coupled variables\n')
    file.write('    u1     = "disp_x"\n')
    file.write('    u2     = "disp_y"\n')
    file.write('    u3     = "disp_z"\n')
    file.write('    phi_11 = "phi_xx"\n')
    file.write('    phi_22 = "phi_yy"\n')
    file.write('    phi_33 = "phi_zz"\n')
    file.write('    phi_23 = "phi_yz"\n')
    file.write('    phi_13 = "phi_xz"\n')
    file.write('    phi_12 = "phi_xy"\n')
    file.write('    phi_32 = "phi_zy"\n')
    file.write('    phi_31 = "phi_zx"\n')
    file.write('    phi_21 = "phi_yx"\n')
    if elem_num is not None:
        file.write(f'    block = "element_{elem_num}"\n')
    file.write('  [../]\n')

    return 0


def write_plastic_material_card(file, yaml_data, elem_num=None):
    '''Write a plastic material block

    :params file file: The file to write to
    :param dict yaml_data: Dictionary containing material information unpacked from a yaml file
    :param int elem_num: Optional element number to make the material card unique to a specific element block
    '''

    mat_line_01 = yaml_data['line 01']
    mat_line_02 = yaml_data['line 02']
    mat_line_03 = yaml_data['line 03']
    mat_line_04 = yaml_data['line 04']
    mat_line_05 = yaml_data['line 05']
    mat_line_06 = yaml_data['line 06']
    mat_line_07 = yaml_data['line 07']
    mat_line_08 = yaml_data['line 08']
    mat_line_09 = yaml_data['line 09']
    mat_line_10 = yaml_data['line 10']
    mat_line_11 = yaml_data['line 11']
    mat_line_12 = yaml_data['line 12']
    mat_line_13 = yaml_data['line 13']
    mat_line_14 = yaml_data['line 14']
    # Write in material info
    if elem_num is not None:
        file.write(f'  [./linear_elastic_{elem_num}]\n')
    else:
        file.write(f'  [./linear_elastic]\n')
    file.write('    type = MicromorphicMaterial\n')
    file.write(f'    material_fparameters = "{mat_line_01}\n')
    file.write(f'                            {mat_line_02}\n')
    file.write(f'                            {mat_line_03}\n')
    file.write(f'                            {mat_line_04}\n')
    file.write(f'                            {mat_line_05}\n')
    file.write(f'                            {mat_line_06}\n')
    file.write(f'                            {mat_line_07}\n')
    file.write(f'                            {mat_line_08}\n')
    file.write(f'                            {mat_line_09}\n')
    file.write(f'                            {mat_line_10}\n')
    file.write(f'                            {mat_line_11}\n')
    file.write(f'                            {mat_line_12}\n')
    file.write(f'                            {mat_line_13}\n')
    file.write(f'                            {mat_line_14}"\n')
    file.write('    number_SDVS = 55\n')
    file.write(f'    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
    file.write('\n')
    file.write('    #Coupled variables\n')
    file.write('    u1     = "disp_x"\n')
    file.write('    u2     = "disp_y"\n')
    file.write('    u3     = "disp_z"\n')
    file.write('    phi_11 = "phi_xx"\n')
    file.write('    phi_22 = "phi_yy"\n')
    file.write('    phi_33 = "phi_zz"\n')
    file.write('    phi_23 = "phi_yz"\n')
    file.write('    phi_13 = "phi_xz"\n')
    file.write('    phi_12 = "phi_xy"\n')
    file.write('    phi_32 = "phi_zy"\n')
    file.write('    phi_31 = "phi_zx"\n')
    file.write('    phi_21 = "phi_yx"\n')
    if elem_num is not None:
        file.write(f'    block = "element_{elem_num}"\n')
    file.write('  [../]\n')

    return 0


def write_preconditioner_block(file):
    '''Write the default SMP preconditioner block

    :params file file: The file to write to
    '''

    file.write('[Preconditioning]\n')
    file.write('  [./SMP]\n')
    file.write('    type = SMP\n')
    file.write('    full = true\n')
    file.write('  [../]\n')
    file.write('[]\n')
    file.write('\n')

    return 0


def write_outputs_block(file):
    '''Write the default outputs block

    :params file file: The file to write to
    '''

    file.write('[Outputs]\n')
    file.write('  exodus = true\n')
    file.write('  perf_graph = true\n')
    file.write('  csv = true\n')
    file.write('  [./console]\n')
    file.write('    type = Console\n')
    file.write('  [../]\n')
    file.write('[]\n')

    return 0
