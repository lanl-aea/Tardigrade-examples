#!python
import argparse
import os
import pathlib
import sys


def build_input(output_file, mesh_file, parameter_csv, BCs, disp, duration):
    '''Write Tardigrade-MOOSE input file for a gradient-enhanced damage plasticity simulation
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param list parameter_csv: The csv file containing unique calibrations for each element
    :param str BCs: The type of boundary conditions, either "slip" or "clamp"
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # Write input file
    with open(output_file, 'w') as f:
        f.write('###############################################################################\n')
        f.write('[Mesh]\n')
        f.write('  type = FileMesh\n')
        f.write('  displacements = "disp_x disp_y disp_z"\n')
        f.write('  dim = 3\n')
        f.write(f'  file = "{mesh_file}"\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[Variables]\n')
        f.write('  [./disp_x]\n')
        f.write('  [../]\n')
        f.write('  [./disp_y]\n')
        f.write('  [../]\n')
        f.write('  [./disp_z]\n')
        f.write('  [../]\n')
        f.write('  [./phi_xx]\n')
        f.write('  [../]\n')
        f.write('  [./phi_yy]\n')
        f.write('  [../]\n')
        f.write('  [./phi_zz]\n')
        f.write('  [../]\n')
        f.write('  [./phi_yz]\n')
        f.write('  [../]\n')
        f.write('  [./phi_xz]\n')
        f.write('  [../]\n')
        f.write('  [./phi_xy]\n')
        f.write('  [../]\n')
        f.write('  [./phi_zy]\n')
        f.write('  [../]\n')
        f.write('  [./phi_zx]\n')
        f.write('  [../]\n')
        f.write('  [./phi_yx]\n')
        f.write('  [../]\n')
        f.write('  [nonlocal_damage]\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[Kernels]\n')
        f.write('  #Define the internal force balance equations\n')
        f.write('  [./force_1]\n')
        f.write('    type = InternalForce\n')
        f.write('    component = 0\n')
        f.write('    dof_num   = 0\n')
        f.write('    variable  = disp_x\n')
        f.write('    save_in = force_x\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./force_2]\n')
        f.write('    type = InternalForce\n')
        f.write('    component = 1\n')
        f.write('    dof_num   = 1\n')
        f.write('    variable  = disp_y\n')
        f.write('    save_in = force_y\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./force_3]\n')
        f.write('    type = InternalForce\n')
        f.write('    component = 2\n')
        f.write('    dof_num   = 2\n')
        f.write('    variable  = disp_z\n')
        f.write('    save_in = force_z\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  #Define the internal couple balance equations\n')
        f.write('  [./couple_11]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 0\n')
        f.write('    component_j = 0\n')
        f.write('    dof_num     = 3\n')
        f.write('    variable    = phi_xx\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_12]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 0\n')
        f.write('    component_j = 1\n')
        f.write('    dof_num     = 4\n')
        f.write('    variable    = phi_xy\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_13]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 0\n')
        f.write('    component_j = 2\n')
        f.write('    dof_num     = 5\n')
        f.write('    variable    = phi_xz\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_21]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 1\n')
        f.write('    component_j = 0\n')
        f.write('    dof_num     = 6\n')
        f.write('    variable    = phi_yx\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_22]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 1\n')
        f.write('    component_j = 1\n')
        f.write('    dof_num     = 7\n')
        f.write('    variable    = phi_yy\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_23]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 1\n')
        f.write('    component_j = 2\n')
        f.write('    dof_num     = 8\n')
        f.write('    variable    = phi_yz\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_31]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 2\n')
        f.write('    component_j = 0\n')
        f.write('    dof_num     = 9\n')
        f.write('    variable    = phi_zx\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_32]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 2\n')
        f.write('    component_j = 1\n')
        f.write('    dof_num     = 10\n')
        f.write('    variable    = phi_zy\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('  [./couple_33]\n')
        f.write('    type = InternalCouple\n')
        f.write('    component_i = 2\n')
        f.write('    component_j = 2\n')
        f.write('    dof_num     = 11\n')
        f.write('    variable    = phi_zz\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = disp_x\n')
        f.write('    u2     = disp_y\n')
        f.write('    u3     = disp_z\n')
        f.write('    phi_11 = phi_xx\n')
        f.write('    phi_22 = phi_yy\n')
        f.write('    phi_33 = phi_zz\n')
        f.write('    phi_23 = phi_yz\n')
        f.write('    phi_13 = phi_xz\n')
        f.write('    phi_12 = phi_xy\n')
        f.write('    phi_32 = phi_zy\n')
        f.write('    phi_31 = phi_zx\n')
        f.write('    phi_21 = phi_yx\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('\n')
        f.write('[AuxVariables]\n')
        f.write('  [force_x][]\n')
        f.write('  [force_y][]\n')
        f.write('  [force_z][]\n')
        f.write('  [./pk2_11]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./pk2_22]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./pk2_33]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./sigma_11]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./sigma_22]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./sigma_33]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./macro_gamma]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gamma]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_gamma_1]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_gamma_2]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_gamma_3]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./macro_isv]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_isv]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_isv_1]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_isv_2]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('  [./micro_gradient_isv_3]\n')
        f.write('    order = CONSTANT\n')
        f.write('    family = MONOMIAL\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./pk2_11]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = PK2\n')
        f.write('    index = 0\n')
        f.write('    variable = pk2_11\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./pk2_22]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = PK2\n')
        f.write('    index = 4\n')
        f.write('    variable = pk2_22\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./pk2_33]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = PK2\n')
        f.write('    index = 8\n')
        f.write('    variable = pk2_33\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./sigma_11]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SIGMA\n')
        f.write('    index = 0\n')
        f.write('    variable = sigma_11\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./sigma_22]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SIGMA\n')
        f.write('    index = 4\n')
        f.write('    variable = sigma_22\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./sigma_33]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SIGMA\n')
        f.write('    index = 8\n')
        f.write('    variable = sigma_33\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./macro_gamma]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 45\n')
        f.write('    variable = macro_gamma\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gamma]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 46\n')
        f.write('    variable = micro_gamma\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_gamma_1]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 47\n')
        f.write('    variable = micro_gradient_gamma_1\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_gamma_2]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 48\n')
        f.write('    variable = micro_gradient_gamma_2\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_gamma_3]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 49\n')
        f.write('    variable = micro_gradient_gamma_3\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./macro_isv]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 50\n')
        f.write('    variable = macro_isv\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_isv]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 51\n')
        f.write('    variable = micro_isv\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_isv_1]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 52\n')
        f.write('    variable = micro_gradient_isv_1\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_isv_2]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 53\n')
        f.write('    variable = micro_gradient_isv_2\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[AuxKernels]\n')
        f.write('  [./micro_gradient_isv_3]\n')
        f.write('    type = MaterialStdVectorAux\n')
        f.write('    property = SDVS\n')
        f.write('    index = 54\n')
        f.write('    variable = micro_gradient_isv_3\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        if BCs == 'brazil':
            sample_force = 'force_y'
            sample_boundary = 'brazil_load'
        else:
            sample_force = 'force_z'
            sample_boundary = 'top'
        f.write('# Do some cool math to get the reaction force\n')
        f.write('[Postprocessors]\n')
        f.write('  [bot_react_z]\n')
        f.write('    type = NodalSum\n')
        f.write(f'    variable = {sample_force}\n')
        f.write(f'    boundary = "{sample_boundary}"\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        if BCs == 'slip':
            f.write('[BCs]\n')
            f.write('  active = "x_symm y_symm bottom_z top_z"\n')
            f.write('  [./x_symm]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "x_plane"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./y_symm]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "y_plane"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_z]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
            f.write('[]\n')
        elif BCs == 'clamp':
            f.write('[BCs]\n')
            f.write('  active = "bottom_x bottom_y bottom_z top_x top_y top_z"\n')
            f.write('  [./bottom_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./top_z]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "top"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
            f.write('[]\n')
        elif BCs == 'brazil':
            f.write('[BCs]\n')
            f.write('  active = "bottom_z brazil_fix_x brazil_fix_y brazil_load_x brazil_load_y"\n')
            f.write('  [./bottom_z]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_z\n')
            f.write('    boundary = "bottom"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_fix_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "brazil_fix"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_fix_y]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "brazil_fix"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_load_x]\n')
            f.write('    type = DirichletBC\n')
            f.write('    variable = disp_x\n')
            f.write('    boundary = "brazil_load"\n')
            f.write('    preset = true\n')
            f.write('    value = 0\n')
            f.write('  [../]\n')
            f.write('  [./brazil_load_y]\n')
            f.write('    type = FunctionDirichletBC\n')
            f.write('    variable = disp_y\n')
            f.write('    boundary = "brazil_load"\n')
            f.write('    preset = true\n')
            f.write('    function = top_bc\n')
            f.write('  [../]\n')
            f.write('[]\n')
        else:
            print('Specify a valid BC type!')
        f.write('\n')
        f.write('[Functions]\n')
        f.write('  [./top_bc]\n')
        f.write('    type  = ParsedFunction\n')
        f.write(f'    expression = -{disp}*t\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[Materials]\n')
        # Write in material info
        f.write('  [./linear_elastic]\n')
        f.write('    type = GradientEnhancedDamagedMicromorphicMaterial\n')
        f.write('    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
        f.write('    material_fparameters = "2 1.0 2.0\n')
        f.write('                            2 4.0 5.0\n')
        f.write('                            2 6.0 7.0\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 0. 0.\n')
        f.write('                            2 28. 29.\n')
        f.write('                            5 31. 32. 33. 34. 35.\n')
        f.write('                            11 37. 38. 39. 40. 41. 42. 43. 44. 45. 46. 47.\n')
        f.write('                            2 32. 35.\n')
        f.write('                            0.5 0.5 0.5 1e-9 1e-9"\n')
        f.write('\n')
        f.write('    user_material_prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11"\n')
        f.write('    user_material_prop_indices = "1 2 4 5 7 8 28 29 31 32 33 34 35 37 38 39 40 41 42 43 44 45 46 47"\n')
        f.write('\n')
        f.write('    number_SDVS = 55\n')
        #f.write('    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
        f.write('\n')
        f.write(f'    gradient_enhanced_damage_fparameters = "15 0.095 1.0 1"\n')
        #f.write('    #Coupled variables\n')
        #f.write('    u1     = "disp_x"\n')
        #f.write('    u2     = "disp_y"\n')
        #f.write('    u3     = "disp_z"\n')
        #f.write('    phi_11 = "phi_xx"\n')
        #f.write('    phi_22 = "phi_yy"\n')
        #f.write('    phi_33 = "phi_zz"\n')
        #f.write('    phi_23 = "phi_yz"\n')
        #f.write('    phi_13 = "phi_xz"\n')
        #f.write('    phi_12 = "phi_xy"\n')
        #f.write('    phi_32 = "phi_zy"\n')
        #f.write('    phi_31 = "phi_zx"\n')
        #f.write('    phi_21 = "phi_yx"\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        # Specify csv file
        f.write('[UserObjects]\n')
        f.write('  [reader_element]\n')
        f.write('    type = PropertyReadFile\n')
        f.write(f'    prop_file_name = "{parameter_csv}"\n')
        f.write('    read_type = "element"\n')
        f.write('    nprop = 24\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        # Set up function material map
        f.write('[Materials]\n')
        f.write('  [E_nu]\n')
        f.write('    type = GenericFunctionMaterial\n')
        f.write('    prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11"\n')
        f.write('    prop_values = "func_cu0 func_Hu func_cchi0 func_Hchi func_cnablachi0 func_Hnablachi func_lambda func_mu func_eta func_tau func_kappa func_nu func_sigma func_tau1 func_tau2 func_tau3 func_tau4 func_tau5 func_tau6 func_tau7 func_tau8 func_tau9 func_tau10 func_tau11"\n')
        f.write('    outputs = exodus\n')
        f.write('    []\n')
        f.write('[]\n')
        f.write('\n')
        # assign parameters
        f.write('[Functions]\n')
        parameters = ['cu0', 'Hu', 'cchi0', 'Hchi', 'cnablachi0', 'Hnablachi',
                      'lambda', 'mu',
                       'eta', 'tau', 'kappa', 'nu', 'sigma',
                       'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7', 'tau8', 'tau9', 'tau10', 'tau11']
        i = 0
        for p in parameters:
            f.write(f'  [func_{p}]\n')
            f.write('    type = PiecewiseConstantFromCSV\n')
            f.write('    read_prop_user_object = "reader_element"\n')
            f.write('    read_type = "element"\n')
            f.write(f'    column_number = "{i}"\n')
            f.write('  []\n')
            i = i + 1
        f.write('[]\n')
        # Solver setup
        f.write('\n')
        f.write('[Preconditioning]\n')
        f.write('  [./SMP]\n')
        f.write('    type = SMP\n')
        f.write('#    type = FDP\n')
        f.write('    full = true\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        dt = duration / 100
        f.write('[Executioner]\n')
        f.write('  type = Transient\n')
        f.write('  solve_type = NEWTON\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package"\n')
        f.write('  petsc_options_value = "lu       superlu_dist                 "\n')
        f.write('  line_search = none\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-5\n')
        f.write('  nl_abs_tol = 1e-7\n')
        f.write('  nl_max_its = 50\n')
        f.write('  start_time = 0.0\n')
        f.write('  end_time = 1.0\n')
        f.write('  dtmin = 1e-12\n')
        f.write('  dtmax= 0.1\n')
        f.write('  [TimeStepper]\n')
        f.write('    type = IterationAdaptiveDT\n')
        f.write('    optimal_iterations = 4\n')
        f.write('    iteration_window = 3\n')
        f.write('    linear_iteration_ratio = 1000\n')
        f.write('    growth_factor=1.1\n')
        f.write('    cutback_factor=0.5\n')
        f.write(f'    dt = {dt}\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('[Outputs]\n')
        f.write('  exodus = true\n')
        f.write('  perf_graph = true\n')
        f.write('  csv = true\n')
        f.write('  [./console]\n')
        f.write('    type = Console\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')

    return 0


def get_parser():

    script_name = pathlib.Path(__file__)

    prog = f"python {script_name.name} "
    cli_description = "Write Tardigrade-MOOSE input file for a gradient-enhanced damage plasticity simulation"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--parameter-csv', nargs="+", required=True,
        help='The csv file containing unique calibrations for each element')
    parser.add_argument('--BCs', type=str, required=True,
        help='Specify the type of boundary conditions, either "slip" or "clamp"')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         parameter_csv=args.parameter_csv,
                         BCs=args.BCs,
                         disp=args.disp,
                         duration=args.duration,
                         ))
