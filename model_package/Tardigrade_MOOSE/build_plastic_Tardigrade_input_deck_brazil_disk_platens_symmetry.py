#!python
import argparse
import inspect
import os
import pathlib
import sys
import yaml


def build_input(output_file, mesh_file, parameter_sets, disp, duration,
                specimen_bottom_surface, bottom_platen_contact,
                top_symmetry, back_symmetry, side_symmetry,
                bottom_platen_fixture):
    '''Write Tardigrade-MOOSE input file for eighth symmetry Brazilian disk simulation with platens
    
    :param str output_file: The name of Tardigrade-MOOSE file to write
    :param str mesh_file: The name of the mesh file
    :param list parameter_sets: The list of yaml files containing calibration results
    :param float disp: The compressive displacement to be applied
    :param float duration: The duration of the simulation
    :param str specimen_bottom_surface: The name of the specimen bottom contact surface
    :param str bottom_platen_contact: The name of the bottom platen contact surface
    :param str top_symmetry: The name of the top symmetry surface(s)
    :param str back_symmetry: The name of the back symmetry surface(s)
    :param str side_symmetry: The name of the side symmetry surface(s)
    :param str bottom_platen_fixture: The name of the bottom platen fixture surface

    :returns: ``output_file``
    '''

    assert os.path.exists(mesh_file), f"Mesh file not found: {mesh_file}"

    # Write input file
    with open(output_file, 'w') as f:
        f.write('###############################################################################\n')
        f.write('[Mesh]\n')
        f.write('  type = FileMesh\n')
        f.write(f'  file = "{mesh_file}"\n')
        f.write('  patch_update_strategy = iteration\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[GlobalParams]\n')
        f.write('  displacements = "disp_x disp_y disp_z"\n')
        f.write('[]\n')
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
        f.write('  [./inc_slip_bottom_x]\n')
        f.write('  [../]\n')
        f.write('  [./inc_slip_bottom_y]\n')
        f.write('  [../]\n')
        f.write('  [./inc_slip_bottom_z]\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_x]\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_y]\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_z]\n')
        f.write('  [../]\n')
        f.write('  [./tang_force_bottom_x]\n')
        f.write('  [../]\n')
        f.write('  [./tang_force_bottom_y]\n')
        f.write('  [../]\n')
        f.write('  [./tang_force_bottom_z]\n')
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
        f.write('  # Contact Kernels\n')
        f.write('  [./inc_slip_bottom_x]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = inc_slip_bottom_x\n')
        f.write('    quantity = incremental_slip_x\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./inc_slip_bottom_y]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = inc_slip_bottom_y\n')
        f.write('    quantity = incremental_slip_y\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./inc_slip_bottom_z]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = inc_slip_bottom_z\n')
        f.write('    quantity = incremental_slip_z\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./tangential_force_bottom_x]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = tang_force_bottom_x\n')
        f.write('    execute_on = timestep_end\n')
        f.write('    quantity = tangential_force_x\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./tangential_force_bottom_y]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = tang_force_bottom_y\n')
        f.write('    execute_on = timestep_end\n')
        f.write('    quantity = tangential_force_y\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./tangential_force_bottom_z]\n')
        f.write('    type = PenetrationAux\n')
        f.write('    variable = tang_force_bottom_z\n')
        f.write('    execute_on = timestep_end\n')
        f.write('    quantity = tangential_force_z\n')
        f.write(f'    boundary = "{specimen_bottom_surface}"\n')
        f.write(f'    paired_boundary = "{bottom_platen_contact}"\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_x]\n')
        f.write('    type = AccumulateAux\n')
        f.write('    variable = accum_slip_bottom_x\n')
        f.write('    accumulate_from_variable = inc_slip_bottom_x\n')
        f.write('    execute_on = timestep_end\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_y]\n')
        f.write('    type = AccumulateAux\n')
        f.write('    variable = accum_slip_bottom_y\n')
        f.write('    accumulate_from_variable = inc_slip_bottom_y\n')
        f.write('    execute_on = timestep_end\n')
        f.write('  [../]\n')
        f.write('  [./accum_slip_bottom_z]\n')
        f.write('    type = AccumulateAux\n')
        f.write('    variable = accum_slip_bottom_z\n')
        f.write('    accumulate_from_variable = inc_slip_bottom_z\n')
        f.write('    execute_on = timestep_end\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        # Reaction Force
        sample_force = 'force_y'
        f.write('# Do some cool math to get the reaction force\n')
        f.write('[Postprocessors]\n')
        f.write('  [bot_react_y]\n')
        f.write('    type = NodalSum\n')
        f.write(f'    variable = {sample_force}\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('  []\n')
        f.write('[]\n')
        f.write('\n')
        # BCs
        f.write('[BCs]\n')
        f.write('  active = "top_symmetry back_symmetry side_symmetry bottom_y"\n')
        f.write('  [./top_symmetry]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_y\n')
        f.write(f'    boundary = "{top_symmetry}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./back_symmetry]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_z\n')
        f.write(f'    boundary = "{back_symmetry}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./side_symmetry]\n')
        f.write('    type = DirichletBC\n')
        f.write('    variable = disp_x\n')
        f.write(f'    boundary = "{side_symmetry}"\n')
        f.write('    preset = true\n')
        f.write('    value = 0\n')
        f.write('  [../]\n')
        f.write('  [./bottom_y]\n')
        f.write('    type = FunctionDirichletBC\n')
        f.write('    variable = disp_y\n')
        f.write(f'    boundary = "{bottom_platen_fixture}"\n')
        f.write('    preset = true\n')
        f.write('    function = top_bc\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        # Loading function
        f.write('[Functions]\n')
        f.write('  [./top_bc]\n')
        f.write('    type  = ParsedFunction\n')
        f.write(f'    expression = {disp}*t\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        # Contact
        f.write('[Contact]\n')
        f.write('  [./bottom_center_cont]\n')
        f.write(f'    secondary = "{bottom_platen_contact}"\n')
        f.write(f'    primary = "{specimen_bottom_surface}"\n')
        f.write('    model = coulomb\n')
        f.write('    formulation = tangential_penalty\n')
        f.write('    friction_coefficient = "0.2"\n')
        f.write('    penalty = 1e4\n')
        f.write('    normalize_penalty = true\n')
        f.write('    tangential_tolerance = 1.e-1\n')
        f.write('    normal_smoothing_distance = 0.001\n')
        f.write('  [../]\n')
        f.write('[]\n')
        f.write('\n')
        f.write('[Dampers]\n')
        f.write('  [./contact_slip_bottom]\n')
        f.write('    type = ContactSlipDamper\n')
        f.write(f'    secondary = "{bottom_platen_contact}"\n')
        f.write(f'    primary = "{specimen_bottom_surface}"\n')
        f.write('    min_damping_factor = 1.e-2\n')
        f.write('  [../]\n')
        f.write('[]\n')
        # Materials
        f.write('[Materials]\n')
        # Load in parameter data for each filter domain / element
        if len(parameter_sets) > 1:
            for i, set in enumerate(parameter_sets):
                # Load yaml file
                stream = open(set, 'r')
                UI = yaml.load(stream, Loader=yaml.FullLoader)
                stream.close()
                mat_line_01 = UI['line 01']
                mat_line_02 = UI['line 02']
                mat_line_03 = UI['line 03']
                mat_line_04 = UI['line 04']
                mat_line_05 = UI['line 05']
                mat_line_06 = UI['line 06']
                mat_line_07 = UI['line 07']
                mat_line_08 = UI['line 08']
                mat_line_09 = UI['line 09']
                mat_line_10 = UI['line 10']
                mat_line_11 = UI['line 11']
                mat_line_12 = UI['line 12']
                mat_line_13 = UI['line 13']
                mat_line_14 = UI['line 14']
                # Write in material info
                f.write(f'  [./linear_elastic_{i}]\n')
                f.write('    type = MicromorphicMaterial\n')
                f.write(f'    material_fparameters = "{mat_line_01}\n')
                f.write(f'                            {mat_line_02}\n')
                f.write(f'                            {mat_line_03}\n')
                f.write(f'                            {mat_line_04}\n')
                f.write(f'                            {mat_line_05}\n')
                f.write(f'                            {mat_line_06}\n')
                f.write(f'                            {mat_line_07}\n')
                f.write(f'                            {mat_line_08}\n')
                f.write(f'                            {mat_line_09}\n')
                f.write(f'                            {mat_line_10}\n')
                f.write(f'                            {mat_line_11}\n')
                f.write(f'                            {mat_line_12}\n')
                f.write(f'                            {mat_line_13}\n')
                f.write(f'                            {mat_line_14}"\n')
                f.write('    number_SDVS = 55\n')
                f.write(f'    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
                f.write('\n')
                f.write('    #Coupled variables\n')
                f.write('    u1     = "disp_x"\n')
                f.write('    u2     = "disp_y"\n')
                f.write('    u3     = "disp_z"\n')
                f.write('    phi_11 = "phi_xx"\n')
                f.write('    phi_22 = "phi_yy"\n')
                f.write('    phi_33 = "phi_zz"\n')
                f.write('    phi_23 = "phi_yz"\n')
                f.write('    phi_13 = "phi_xz"\n')
                f.write('    phi_12 = "phi_xy"\n')
                f.write('    phi_32 = "phi_zy"\n')
                f.write('    phi_31 = "phi_zx"\n')
                f.write('    phi_21 = "phi_yx"\n')
                f.write(f'    block = "element_{i}"\n')
                f.write('  [../]\n')
        else:
            # Load yaml file
            set = parameter_sets[0]
            stream = open(set, 'r')
            UI = yaml.load(stream, Loader=yaml.FullLoader)
            stream.close()
            mat_line_01 = UI['line 01']
            mat_line_02 = UI['line 02']
            mat_line_03 = UI['line 03']
            mat_line_04 = UI['line 04']
            mat_line_05 = UI['line 05']
            mat_line_06 = UI['line 06']
            mat_line_07 = UI['line 07']
            mat_line_08 = UI['line 08']
            mat_line_09 = UI['line 09']
            mat_line_10 = UI['line 10']
            mat_line_11 = UI['line 11']
            mat_line_12 = UI['line 12']
            mat_line_13 = UI['line 13']
            mat_line_14 = UI['line 14']
            # Write in material info
            f.write(f'  [./linear_elastic]\n')
            f.write('    type = MicromorphicMaterial\n')
            f.write(f'    material_fparameters = "{mat_line_01}\n')
            f.write(f'                            {mat_line_02}\n')
            f.write(f'                            {mat_line_03}\n')
            f.write(f'                            {mat_line_04}\n')
            f.write(f'                            {mat_line_05}\n')
            f.write(f'                            {mat_line_06}\n')
            f.write(f'                            {mat_line_07}\n')
            f.write(f'                            {mat_line_08}\n')
            f.write(f'                            {mat_line_09}\n')
            f.write(f'                            {mat_line_10}\n')
            f.write(f'                            {mat_line_11}\n')
            f.write(f'                            {mat_line_12}\n')
            f.write(f'                            {mat_line_13}\n')
            f.write(f'                            {mat_line_14}"\n')
            f.write('    number_SDVS = 55\n')
            f.write(f'    model_name = "LinearElasticityDruckerPragerPlasticity"\n')
            f.write('\n')
            f.write('    #Coupled variables\n')
            f.write('    u1     = "disp_x"\n')
            f.write('    u2     = "disp_y"\n')
            f.write('    u3     = "disp_z"\n')
            f.write('    phi_11 = "phi_xx"\n')
            f.write('    phi_22 = "phi_yy"\n')
            f.write('    phi_33 = "phi_zz"\n')
            f.write('    phi_23 = "phi_yz"\n')
            f.write('    phi_13 = "phi_xz"\n')
            f.write('    phi_12 = "phi_xy"\n')
            f.write('    phi_32 = "phi_zy"\n')
            f.write('    phi_31 = "phi_zx"\n')
            f.write('    phi_21 = "phi_yx"\n')
            f.write('    block = "specimen"\n')
            f.write('  [../]\n')
        f.write('  [./platen]\n')
        f.write('    type = MicromorphicMaterial\n')
        f.write('    material_fparameters = "2 121200. 80770.\n')
        f.write('                            5 0. 0. 0. 0. 0.\n')
        f.write('                            11 0. 0. 0. 0. 0. 0. 0.001 0. 0. 0. 0.\n')
        f.write('                            2 0. 0."\n')
        f.write('    model_name = "LinearElasticity"\n')
        f.write('\n')
        f.write('    #Coupled variables\n')
        f.write('    u1     = "disp_x"\n')
        f.write('    u2     = "disp_y"\n')
        f.write('    u3     = "disp_z"\n')
        f.write('    phi_11 = "phi_xx"\n')
        f.write('    phi_22 = "phi_yy"\n')
        f.write('    phi_33 = "phi_zz"\n')
        f.write('    phi_23 = "phi_yz"\n')
        f.write('    phi_13 = "phi_xz"\n')
        f.write('    phi_12 = "phi_xy"\n')
        f.write('    phi_32 = "phi_zy"\n')
        f.write('    phi_31 = "phi_zx"\n')
        f.write('    phi_21 = "phi_yx"\n')
        f.write('    block = "bottom_platen"\n')
        f.write('  [../]\n')
        f.write('[]\n')
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
        #f.write('  solve_type = NEWTON\n')
        f.write('  solve_type = PJFNK\n')
        #f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package"\n')
        #f.write('  petsc_options_value = "lu       superlu_dist                 "\n')
        f.write('  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package -pc_factor_shift_type"\n')
        f.write('  petsc_options_value = "lu     superlu_dist NONZERO"\n')
        f.write('  line_search = none\n')
        f.write('  automatic_scaling = true\n')
        f.write('  nl_rel_tol = 1e-5\n')
        f.write('  nl_abs_tol = 1e-5\n')
        f.write('  l_tol = 5e-3\n')
        f.write('  l_max_its = 250\n')
        f.write('  nl_max_its = 150\n')
        f.write('  start_time = 0.0\n')
        f.write(f'  end_time = {duration}\n')
        f.write('  dtmin = 1e-10\n')
        f.write('  dtmax= 0.1\n')
        f.write('  [TimeStepper]\n')
        f.write('    type = IterationAdaptiveDT\n')
        #f.write('    optimal_iterations = 4\n')
        #f.write('    iteration_window = 3\n')
        #f.write('    linear_iteration_ratio = 1000\n')
        f.write('    growth_factor=1.5\n')
        #f.write('    cutback_factor=0.5\n')
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
    cli_description = "Write Tardigrade-MOOSE input file for eighth symmetry Brazilian disk simulation with platens"
    parser = argparse.ArgumentParser(description=cli_description, prog=prog)
    parser.add_argument('-o', '--output-file', type=str, required=True,
        help="Specify the name of Tardigrade-MOOSE file to write")
    parser.add_argument('--mesh', type=str, required=True,
        help='Specify the mesh file')
    parser.add_argument('--parameter-sets', nargs="+", required=True,
        help='Specify the list of yaml files containing calibration results')
    parser.add_argument('--disp', type=float, required=True,
        help='Specify the compressive displacement to be applied')
    parser.add_argument('--duration', type=float, required=True,
        help='Specify the duration of the simulation')
    parser.add_argument('--specimen-bottom-surface', type=str, required=True,
        help='Specify the name of the specimen bottom contact surface')
    parser.add_argument('--bottom-platen-contact', type=str, required=True,
        help='Specify the name of the bottom platen contact surface')
    parser.add_argument('--top-symmetry', type=str, required=True,
        help='Specify the name of the top symmetry surface(s)')
    parser.add_argument('--back-symmetry', type=str, required=True,
        help='Specify the name of the back symmetry surface(s)')
    parser.add_argument('--side-symmetry', type=str, required=True,
        help='Specify the name of the side symmetry surface(s)')
    parser.add_argument('--bottom-platen-fixture', type=str, required=True,
        help='Specify the name of the bottom platen fixture surface')

    return parser


if __name__ == '__main__':
    parser = get_parser()
    
    args, unknown = parser.parse_known_args()
    sys.exit(build_input(output_file=args.output_file,
                         mesh_file=args.mesh,
                         parameter_sets=args.parameter_sets,
                         disp=args.disp,
                         duration=args.duration,
                         specimen_bottom_surface=args.specimen_bottom_surface,
                         bottom_platen_contact=args.bottom_platen_contact,
                         top_symmetry=args.top_symmetry,
                         back_symmetry=args.back_symmetry,
                         side_symmetry=args.side_symmetry,
                         bottom_platen_fixture=args.bottom_platen_fixture,
                         ))
