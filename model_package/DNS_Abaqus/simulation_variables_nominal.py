elastic_cylinder = {
    'diam' : 5.0,
    'height' : 5.0,
    'seed' : 0.5,
    'material_E' : 165.0,
    'material_nu' : 0.39,
    'material_rho' : 2.00,
    'disp' : 0.01,
    'num_steps' : 5,
    'block_name' : 'CYLINDER-1',
    'collocation_option' : 'ip',
    'cut': True,
    'micro_BC': 'slip',
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'cylinder_5_5',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 1,
    'UQ_file': False,
    'ignore_boundary': False,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
}

elastic_cylinder_clamp = {
    'diam' : 5.0,
    'height' : 5.0,
    'seed' : 0.25,
    'material_E' : 165.0,
    'material_nu' : 0.39,
    'material_rho' : 2.00,
    'disp' : 0.01,
    'num_steps' : 5,
    'block_name' : 'CYLINDER-1',
    'collocation_option' : 'ip',
    'cut': True,
    'micro_BC': 'clamp',
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'cylinder_5_5',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 1,
    'UQ_file': False,
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}

dynamic_elastic_cylinder = {
    'diam' : 5.0,
    'height' : 5.0,
    'seed' : 0.25,
    'material_E' : 302.4,
    'material_nu' : 0.0,
    'material_rho' : 1.890,
    'total_force' : 29.6881,
    'duration' : 1.5e-4,
    'num_steps' : 600,
    'block_name' : 'CYLINDER-1',
    'collocation_option' : 'ip',
    'cut': True,
    'finite_rise': 5,
    'field_frame': 100,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'cylinder_5_5',
    # parameters for micromorphic filter
    'acceleration': True,
    'velocity': True,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 1,
    'calibration_increment': 4,
    # paramters for Tardigrade-MOOSE
    'macro_BC': 'slip',
    'macro_disp': 0.05,
    'macro_duration': 1.0,
}
