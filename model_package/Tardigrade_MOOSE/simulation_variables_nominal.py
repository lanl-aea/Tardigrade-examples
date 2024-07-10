elastic_cylinder = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'material_E': 250.,
    'material_nu': 0.2,
    'material_rho': 2.0e-9,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'cylinder_5_5',
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
}