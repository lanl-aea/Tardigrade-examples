elastic_cylinder = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'material_E': 250.,
    'material_nu': 0.2,
    'material_rho': 2.0e-9,
    'cut': True,
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
    'mesh_copy_root': 'cylinder_5_5',
}

dynamic_elastic_cylinder = {
    'diam': 5.0,
    'height': 5.0,
    'material_E': 250.,
    'material_nu': 0.0,
    'material_rho': 1.935e-9,
    'pressure': 1.273239545,
    'duration' : 1.669251329e-4,
    'num_steps' : 60,
    'finite_rise': 0,
    'macro_BC': 'slip',
    'mesh_copy_root': 'cylinder_5_5',
}
