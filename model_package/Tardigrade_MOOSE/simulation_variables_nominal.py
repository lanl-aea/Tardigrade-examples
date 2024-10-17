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
    'num_steps' : 240,
    'finite_rise': 0,
    'macro_BC': 'slip',
    'mesh_copy_root': 'cylinder_5_5',
}

Brazilian_disk = {
    'height': 10.,
    'width': 36.,
    'chord': 22.,
    'app_rad': 12.5,
    'app_dep': 5.5,
    'spec_rad': 9.55,
    'spec_dep': 5.1,
    'tol': 0.001,
    'seed_size': 0.5,
    # material_E
    'lambda': 696.441593,
    'mu': 126.7138,
    'eta': -18.67498,
    'tau': -37.817315,
    'kappa': 15.177654,
    'nu': -24.071197,
    'sigma': -5.861821,
    'tau7': 792.523471,
    'cu0': 3.192202765,
    'fraction': 0.1,
    # macro_BC
    'macro_BC': 'brazil',
    'macro_disp': 1.0,
    'macro_duration': 1.0,
}