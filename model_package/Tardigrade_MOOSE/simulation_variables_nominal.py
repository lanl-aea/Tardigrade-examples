import numpy
import copy

import model_package.Calibrate.calibration_tools as CT

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

elastic_cube = copy.deepcopy(elastic_cylinder)
del elastic_cube['diam']
del elastic_cube['mesh_copy_root']

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
    'app_dep': 6.0,
    'spec_rad': 9.55,
    'spec_dep': 5.1,
    'rad_diff': 2.95, #I'm not sure why app_rad - spec_rad is giving weird floating point values like 2.9499999993
    'tol': 0.0,
    'seed_size': 0.5,
    'platen_seed_size': 1.0,
    'mesh_copy_root': 'brazilian_disk_apparatus',
    # material_E
    'material_E': 833.75,
    'material_nu': 0.3625,
    'platen_E': 210000.,
    'platen_nu': 0.30,
    'cu0': 3.0,
    'fraction': 0.1,
    # macro parameters
    'macro_BC': 'brazil',
    'macro_disp': 0.5,
    'macro_duration': 1.0,
    'friction_coefficient': 0.1,
}

elastic_parameter_ordering = ['lambda', 'mu', 'eta', 'tau', 'kappa', 'nu', 'sigma',\
                              'tau1', 'tau2', 'tau3', 'tau4', 'tau5', 'tau6', 'tau7',\
                              'tau8', 'tau9', 'tau10', 'tau11']

elasticity_parameters = CT.Isbuga_micrormorphic_elasticity_parameters(Brazilian_disk['material_E'], Brazilian_disk['material_nu'], 0.5)
# Use a try/except because unittest.mock.Mock() has been problematic for calibration_tools.py
try:
    for key, param in zip(elastic_parameter_ordering, elasticity_parameters):
        Brazilian_disk[key] = param
except:
    pass

Brazilian_disk_platens = copy.deepcopy(Brazilian_disk)
