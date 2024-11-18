import copy

elastic_cylinder = {
    # DNS parameters
    'diam': 2.0,
    'height': 2.0,
    'material_E': 450.0,
    'material_nu': 0.25,
    'rho_binder': 1.860e-9,
    'rho_grain': 2.001e-9,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'GEOS_cylinder',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 3,
    'calibration_increment': 5,
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.2729,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
# Location of DNS files in petalibrary
elastic_cylinder['DNS_file'] = 'vtkOutput.pvd'
elastic_cylinder['DNS_fileroot'] = '/projects/tea/PSAAP/GEOS_new_cylinder_from_Jay/'
elastic_cylinder['DNS_forces'] = '.'