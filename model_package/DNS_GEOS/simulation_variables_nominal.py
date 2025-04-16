import copy

elastic_cylinder = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'material_E': 250.,
    'material_nu': 0.2,
    'rho_binder': 1.935e-9,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'GEOS_cylinder',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 3,
    'calibration_increment': [5],
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.1250143659720484,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
}
# Handle DNS files for petalibrary copy
root = 'GEOS-MPM/MiscellaneousSimulations/UniaxialCylinderCompression'
levels = ['000000', '000023', '000045', '000067', '000089', '000111', '000133', '000155',
          '000177', '000199', '000221', '000243', '000265', '000266']
ranks = 8
elastic_cylinder['DNS_files'] = {
    'main': [f'{root}/vtkOutput.pvd'] + [f'{root}/reactionHistory.csv'],
    'vtms': [f'{root}/vtkOutput/{level}.vtm' for level in levels]}
for level in levels:
    elastic_cylinder['DNS_files'][level] = [f'{root}/vtkOutput/{level}/particles/Level0/ParticleRegion1/rank_{i}.vtu' for i in range(0, ranks)]

# After files have been copied to peta_data_copy/, only the following items are needed to run the workflow(s)
elastic_cylinder['DNS_file'] = 'vtkOutput.pvd'
elastic_cylinder['DNS_fileroot'] = 'GEOS_elastic_cylinder'
elastic_cylinder['DNS_forces'] = 'reactionHistory.csv'

I43_01 = {
    # DNS parameters
    'diam': 5.0501088,
    'height': 4.6736,
    'material_E': 250.,
    'material_nu': 0.2,
    'rho_binder': 1.935e-9,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'GEOS_cylinder',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 1,
    # parameters for calibration
    'calibration_case': 3,
    'calibration_increment': [5],
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.1250143659720484,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
}

I43_01['DNS_fileroot'] = 'GEOS_2025-03-02_I43-01_take2'
I43_01['DNS_file'] = 'vtkOutput_short.pvd'
I43_01['DNS_forces'] = 'reactionHistory.csv'