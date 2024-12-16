import copy

elastic_cylinder = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'seed': 0.25,
    'material_E': 165.0,
    'material_nu': 0.39,
    'material_rho': 2.0e-9,
    'disp': 0.01,
    'num_steps': 5,
    'top_surface_id': 1,
    'bottom_surface_id': 2,
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
    'calibration_increment': [1, 2, 3, 4, 5],
    'UQ_file': False,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'slip',
}

elastic_cylinder_clamp = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'seed': 0.25,
    'material_E': 165.0,
    'material_nu': 0.39,
    'material_rho': 2.0e-9,
    'disp': 0.01,
    'num_steps': 5,
    'top_surface_id': 1,
    'bottom_surface_id': 2,
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
    'calibration_increment': [1, 2, 3, 4, 5],
    'UQ_file': False,
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}

F83 = {
    # DNS parameters
    'diam': 5.0,
    'height': 5.0,
    'material_E': 23520,
    'material_nu': 0.25,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'Ratel_F83',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 1,
    'calibration_increment': [4],
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.05,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
# Location of DNS files in petalibrary
F83['DNS_files'] = [
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_0.vtu',
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_1.vtu',
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_2.vtu',
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_3.vtu',
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_4.vtu',
    'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/F83_01_er10_sm10_nt50_nsf_5.vtu',
    ]
# DNS force file
F83['DNS_forces'] = 'RatelDNS/SAND-KELF/F83_01_2strain_clamped_linear_el_er10_sm10_nt50/'

I41_02 = {
    # DNS parameters
    'diam': 6.0,
    'height': 5.485,
    'material_E': 450.0,
    'material_nu': 0.25,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'Ratel_I41_02',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 3,
    'calibration_increment': [5],
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.2729,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
# Location of DNS files in petalibrary
I41_02['DNS_files'] = [
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_0.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_20.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_40.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_60.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_80.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_100.vtu',
    ]
# DNS force file
I41_02['DNS_forces'] = 'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3.csv'

# I41 RVEs
I41_RVEs = copy.deepcopy(I41_02)
additional_files = [
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_10.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_30.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_50.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_70.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_df20_clamped_elast3/I41_02_sm3_fd20_l200_Ge01_Ee250_clamped_elast3_90.vtu',
    ]

# softening study
I41_02_clamp = {
    # DNS parameters
    'diam': 6.0,
    'height': 5.485,
    'material_E': 450.0,
    'material_nu': 0.25,
    'cut': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.30,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
    'lambda': 696.441593,
    'mu': 126.7138,
    'eta': -18.67498,
    'tau': -37.817315,
    'kappa': 15.177654,
    'nu': -24.071197,
    'sigma': -5.861821,
    'tau7': 792.523471,
    'cchi0': 3.192202765,
}

# Combine file lists, force the order
all_files = []
for file1, file2 in zip(I41_02['DNS_files'], additional_files):
    all_files.append(file1)
    all_files.append(file2)
all_files.append(I41_02['DNS_files'][-1])

# I41_RVEs
I41_RVEs['DNS_files'] = all_files
I41_RVEs['rho_binder'] = 1.910e-9
I41_RVEs['rho_grain'] = 2.001e-9

# I43_damage
I43_damage_coarse_finetime = {
    # DNS parameters
    'diam': 5.059347391457650,
    'height': 4.400976867675780,
    'center': [2.529673695728830, 2.529673695728830, 2.200488433837890],
    'material_E': 250.0,
    'material_nu': 0.25,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'Ratel_I43_09',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for elastic calibration
    'plastic_calibration_case': 6,
    'calibration_increment_elastic': [0, 1, 2, 3, 4],
    'calibration_increment_plastic': [3, 5, 7, 9, 10, 11, 12, 13],
    # parameters for plastic calibration
    'cohesion_case': 1,
    'cohesion_increment': [3, 4, 5],
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.210*0.72221649302258428, #timestep for Monitor_000124.vtu
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
I43_damage_coarse_finetime['DNS_files'] = [
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000000.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000005.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000010.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000015.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000020.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000025.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000030.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000035.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000040.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000045.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000050.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000055.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000060.vtu',
    'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Monitor_000065.vtu',
    ]
I43_damage_coarse_finetime['DNS_forces'] = 'RatelDNS/IDOX-ESTANE/I43_09/BaseLine-BIN4_SM10_DF20_MI10_MD0.5/Forces_OUT_I43-09_BIN4_SM10_DF20_MI10_MD0.5.csv'
