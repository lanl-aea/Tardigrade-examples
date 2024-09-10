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
    'calibration_increment': 4,
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
    'calibration_increment': 5,
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

# I41_damage
I41_damage = {
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
    'calibration_increment': 5,
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.2729,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
# Location of DNS files in petalibrary
I41_damage['DNS_files'] = [
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_fd20_damage_AT2_viscosity01/I41_02_sm3_fd20_l200_Ge01_Ee230_AT2_xi01_clampedfull_0.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_fd20_damage_AT2_viscosity01/I41_02_sm3_fd20_l200_Ge01_Ee230_AT2_xi01_clampedfull_300.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_fd20_damage_AT2_viscosity01/I41_02_sm3_fd20_l200_Ge01_Ee230_AT2_xi01_clampedfull_600.vtu',
    'RatelDNS/IDOX-ESTANE/I41_02/I41_02_sm3_fd20_damage_AT2_viscosity01/I41_02_sm3_fd20_l200_Ge01_Ee230_AT2_xi01_clampedfull_900.vtu',
    ]

# I41_damage
I43_damage = {
    # DNS parameters
    'diam': 5.05934739145766,
    'height': 4.400976867675782,
    'center': [3.0664815594943176, 3.065754248753459, 2.7645057067871095],
    'material_E': 450.0,
    'material_nu': 0.25,
    'cut': True,
    # Mesh file root to copy if Cubit is not found
    'mesh_copy_root': 'Ratel_I43_09',
    # parameters for micromorphic filter
    'acceleration': False,
    'velocity': False,
    'filter_parallel': 8,
    # parameters for calibration
    'calibration_case': 3,
    'calibration_increment': [1, 2, 3],
    'ignore_boundary': True,
    # paramters for Tardigrade-MOOSE
    'macro_disp': 0.156857,
    'macro_duration': 1.0,
    'macro_BC': 'clamp',
}
I43_damage_coarse = copy.deepcopy(I43_damage)
I43_damage['DNS_files'] = [
    'from_Erik_9-6-2024/Monitor_0.vtu',
    'from_Erik_9-6-2024/Monitor_15.vtu',
    'from_Erik_9-6-2024/Monitor_30.vtu',
    'from_Erik_9-6-2024/Monitor_45.vtu',
    'from_Erik_9-6-2024/Monitor_60.vtu',
    'from_Erik_9-6-2024/Monitor_70.vtu',
    'from_Erik_9-6-2024/Monitor_90.vtu',
    'from_Erik_9-6-2024/Monitor_100.vtu',
    'from_Erik_9-6-2024/Monitor_125.vtu',
    'from_Erik_9-6-2024/Monitor_150.vtu',
    'from_Erik_9-6-2024/Monitor_165.vtu',
    ]
I43_damage['DNS_forces'] = 'from_Erik_9-6-2024/Export_I43-09_QOI-DF05.csv'

I43_damage_coarse['DNS_files'] = [
    'from_Erik_9-10-2024/Monitor_0.vtu',
    'from_Erik_9-10-2024/Monitor_15.vtu',
    'from_Erik_9-10-2024/Monitor_30.vtu',
    'from_Erik_9-10-2024/Monitor_45.vtu',
    'from_Erik_9-10-2024/Monitor_65.vtu',
    'from_Erik_9-10-2024/Monitor_75.vtu',
    'from_Erik_9-10-2024/Monitor_85.vtu',
    'from_Erik_9-10-2024/Monitor_100.vtu',
    'from_Erik_9-10-2024/Monitor_125.vtu',
    'from_Erik_9-10-2024/Monitor_150.vtu',
    'from_Erik_9-10-2024/Monitor_180.vtu',
    'from_Erik_9-10-2024/Monitor_215.vtu',
    ]
I43_damage_coarse['DNS_forces'] = 'from_Erik_9-10-2024/Export_I43-09_QOI-DF20.csv'