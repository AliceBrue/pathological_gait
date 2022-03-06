"""
 This script generates scone folders with multiple combination of controller files
 with altered parameter values and the remaining parameter values imported from the healthy **.par** file.
 It does so by parsing the **.scone** controller file.
 The exports folder where are generated the scone folder has the following structure:
 exports/1D/biomechanical/ biomechanical parameters folders
           /stance/ stance reflex gains folders
           /swing/ swing reflex gains folders
        /2D/ 2D parameters folders
"""
import scone_recursive_parser as srp
import par_to_controller as ptc
import os


# Set to generate scone folders for 1D or 2D parameters
export_folder = '../exports'
case = '1D'  # '1D' or '2D'
# Altered parameter
param_1 = 'KS'  # 'max_isometric_force', 'optimal_fiber_length', 'KS', 'KL', 'KF', 'TA_KS', 'weak_KS',
# 'weak_KL', or 'weak_KF'
phase_1 = 'stance'  # 'stance' or 'swing'
key_values_1 = range(150, 450, 50)
target_1 = ['soleus', 'gastroc']
# Second altered parameter if 2D
param_2 = 'KS'  # 'KS', 'KL', 'KF', 'TA_KS', 'weak_KS', 'weak_KL' or 'weak_KF'
phase_2 = 'swing'  # 'stance' or 'swing'
key_values_2 = range(150, 500, 100)
target_2 = ['soleus', 'gastroc']

# import controller file to string
controller = srp.controller_file_to_string("../models/gait_controller/reflex_controller.txt")

# convert controller string to recursive dictionnary
dict = srp.dict_recursive_split(controller)

# convert controller string to recursive pandas Series
series = srp.series_recursive_split(controller)

# apply optimised parameter values to controller with 5% varation constraint imposed by new value ranges.
par_file = '../models/scone_base_files/optimisation.par'
min_max_ratio = 0.05  # 5% variation constraint imposed by new max/min = value*(1 +/- min_max_ratio)
series_new_range = ptc.apply_par_file_to_series_with_new_range(par_file, series, min_max_ratio)

# create folder structure
if not os.path.isdir(export_folder):
    os.makedirs(export_folder)
if not os.path.isdir(export_folder+'/1D'):
    os.makedirs(export_folder+'/1D')
if not os.path.isdir(export_folder+'/1D/biomechanical'):
    os.makedirs(export_folder+'/1D/biomechanical')
if not os.path.isdir(export_folder+'/1D/stance'):
    os.makedirs(export_folder+'/1D/stance')
if not os.path.isdir(export_folder+'/1D/swing'):
    os.makedirs(export_folder+'/1D/swing')
if not os.path.isdir(export_folder+'/2D'):
    os.makedirs(export_folder+'/2D')

# generate controllers with altered reflex parameters
if param_1 not in ['max_isometric_force', 'optimal_fiber_length']:
    if phase_1 == 'stance':
        states_1 = 'LateStance EarlyStance Liftoff'
    elif phase_1 == 'swing':
        states_1 = 'Swing Landing'
    if param_1 == 'TA_KS':
        source_1 = 'tib_ant'
    else:
        source_1 = None
    key_dict_1 = {'key_name': param_1.split('_')[-1], 'states': states_1, 'target': target_1,
                  'source': source_1, 'key_values': key_values_1}

    if case == '1D':
        key_dict_list = [key_dict_1]
        output_folder_path = export_folder + '/1D/' + phase_1 + '/' + phase_1 + '_' + param_1

    elif case == '2D':
        if phase_2 == 'stance':
            states_2 = 'LateStance EarlyStance Liftoff'
        elif phase_2 == 'swing':
            states_2 = 'Swing Landing'
        if param_2 == 'TA_KS':
            source_2 = 'tib_ant'
        else:
            source_2 = None
        key_dict_2 = {'key_name': param_2.split('_')[-1], 'states': states_2, 'target': ['soleus', 'gastroc'],
                      "source": source_2, 'key_values': key_values_2}
        key_dict_list = [key_dict_1, key_dict_2]
        output_folder_path = '../new_exports/2D/' + phase_1 + '_' + param_1 + '_' + phase_2 + '_' + param_2

    scone_base_files_path = '../models/scone_base_files'
    scone_main_file_name = 'main_scone.scone'

    srp.generate_controllers_with_key_combinations(series_new_range, key_dict_list, output_folder_path,
                                                     scone_base_files_path, scone_main_file_name)

# generate scone main files with altered biomechanical parameters
else:
    key_dict_1 = {'key_name': param_1, 'target': target_1, 'key_values': key_values_1}

    if case == '1D':
        key_dict_list = [key_dict_1]
        output_folder_path = export_folder + '/1D/biomechanical/' + param_1

    elif case == '2D':
        key_dict_2 = {'key_name': param_2, 'target': target_2, 'key_values': key_values_2}
        key_dict_list = [key_dict_1, key_dict_2]
        output_folder_path = '../new_exports/2D/' + param_1 + '_' + param_2

    scone_base_files_path = '../models/scone_base_files'
    scone_main_file_name = 'main_scone.scone'

    srp.generate_scone_main_files(series_new_range, key_dict_list, output_folder_path,
                                                   scone_base_files_path, scone_main_file_name)

