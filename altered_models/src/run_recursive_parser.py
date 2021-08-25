"""
 This script is an example showing how to parse a **.scone** controller file,
 generate scone folders with multiple combination of controller files with
 altered parameter values and other parameter values imported from the healthy **.par** file.
 Code blocks are commented, displaying the functionalities available.
"""
import scone_recursive_parser as srp
import par_to_controller as ptc

# import controller file to string
controller = srp.controller_file_to_string("../models/gait_controller/reflex_controller.txt")

# convert controller string to recursive dictionnary
dict = srp.dict_recursive_split(controller)

# convert controller string to recursive pandas Series
series = srp.series_recursive_split(controller)

# apply optimised parameter values to controller with new value ranges.
par_file = '../models/scone_base_files/optimisation.par'
min_max_ratio = 0.05  # new max/min = value*(1 +/- min_max_ratio)
std_ratio = 1  # new std = value*std_ratio

series_new_range = ptc.apply_par_file_to_series_with_new_range_2(par_file, series, min_max_ratio)

# generate scone folders with multiple combination of controller files.
case = '1D'  # '1D' or '2D'
# Altered parameter
phase_1 = 'stance'  # 'stance' or 'swing'
param_1 = 'KS'  # 'KS', 'KL', 'KF', 'TA_KS', 'weak_KS', 'weak_KL' or 'weak_KF'
key_values_1 = range(150, 400, 50)
# Second altered parameter if 2D
phase_2 = 'swing'   # 'stance' or 'swing'
param_2 = 'KS'  # 'KS', 'KL', 'KF', 'TA_KS', 'weak_KS', 'weak_KL' or 'weak_KF'
key_values_2 = range(150, 400, 50)

if phase_1 == 'stance':
    states_1 = 'LateStance EarlyStance Liftoff'
elif phase_1 == 'swing':
    states_1 = 'Swing Landing'
if param_1 == 'TA_KS':
    source_1 = 'tib_ant'
else:
    source_1 = None
key_dict_0 = {'key_name': param_1.split('_')[-1], 'states': states_1, 'target': ['soleus', 'gastroc'],
              'source': source_1, 'key_values': key_values_1}

if case == '1D':
    key_dict_list = [key_dict_0]
    output_folder_path = '../exports/1D/' + phase_1 + '/' + 'stance_' + param_1

elif case == '2D':
    if phase_2 == 'stance':
        states_2 = 'LateStance EarlyStance Liftoff'
    elif phase_2 == 'swing':
        states_2 = 'Swing Landing'
    if param_2 == 'TA_KS':
        source_2 = 'tib_ant'
    else:
        source_2 = None
    key_dict_1 = {'key_name': param_1.split('_')[-1], 'states': states_2, 'target': ['soleus', 'gastroc'],
                  "source": source_2, 'key_values': key_values_2}
    key_dict_list = [key_dict_0, key_dict_1]
    output_folder_path = '../exports/2D/' + phase_1 + '_' + param_1 + '_' + phase_2 + '_' + param_2

scone_base_files_path = '../models/scone_base_files'
scone_main_file_name = 'scone_main.scone'

srp.generate_controllers_with_key_combinations(series_new_range, key_dict_list, output_folder_path,
                                                 scone_base_files_path, scone_main_file_name)
