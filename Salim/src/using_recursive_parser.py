import scone_recursive_parser as srp
import par_to_controller as ptc
import numpy as np

# import controller file to string
controller = srp.controller_file_to_string("../models/gait_controllers/alice_spindle_add_recip_controller.txt")

# convert controller string to recursive dictionnary
dict = srp.dict_recursive_split(controller)

# convert controller string to recursive pandas Series
series = srp.series_recursive_split(controller)

# apply par file to controller with new value ranges.
par_file = '../models/alice_spindle_add_recip_base_files/add_recip_optim.par'
min_max_ratio = 0.01  # new max/min = value*(1 +/- min_max_ratio)
std_ratio = 1  # new std = value*std_ratio

series_new_range = ptc.apply_par_file_to_series_with_new_range_2(par_file, series, min_max_ratio)  # same std
#series_fixed_values = ptc.apply_par_file_to_series_with_fixed_values(par_file, series)

# generate scone folders with multiple combination of controller files.
key_dict_0 = {'key_name': 'KS', 'states': 'LateStance EarlyStance Liftoff', 'target': ['soleus', 'gastroc'],
              'source': '', 'key_values': range(200, 1000, 100)}

key_dict_1 = {'key_name': 'KS', 'states': 'Swing Landing', 'target': ['soleus', 'gastroc'], 'source': '',
              'key_values': range(200, 1000, 100)}

key_dict_list = [key_dict_0, key_dict_1] # [key_dict_0] for 1D, [key_dict_0, key_dict_1] for 2D

output_folder_path = '../exports/alice_spindle_add_recip_sol_gas_KS_2D'
scone_base_files_path = '../models/alice_spindle_add_recip_base_files'
scone_main_file_name = 'scone_main.scone'

srp.generate_controllers_with_key_combinations_2(series_new_range, key_dict_list,
                                               output_folder_path, scone_base_files_path, 
                                               scone_main_file_name)
