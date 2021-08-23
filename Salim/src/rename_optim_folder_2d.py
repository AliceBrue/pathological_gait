import sys, os
from pathlib import Path
import numpy as np

# change this if suffix of folder changes
suffix = '.f0914m.GH2010v8.S10W.D10.I'

'''
This script takes as input a folder path of a 2d SCONE optimization and renames it with the parameter gain value.
example input: /results/folder/path/geyer_solLV_multiobjective_stance_KV_swing_KV_KV_0.0_KV_5.0.f0914m.GH2010v8.S10W.D10.I
        output: /results/folder/path/0.0_5.0
To use this script, run the following terminal command:\
    python3 rename_optim_folder_1d.py *path/to/2d/experiment/result/folder/to/rename*
'''
if __name__ == "__main__":
    folder_path = sys.argv[1] # path of folder containing optim_results
    
    folder_path_no_suffix = folder_path.replace(suffix,'') # remove suffix
    first_parameter_gain = folder_path_no_suffix.split(sep='_')[-3]
    second_parameter_gain = folder_path_no_suffix.split(sep='_')[-1]

    parent_folder = folder_path_no_suffix.rsplit(sep='/',maxsplit=1)[0] # extract score
    
    os.rename(folder_path,os.path.join(parent_folder, '_'.join([first_parameter_gain,second_parameter_gain])))	
    