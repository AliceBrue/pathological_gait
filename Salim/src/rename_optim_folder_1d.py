import sys, os
from pathlib import Path
import numpy as np

# change this if suffix of folder changes
suffix = '.f0914m.GH2010v8.S10W.D10.I'

'''
This script takes as input a folder path of a SCONE optimization and renames it with the parameter gain value.
example input: /results/folder/path/geyer_couplings_multiobjective_soleus-tib_ant_stance_KL_-7.f0914m.GH2010v8.S10W.D10.I
        output: /results/folder/path/-7
To use this script, run the following terminal command:\
    python3 rename_optim_folder_1d.py *path/to/1d/experiment/result/folder/to/rename*
'''
if __name__ == "__main__":
    folder_path = sys.argv[1] # path of folder containing optim_results
    
    folder_path_no_suffix = folder_path.replace(suffix,'') # remove suffix
    parameter_gain = folder_path_no_suffix.rsplit(sep='_',maxsplit=1)[-1] # extract score
    parent_folder = folder_path_no_suffix.rsplit(sep='/',maxsplit=1)[0] # extract score
    
    os.rename(folder_path,os.path.join(parent_folder, parameter_gain))
    