import sys, os
from pathlib import Path
import numpy as np

# change this if suffix of folder changes
suffix = '.f0914m.GH2010v8.S10W.D10.I'

'''
This script takes as input a folder path of a SCONE optimization and returns the parameter value tested.
example input: geyer_couplings_multiobjective_soleus-tib_ant_stance_KL_-7.f0914m.GH2010v8.S10W.D10.I
        output: -7
To use this script, run the following terminal command:\
    python3 get_parameter_gain.py *path/to/experiment/result/folder*
'''
if __name__ == "__main__":
    folder_path = sys.argv[1] # path of folder containing optim_results
    folder_path_no_suffix = folder_path.replace(suffix,'') # remove suffix
    parameter_gain = folder_path_no_suffix.rsplit(sep='_',maxsplit=1)[-1] # extract score
    print(parameter_gain)
    