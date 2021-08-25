"""
This script is called with a folder containing SCONE experiment result folders
as first argument. It will identify all the .par files found recursively under
the source folder and returns the paths of the par files which present the best
score for each optimization.
"""
import sys, os
from pathlib import Path
import numpy as np


def str_to_number(s):
    """Returns number  when a string is a number otherwise returns inf
    Parameters
    ---------
    s: (string)

    Returns
    ---------
    float or inf
    """
    try:
        value = float(s)
        if value <= 0:
            return float('inf')
        else:
            return value
    except ValueError:
        return float('inf')

par_files_to_evaluate = []
if __name__ == "__main__":
    folder_path = sys.argv[1] # path of folder containing optim_results
    # extract list of folders from folder_path
    optim_folders = [ os.path.join(folder_path,name) for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name)) ]
    # loop through all optim folders
    for optim_folder in optim_folders:
        # get all par files
        par_files =  [str(par_file) for par_file in Path(optim_folder).rglob('*.par')]
        # sort to get higher generations first
        par_files.sort(reverse=True)
        # extract scores  from par file names
        scores_str = [par_file.rstrip('.par').split('_')[-1] for par_file in par_files]
        scores = [str_to_number(score_str) for score_str in scores_str]
        # get files with best scores
        best_score_index = np.argmin(scores)
        best_par_file = par_files[best_score_index]
        # only add par files who have not been evaluated yet
        if not os.path.isfile(best_par_file+'.sto'):
            par_files_to_evaluate.append(best_par_file)
        
    print('\n'.join(par_files_to_evaluate))