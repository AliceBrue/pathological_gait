"""
This script takes as input a path of a results folder with an organised folder structure to
generate a report folder containing plots of various metrics for each of the parameters explored.
"""
import os
from pathlib import Path
import assess_results_folder as arf

# Set this to '1D' or '2D' to analyze 1D or 2D experiment results
case = '2D'
side = 'l'  # 'l' or 'r' for left or right leg plots

results_folder = "../results/"
results_folder_1d = os.path.join(results_folder, '1D')
results_folder_2d = os.path.join(results_folder, '2D')

# extract folders for each 1d optimization
if case == '1D':
    sto_files_1d = Path(results_folder_1d).rglob('*.par.sto')
    parameter_folders_1d = [sto_file.parent.parent for sto_file in sto_files_1d]
    parameter_folders_1d = list(set(parameter_folders_1d))
    for parameter_folder in parameter_folders_1d:
        arf.assess_parameter_folder_1d(parameter_folder, side)

# extract folders for each 2d optimization
elif case == '2D':
    sto_files_2d = Path(results_folder_2d).rglob('*.par.sto')
    parameter_folders_2d = [sto_file.parent.parent for sto_file in sto_files_2d]
    parameter_folders_2d = list(set(parameter_folders_2d))
    for parameter_folder in parameter_folders_2d:
        arf.assess_parameter_folder_2d(parameter_folder, side)

