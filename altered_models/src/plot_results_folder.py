"""
This script takes as input a path of a results folder and generates a report folder
containing plots of various metrics for each of the explored parameters.
The results folder must have the same struture as the exports folder:
results/1D/biomechanical/ biomechanical parameters folders
           /stance/ stance reflex gains folders
           /swing/ swing reflex gains folders
        /2D/ 2D parameters folders
"""
import os
from pathlib import Path
import assess_results_folder as arf


# Set to plot 1D, 2D or repeated results
results_folder = "../results"  # to adapt
case = '2D'  # set to '1D', '2D' ro 'repeated' to analyze 1D, 2D or repeated experiments results
side = 'l'  # 'l' or 'r' for left or right leg plots
scone_folder = '.f0914m.GH2010v8.S05W.D15.I'  # depends on SCONE installation

# extract folders for each 1d optimization
if case == '1D':
    results_folder_1d = os.path.join(results_folder+'/', '1D')
    sto_files_1d = Path(results_folder_1d).rglob('*.par.sto')
    parameter_folders_1d = [sto_file.parent.parent for sto_file in sto_files_1d]
    parameter_folders_1d = list(set(parameter_folders_1d))
    for parameter_folder in parameter_folders_1d:
        arf.assess_parameter_folder_1d(parameter_folder, side, scone_folder)

# extract folders for each 2d optimization
elif case == '2D':
    results_folder_2d = os.path.join(results_folder+'/', '2D')
    sto_files_2d = Path(results_folder_2d).rglob('*.par.sto')
    parameter_folders_2d = [sto_file.parent.parent for sto_file in sto_files_2d]
    parameter_folders_2d = list(set(parameter_folders_2d))
    for parameter_folder in parameter_folders_2d:
        arf.assess_parameter_folder_2d(parameter_folder, side, scone_folder)

elif case == 'repeated':
    results_folder_rep = os.path.join(results_folder+'/', 'repeated')
    sto_files_rep = Path(results_folder_rep).rglob('*.par.sto')
    parameter_folders_rep = [sto_file.parent.parent for sto_file in sto_files_rep]
    parameter_folders_rep = list(set(parameter_folders_rep))
    for parameter_folder in parameter_folders_rep:
        arf.assess_parameter_folder_rep(parameter_folder, side, scone_folder)
