## Description
This folder contains the different python scripts used to process and generate desired **SCONE** files.

## Folder Content
### alice_thesis_scripts 
This folder contains the different scripts authored by Alice Bruel as part of her Master Thesis work. The methods are used by the scripts present in this folder to extract metrics from  **SCONE .sto** files.
### assess_results_folder.py
This script takes as input a path of a results folder with an organized folder structure. It evaluates different metrics for each of the parameters explored and generates a report folder containing visualization of these metrics. To use this script, run the following terminal command:\
`python3  assess_results_folder.py path/to/organized/results/folder`
### batch_evaluate_par.py
This script is called witha folder containing SCONE experiment result folders as first argument. It will identify all the **.par** files found recursively under the source folder and returns the paths of the par files which present the best score for each optimization. To use this script, run the following terminal command:\
`python3 batch_evaluate_par.py path/to/results/folder/to/evaluate`
### get_parameter_gain.py
This script takes as input a folder path of a SCONE optimization and returns the parameter value tested.Example:\
**input:** geyer_couplings_multiobjective_soleus-tib_ant_stance_KL_-7.f0914m.GH2010v8.S10W.D10.I\
**output:** -7\
To use this script, run the following terminal command:\
`python3 get_parameter_gain.py path/to/experiment/result/folder`
### par_sto_metrics.py
This file contains different methods used to extract and plot metrics of SCONE simulation results.
### par_to_controller.py
This file contains multiple methods aiming to parse SCONE .par files, to modify them and apply them on a controller file.
### rename_optim_folder_1d.py
This script takes as input a folder path of a SCONE optimization and renames it with the parameter gain value. Example:\
**input:** results/folder/path/geyer_couplings_multiobjective_soleus-tib_ant_stance_KL_-7.f0914m.GH2010v8.S10W.D10.I\
**output:** results/folder/path/-7\
To use this script, run the following terminal command:\
`python3 rename_optim_folder_1d.py path/to/1d/experiment/result/folder/to/rename`
### rename_optim_folder_2d.py
This script takes as input a folder path of a SCONE optimization and renames it with the parameter gain value. Example:\
**input:** results/folder/path/geyer_solLV_multiobjective_stance_KV_swing_KV_KV_0.0_KV_5.0.f0914m.GH2010v8.S10W.D10.I\
**output:** results/folder/path/0.0_5.0\
To use this script, run the following terminal command:\
`python3 rename_optim_folder_1d.py path/to/2d/experiment/result/folder/to/rename`
### scone_recursive_parser.py
This file contains multiple methods to handle SCONE controller files, namely parsing SCONE files to PANDAS Series, filtering, getter and setter methods, exporting to SCONE files, etc.
### using_recursive_parser_simple.py, using_recursive_parser_complex.py, using_recursive_parser_simple_cluster_dev.py, using_recursive_parser_simple_cluster_prod.py
These files are example files, showing how to parse a **.scone** controller file, import new data from a **.par** and generate the necessary folders to run experiments with desired parameter values. Code blocks are commented, displaying the functionalities available. To use these scripts, run the following terminal command:\
`python3 using_recursive_parser*.py`
