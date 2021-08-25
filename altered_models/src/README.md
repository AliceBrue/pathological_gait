# src
This folder contains the python scripts used to process and generate SCONE files for multiple altered models, and to plot optimisation results.

## Folder Content
+ **assess_results_folder.py:** This script contains multiple functions to evaluate various metrics for each of the parameters explored and generates a report folder containing plots of these metrics.
+ **batch_evaluate_par.py:** This script is called with a folder containing SCONE experiment result folders as first argument. It will identify all the **.par** files found recursively under the source folder and returns the paths of the par files which present the best score for each optimisation. 
+ **extract_sto.py:** This script contains functions to extract variables from SCONE **.sto** files and average them over gait phases.
+ **par_sto_metrics.py:** This script contains various methods used to extract and plot metrics of SCONE simulation results.
+ **par_to_controller.py:** This script contains multiple methods aiming to parse SCONE .par files, to modify them and apply them on a controller file.
+ **plot_results_folder.py:** This script takes as input a path of a results folder with an organised folder structure to generate a report folder containing plots of various metrics for each of the parameters explored.
+ **run_recursive_parser.py:**  This script is an example showing how to parse a **.scone** controller file, generate scone folders with multiple combination of controller files with altered parameter values and other parameter values imported from the healthy **.par** file. Code blocks are commented, displaying the functionalities available.
+ **scone_recursive_parser.py:** This script contains multiple methods to handle SCONE controller files, namely parsing SCONE files to PANDAS Series, filtering, getter and setter methods, exporting to SCONE files, etc.
+ **spasticity_index.py.py:** This script contains functions to compute spasticity indexes.