# shell_scripts
This folder containts the shell scripts used to run SCONE optimisations and evaluate par files for multiple altered models.

## Folder Content
+ **run_experiments_folder.sh:** This script can be used to optimise with SCONE all the files named **scone_main.scone** found recursively under `combinations_dir`. `combinations_dir` is a folder path entered as first argument, the results folders will be under SCONE_dir/results.\
Make sure to update `scone_install` accordingly.
To use this script, run the following terminal command:\
`run_experiments_folder.sh path/to/experiment/dir/to/optimise`
+ **evaluate_pars.sh:** This script can be used to evaluate with SCONE the parameter files with the best scores found under the simulation results folders under `results_dir`.
`results_dir` is a folder path entered as first argument. It needs to contain folders which are obtained from SCONE optimisations.\
Make sure to update `scone_install` accordingly.\
To use this script, run the following terminal command:\
`evaluate_pars.sh path/to/results/dir/to/evaluate`


