## Description
This folder containts multiple shell scripts used to transfer **SCONE** files via ssh from and to the computer cluster, deploy jobs on the cluster and later evaluate par files.

## Folder Content
### create_folder_structure_1d.sh
This script can be used to generate proper folder structure from SCONE results folder. First argument is SCONE results folder. Second argument is destination folder. Make sure that it ends with '1d'. Make sure to modify python_script according to where the scripts will be run. To use this script, run the following terminal command:\
`create_folder_structure_1d.sh path/to/unorganized/results/folder/1d desired/path/to/organized/results/folder/1d`
### create_folder_structure_2d.sh
This script can be used to generate proper folder structure from SCONE results folder. First argument is SCONE results folder. Second argument is destination folder. Make sure that it ends with '2d'. Make sure to modify python_script according to where the scripts will be run. To use this script, run the following terminal command:\
`create_folder_structure_2d.sh path/to/unorganized/results/folder/2d desired/path/to/organized/results/folder/2d`
### cluster_run_1d.sbatch
This script is used to deploy SCONE optimization folders for 1d experiments on the Biorob computer cluster.\
Make sure to update the directory definitions in the script according to your personal files.\
To deploy this file on a production node, run: `sbatch -p prod cluster_run_1d.sbatch`\
To deploy this file on a development node, run: `sbatch cluster_run_1d.sbatch` (make sure to update the number of cores and memory pool.)
### cluster_run_2d.sbatch
This script is used to deploy SCONE optimization folders for 2d experiments on the Biorob computer cluster.\
Make sure to update the directory definitions in the script according to your personal files.\
To deploy this file on a production node, run: `sbatch -p prod cluster_run_2d.sbatch`\
To deploy this file on a development node, run: `sbatch cluster_run_2d.sbatch` (make sure to update the number of cores and memory pool.)
### evaluate_pars.sh
This script can be used to evaluate with SCONE the parameter files with the best scores found under the simulation results folders under `results_dir`.
`results_dir` is a folder path entered as first argument. It needs to contain folders which are obtained from SCONE optimizations.\
Make sure to update `scone_install` accordingly.\
To use this script, run the following terminal command:\
`evaluate_pars.sh path/to/results/dir/to/evaluate`
### run_experiments_folder.sh
This script can be used to optimize with SCONE all the files named **scone_main.scone** found recursively under `combinations_dir`. `combinations_dir` is a folder path entered as first argument.\
Make sure to update `scone_install` accordingly.
To use this script, run the following terminal command:\
`run_experiments_folder.sh path/to/experiment/dir/to/optimize`
### rsync_results_folder.txt
This script can be used to synchronize two folders via ssh. Modified and new files on the first folder will be updated on the second folder.
### scp_exports_folder.txt
This script can be used to secure copy via ssh an experiment folder from personal computer to a linux folder.
