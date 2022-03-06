# altered_models
This repository contains the shell scripts and python scripts used to generate and optimise SCONE files for multiple altered models, and plot optimisation results.
Schematically:
1. Experiment folders with SCONE files used to run optimisation are generated for each altered model of interest.
2. SCONE optimisation is run for each altered model.
3. SCONE evaluation of the best result is run for each altered model.
4. Various metrics are plotted for each altered model.

## Folder Content
+ **exports:** This folder contains structured experiment folders with SCONE files used to run optimisation for each altered model of interest.
	       It has the following structure: exports/1D/biomechanical/ biomechanical parameters folders
							  stance/ stance reflex gains folders
							  swing/ swing reflex gains folders
						       2D/ 2D parameters folders
+ **models:** This folder contains the SCONE files used to generate the previous experiment folders.
+ **results:** This folder contains structured folders with optimisation results of each altered models of interest. It structure is the same as the exports folder.
+ **shell_scripts:** This folder containts the shell scripts used to run SCONE optimisations and evaluate par files for multiple altered models.
+ **src:** This folder contains the python scripts used to process and generate SCONE files for multiple altered models, and to plot optimisation results.
