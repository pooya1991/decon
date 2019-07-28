## Concise Siren pipeline

This repo is a concise version of Siren. Specifically, all the files in the C++ folder, the java folder and the scripts folder are thrown away. Also some files in the python folder that were not used in the pipeline are thrown away as well. Instead some files are added to this repo that are not present in the original Siren repo. Specifically there are some R scripts, a sample data and the results that siren generates for that data and some compiled C++ codes. Below is the list of all files and their description:
* **15b.ms1:** the sample data
* **ComputeMS1:** a Linux executable file which is obtained by compiling the contents of Mercury8 folder.
* **ISOTOPE.DAT:** this file is used by ComputeMS1 internally
* **ratest_batch:** this is a Linux executable file which is obtained by compiling the contents of the C++ folder which were related to regression. This file is used for the OLS regression step in Siren
* The python scripts are the same as python scripts in the original Siren repo. There are two main differences: 1) in comupute_paths_job.py which does the Lasso step of Siren, compute_paths.py is evoked as a system call, but I modified compute_paths_job.py to import the functions in compute_paths.py. 2) In regression_justms1_sparse_job_sgd.py, instead of using a shell command to do the regression, I directly called ratest_batch with the arguments that were stored in the shell script
* There are R scripts in this repo which have python counterparts. They do the same thing that the corresponding python script does. Since these files are not well tested, it's better to work with the python scripts instead. There's only one important feature that one has to be aware of. **analysis.R** script does the step of Siren concerning generating X and Y matrices. At first I intended to write my R version of Siren to make it more functional. But later I abandoned the idea due to lack of time.
* **matrix_operations.R:** include functions for working with sparse matrices defined in Siren
* **regression_analysis.R:** this scripts intends to substitute the OLS step of siren with Lasso regression per scan
* **15b_ms1_sparsebinned_singlescans_tannotated:** this folder contains the results concerning making X and Y matrices, decoy X matrices, annotation files and joint annotation file
* **15b_peptide_abundances:** this folder contains the results regarding combining OLS regression coefficients, building a combined B matrix, smoothing and Lasso regression
* **15b_precursorprofiles.txt:** this is the final elution profiles extracted by Siren
* **15b_scannums.txt:** real scan numbers