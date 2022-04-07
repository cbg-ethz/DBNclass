Simulated data:

- script dbnsimfns.R contains functions needed for generating a random DBN structure, parameters and data
- script DBNCVparopt.R contains functions needed for computing MAE
- script DBNsimulation.R includes the function needed for one simulation replicate DBNsimulation() and the code for parallelization of 50 replicates
- the results of simulation can be found in the file “DBNsimres/DBNsimres_4slices_30samples_100n.rds” 

GSE5462, breast cancer dataset:

- raw data can be found here
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5462
- please note, that the data cannot be published in the GitHub folder, so scripts
show only major steps, but the data needs to be downloaded from the GEO server
- script GSE5462RMAnorm.R contains the code for rma normalization 
- script Gene_filtering_GSE5462.R contains script with major steps of gene filtering
- script GSE5462struct.R  demonstrates how to construct blacklists and penalization matrices, folder data contains .rds files with pre-constructed matrices
- script CV_GSE5462_parallel.R contains a function that runs one CV iteration and the code for parallel runs of all 52 samples, it requires scripts DBNCVopt.R, DBNCVparopt.R and DBNpreprocopt.R 
- script DBN_learn_best_GSE5462.R  contains the learning of the final network that was used in the downstream analysis
- the results of CV are contained in the file biological_data_results/all_models_CV_MAE_GSE5462.rds     
- file DBNotheralgos.R contains the script for classification with random forest and naïve Bayes

GSE37182, colon cancer dataset:

- raw data can be found here
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5462
- please note, that the data cannot be published in the GitHub folder, so scripts
show only major steps, but the data needs to be downloaded from the GEO server
- script GSE37182norm.R contains the code with normalization code
- script Gene_filtering_GSE37182.R contains script with major steps of gene filtering
- script GSE37182struct.R  demonstrates how to construct blacklists and penalization matrices, folder data contains .rds files with pre-constructed matrices
- script CV_GSE37182_parallel.R  contains a function that runs one CV iteration and the code for parallel runs of all 29 iterations, it requires scripts DBNCVcolonMAEopt.R, DBNCVparopt.R and DBNpreprocopt.R 
- script DBN_learn_best_GSE37182.R  contains the learning of the final network that was used in the the downstream analysis.
- the results of CV are contained in the file
biological_data_results/all_models_CV_MAE_37182.rds       


