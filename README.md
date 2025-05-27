## ocrRBBR: an R package for the inference of Regression-Based Boolean Rules in multiomic data.
- This repository contains code and tutorials for executing ocrRBBR.
- Sample datasets to execute ocrRBBR are stored within the example_data directory.
- The ocrRBBR package supports parallel execution on multi-CPU platforms, enhancing accessibility for real-world Boolean rule inference applications.

## Step 1. ocrRBBR installation
The ocrRBBR codes are written in R version 4.1.3 and have been tested in both Windows and Linux environments. 

### Installation
1. Download the compiled package file `ocrRBBR_0.1.0.tar.gz` from this GitHub page.
2. Install the ocrRBBR package by running the following command in R:
   
```R
install.packages("path/to/RBBR_0.1.0.tar.gz", repos = NULL, type = "source")
```
<br>

### Dependencies  
Please ensure that you have the following packages installed. The glmnet package is required to fit ridge regressions. In order to run RBBR with parallel computing, the packages doParallel, foreach, and doSNOW need to be installed.

```R
install.packages("glmnet")
install.packages("doParallel")  
install.packages("foreach")
install.packages("doSNOW")
```  

<br>

## Step 2. Prepare input files
### Preprocessing input data
To preprocess raw data, including steps such as rescaling to bring each input feature within the [0,1] range, you can use the rbbr_scaling() function from the RBBR package.

```R
# Preprocessing input data
data_scaled   <- rbbr_scaling(data)
```

## Step 3. Train RBBR and use it for prediction purpose on new dataset
The RBBR package offers the rbbr_train() function for training the model on a dataset to extract Boolean rules, and the rbbr_predictor() function for utilizing the trained model to predict target values or labels on a new dataset.

### Train RBBR
For training the RBBR model on a dataset to extract Boolean rules, you can use the `rbbr_train()` function.

```R
# For training the RBBR model
trained_model <- rbbr_train(data, max_feature = NA, mode = NA, slope = NA, penalty = NA, weight_threshold = NA, balancing = NA, num_cores = NA)
```

```bash
# Required input arguments
# data              The dataset with rescaled features within the [0,1] interval.
#                   Each row represents a sample and each column represents a feature. The target variable (class) should be in the last column.  

# Optional input arguments  
# max_feature       The maximum number of input features allowed in a Boolean rule.
#                   The default value is 3.
 
# mode              Choose between "1L" for fitting 1-layered models or "2L" for fitting 2-layered models.
#                   The default value is "1L".
 
# slope             The slope parameter used in the Sigmoid activation function.
#                   The default value is 10.
 
# penalty           The penalty for the number of parameters in the BIC function.
#                   The default value is 1, but it can be adjusted (e.g., to 10) to favor simpler Boolean rules with fewer features.

# weight_threshold  Conjunctions with weights above this threshold in the fitted ridge regression models will be printed as active conjunctions in the output.
#                   The default value is 0.

# balancing         This is for adjusting the distribution of classes or categories within a dataset to ensure that each class is adequately represented.
#                   The default value is "True". Set it to "False", if you don't need to perform the data balancing.

# num_cores         Specify the number of parallel workers (adjust according to your system)
```
