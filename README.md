## ocrRBBR: an R package for the inference of Regression-Based Boolean Rules in multiomic data.
- This repository contains code and tutorials for executing ocrRBBR.
- Sample datasets to execute ocrRBBR are stored within the Data directory.
- The ocrRBBR package supports parallel execution on multi-CPU platforms, enhancing accessibility for real-world Boolean rule inference applications.

### Step 1. ocrRBBR installation
The ocrRBBR codes are written in R version 4.1.3 and have been tested in both Windows and Linux environments. 

### Installation
1. Download the compiled package file `RBBR_0.1.0.tar.gz` from this GitHub page.
2. Install the ocrRBBR package by running the following command in R:
   
```R
install.packages("path/to/RBBR_0.1.0.tar.gz", repos = NULL, type = "source")
```
<br>

### Dependencies  
Please ensure that you have the following packages installed. The glmnet package is required to fit ridge regressions. In order to run ocrRBBR with parallel computing, the packages doParallel, foreach, and doSNOW need to be installed.

```R
install.packages("glmnet")
install.packages("doParallel")  
install.packages("foreach")
install.packages("doSNOW")
```  

<br>

### Step 2. Prepare input files

```R
### Load data
library(RBBR)
library(readxl)

atacseq <- as.data.frame(read.csv(file = "ImmGenATAC18_AllOCRsInfo.csv", header= TRUE, check.names = FALSE))
rnaseq <- as.data.frame(read.csv(file = "mmc2.csv", header= TRUE, check.names = FALSE))
mmc1 <- as.data.frame(read_excel("mmc1.xlsx", sheet = 1, col_names = TRUE, col_types = "text"))
cell_type_lineage <- mmc1[ ,c(2,4,5)]

### Extract shared cell types between ATAC-seq and RNA-seq data
cells_types <- intersect(colnames(atacseq), colnames(rnaseq))
```

### Step 3. Extract ATAC-seq signal intensities and normalize per peak
```R
atacseq_data <- atacseq[   ,(colnames(atacseq) %in% cells_types)]
peak_names <- rownames(atacseq_data)

atacseq_data <- t(atacseq_data)
colnames(atacseq_data) <- peak_names

atacseq_data <- log(1+atacseq_data,10)
atacseq_data_scaled <- atacseq_data
for(j in 1:ncol(atacseq_data)){
  x <- ( atacseq_data[ ,j] - mean(atacseq_data[ ,j]) )/sd(atacseq_data[ ,j])
  atacseq_data_scaled[ ,j]<- 1/(1+exp(-x))
}
```

### Step 4. Extract ATAC-seq peaks within ±100 kb of the target gene TSS
```R
gene_id   <- "Rag2"

matched_indices <- grep(paste0("(?<!\\w)", gene_id, "(?!\\w)"), atacseq$genes.within.100Kb, perl = TRUE)
atacseq_gene <- atacseq[matched_indices, ]
atacseq_gene <- atacseq_gene[   ,(colnames(atacseq_gene) %in% cells_types)]
```

### Step 5. Remove ATAC-seq peaks with low signal intensities (based on p-values) or peaks not conserved across the mammalian genome. This step helps reduce potential false positive predictions by ocrRBBR and can be omitted if desired.
```R
peak_info <- atacseq[rownames(atacseq_gene), 1:8]

a <- median(peak_info$mm10.60way.phastCons_scores)
b <- median(peak_info$`_-log10_bestPvalue`)

peak_info <- peak_info[ ((peak_info$mm10.60way.phastCons_scores>=a)&(peak_info$`_-log10_bestPvalue`>=b)), ]
log_atacseq <- atacseq_data_scaled[ ,row.names(peak_info)]
```

### Step 6. Extract and rescale RNA-seq data for the target gene across blood cell lineages.
```R
rnaseq_gene <- rnaseq[(rnaseq$X %in% gene_id), ]
rnaseq_gene <- rnaseq_gene[ ,(colnames(rnaseq_gene) %in% cells_types)]
rownames(rnaseq_gene) <- gene_id

log_rnaseq <- log(1+rnaseq_gene,10) - min(log(1+rnaseq_gene,10))
log_rnaseq <- log_rnaseq/quantile(as.numeric(unlist(log_rnaseq)), probs = 0.975 , na.rm = TRUE)

data_scaled <- cbind( log_atacseq, t(log_rnaseq) )
data_scaled <- replace(data_scaled, data_scaled>=1, 0.9999)
data_scaled <- replace(data_scaled, data_scaled<=0, 0.0001)

head(data_scaled)
```

### Step 7. Train the model and output the predicted Boolean regulatory rules.
```R
rbbr           <- rbbr_train(data_scaled, max_feature = min(3,ncol(data_scaled)-1), mode = "1L", slope = 10, penalty = NA, weight_threshold = NA, num_cores = NA)
training process started with  8  computing cores
  |====================| 100%


head(rbbr$boolean_rules_sorted)

                                                                                                        Boolean_Rule                R2       BIC Input_Size Index             Features
1                              [OR(AND(278352,278381,278384),AND(~278352,278381,278384),AND(278352,~278381,278384))] 0.788633167796687 -306.8806          3   706 278352.278381.278384
2 [OR(AND(278362,278381,278398),AND(~278362,278381,278398),AND(278362,~278381,278398),AND(~278362,~278381,~278398))] 0.788220565514429 -306.7148          3  1435 278362.278381.278398
3 [OR(AND(278381,278390,278398),AND(278381,~278390,278398),AND(~278381,~278390,278398),AND(~278381,278390,~278398))] 0.788150615447657 -306.6867          3  2166 278381.278390.278398
4                              [OR(AND(278355,278381,278384),AND(~278355,278381,278384),AND(278355,~278381,278384))] 0.786059559598788 -305.8519          3  1277 278355.278381.278384
5                                                                                               [AND(278386,278398)] 0.746314992241634 -304.7223          2   275        278386.278398
6                                                                                                           [278384] 0.725155753428715 -304.6730          1    16               278384
  Active_Conjunctions                   Weights Layer1, Sub-Rule1
1                   3 0.46:0.72:0.34:-2.45:-2.17:-0.4:-0.64:-0.05
2                   4  0.49:0.75:0.55:-1.6:-2.02:-0.48:-0.67:0.05
3                   4 0.18:-1.59:0.93:-0.35:0.58:0.13:-1.33:-0.85
4                   3 0.54:0.75:0.34:-1.99:-1.88:-0.49:-0.6:-0.08
5                   1                       0.71:-0.8:-1.04:-0.28
6                   1                                  0.54:-0.54
```






























