# MVMR-covariated-GWAS
This repository contains the code to recover the result in "Bias mitigation in causal inference with covariate-adjusted GWAS using robust multivariable Mendelian randomization"

To recover the simulation result, readers should first generate SNP effects, then use the generated SNP effects for each scenario. The corresponding scenario of the code is indicated in the file name.
To conduct the real data analysis in paper, one need to extract the UK Biobank data, then perform GWAS and prune the SNP before MR estimation.
To conduct the real datat example in supplementary file, one need to download the data first, and run the "prepare_data.R" to reformat the data. Then prune the SNPs and do causal estimate.
