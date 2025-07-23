# Introduction
This R software package used GSVA to calculate the enrichment abundance of 85 metabolic pathways through transcriptome data of new samples, and based on the metabolic-feature-based subtypes (MFSs) and enrichment abundance of 85 metabolic pathways of the built-in training set, the new samples were classified as MFS1, or MFS2 or MFS3.

# Usage
```
##Install package
devtools::install_github("Zhou-hs-lab/MFSs")

##Predict your data
library(MFSs)
pre.res <- MFS85(your.expr.data,gsva.methods = "gsva",ncores = 1,scale = FALSE)
```
