# ExprX - an R package to streamline interspecies differential expression analysis
---
__ExprX__ is an R package to streamline interspecies differential expression
analysis. Taking TPM or FPKM/RPKM files for samples from different species as
input, it provides functions to handle all the necessary steps, including data
loading, ortholog matching, normalization, differential analysis and visualization.

## How to install ExprX
__ExprX__ can be installed using the __install_git__ function from devtools
package. However, it depends on several other R packages, which should be 
installed first.

__To install ExprX and its dependencies:__
\# install several dependent packages from Bioconductor using BiocManager
install.packages("BiocManager")
BiocManager::install(c("biomaRt", "edgeR", "RankProd"))
\# install ExprX from GitHub using devtools
install.packages("devtools")
devtools::install_git("https://github.com/mingansun/ExprX")

__To load the ExprX package:__
library(ExprX)

## How to use ExprX
Please refer to the Vignette about how to use __ExprX__ to perform interspecies differential
expression analysis. The PDF version of Vignette can be downloaded from:

https://github.com/mingansun/ExprX/blob/master/inst/doc/ExprX-vignette.pdf

## Citation:
Ming-an Sun et al., ExprX - an R package to streamline interspecies differential gene expression analysis.
