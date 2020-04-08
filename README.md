# ExprX - an R package to streamline interspecies differential expression analysis
---
__ExprX__ is an R package to streamline interspecies differential expression
analysis. Using RNA-Seq data for human and mouse brain, this vignette
demonstrates how to detect differentially expressed genes among species
using ExprX package. All essential steps (eg. data import, ortholog matching,
data normalization, differential analysis and visualization) are described below.

Use the __install_git__ function from devtools package to install __ExprX__ from
GitHub:
---
library(devtools)

install_git("https://github.com/mingansun/ExprX")

To load the ExprX package:
---
library(ExprX)

Citation:
Ming-an Sun et al., ExprX - an R package to streamline interspe-cies differential gene expression analysis.
