# ExprX - an R package to streamline interspecies differential expression analysis
---
__ExprX__ is an R package to streamline interspecies differential expression
analysis. Taking TPM or FPKM/RPKM files for samples from different species as
input, it provides functions to handle all the necessary steps, including data
loading, ortholog matching, normalization, differential analysis and visualization.

## Install ExprX
To install __ExprX__ from GitHub by using the __install_git__ function from __devtools__
package:

library(devtools)

install_git("https://github.com/mingansun/ExprX")

To load the __ExprX__ package:

library(ExprX)

## How to use ExprX
Please refer to the Vignette about how to use __ExprX__ to perform interspecies differential
expression analysis. The PDF version of Vignette can be downloaded from:

https://github.com/mingansun/ExprX/blob/master/inst/doc/ExprX-vignette.pdf

## Citation:
Ming-an Sun et al., ExprX - an R package to streamline interspecies differential gene expression analysis.
