## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # install several dependent packages from Bioconductor using BiocManager
#  install.packages("BiocManager")
#  BiocManager::install(c("biomaRt", "edgeR", "RankProd"))
#  
#  # install ExprX from GitHub using devtools
#  install.packages("devtools")
#  devtools::install_git("https://github.com/mingansun/ExprX")

## ----setup--------------------------------------------------------------------
library(ExprX)

## -----------------------------------------------------------------------------
# meta table file
hs2mm.meta_file <- paste0(path.package("ExprX"), "/extdata/brain_metatable.csv")

# make ExprX object from meta table
hs2mm.data  <- make_ExprX_dataset(
  hs2mm.meta_file,
  data_dir = paste0(path.package("ExprX"), "/extdata")
  )

## -----------------------------------------------------------------------------
sp.lst <- list_species()
head(sp.lst)

## ---- eval = FALSE------------------------------------------------------------
#  # Match 1:1 orthologs between human and mouse
#  hs2mm.orth <- ortholog_match("human", "mouse")
#  
#  # Save the ortholog results on hard disk for later use
#  saveRDS(hs2mm.orth, "hs2mm.orth.rds")

## -----------------------------------------------------------------------------
# Load ortholog result with readRDS
hs2mm.orth <- readRDS(paste0(path.package("ExprX"), "/data/hs2mm.orth.rds"))

# View the structure of hs2mm.orth
str(hs2mm.orth)

## -----------------------------------------------------------------------------
hs2mm.orth.genetype <- summarize_ortholog_gene(hs2mm.orth, group = "genetype")
head(hs2mm.orth.genetype)

## -----------------------------------------------------------------------------
hs2mm.orth.chrom <- summarize_ortholog_gene(hs2mm.orth, group = "chrom")
head(hs2mm.orth.chrom)

## -----------------------------------------------------------------------------
# Filter orthologs by excluding genes from chromosomes X, Y and MT, and only keep
# protein_coding genes
hs2mm.orth.flt <- ortholog_filter(
  hs2mm.orth,
  genetype_include = "protein_coding",
  chrom_exclude = c("X", "Y", "MT")
)

## -----------------------------------------------------------------------------
# Before filtering
dim(hs2mm.orth[[1]])

# After filtering
dim(hs2mm.orth.flt[[1]])


## -----------------------------------------------------------------------------
# Merge data
hs2mm.data <- ortholog_expression_merge(
  expr_data = hs2mm.data, orth_data = hs2mm.orth.flt
)

## -----------------------------------------------------------------------------
# Load required package
library(edgeR)

# Perform data normalization
hs2mm.data <- ortholog_expression_normalize(
  expr_data = hs2mm.data, method = "TMM"
)

## -----------------------------------------------------------------------------
# Load required package
library(RankProd)

# Perform differential analysis and then sort by p-value
hs2mm.deg <- ortholog_expression_compare(
  hs2mm.data, method = "RankProd", p_adjust = "fdr"
  )
hs2mm.deg <- hs2mm.deg[order(hs2mm.deg$P_value),]

# Check what the differential analysis result looks like
head(hs2mm.deg)

## -----------------------------------------------------------------------------
# Determine significant differential genes (human>mouse) based on p-values and log2foldChange
hs2mm.hsHigh <- subset(hs2mm.deg, subset = log2foldChange > 1  & P_value < 0.05)
head(hs2mm.hsHigh)

# Determine significant differential genes (human<mouse) based on p-values and log2foldChange
hs2mm.mmHigh <- subset(hs2mm.deg, subset = log2foldChange < -1 & P_value < 0.05)
head(hs2mm.mmHigh)

## ---- fig.width = 4, fig.height = 4-------------------------------------------
# Generate MA-plot
ortholog_expression_plot(
  hs2mm.deg, "MA",
  main = "Brain", xlim = c(0,120000), ylim = c(-18,18)
  )

## ---- fig.width = 4, fig.height = 4-------------------------------------------
# Generate Volcano-plot
ortholog_expression_plot(
  hs2mm.deg, "volcano",
  main = "Brain", xlim = c(-18,18), ylim = c(0, 3)
  )

