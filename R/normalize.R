################################################################################
# update ExprX object by integrating normalized expression data
#------------------------------------------------------------------------------
# ortholog_expression_normalize
# normalize_matrix
################################################################################

## --------- expression matrix normalization -------------

#' Normalize the expression data of orthologs among species.
#'
#' Take ExprX object as input, it normalizes the expression data for 1:1
#' orthologs among samples of compared species. Different normalization
#' approaches are supported, including TMM, TMMwsp, RLE, upperquartile and
#' quantile. The implement of these approaches are based on edgeR and limma.
#' The normalized data matrix will be appended to the original ExprX object and
#' returned as the updated ExprX object.
#'
#' @param expr_data
#' ExprX object. Contains the expression data for the 1:1 orthologs for
#' all samples of compared species. It can be generated using
#' \code{\link{ortholog_expression_merge}}.
#' @param method
#' Normalization method. Supported methods include: TMM, TMMwsp, RLE,
#' upperquartile, quantile, none.
#'
#' @return
#' Updated ExprX object. It is of the same structure of expr_data, but with
#' normalized expression data appended.
#'
#' @examples
#' # update the ExprX dataset to included normalized expression data matrix
#' # for 1:1 orthologs. The updated object can then be used for differential
#' # expression analysis
#' hs2mm.data <- ortholog_expression_normalize(hs2mm.data, method = "TMM")
#'
#' @export
ortholog_expression_normalize <- function(expr_data, method) {

  # expr_data
  if(missing(expr_data)){
    stop("expr_data not specified.")
  }
  expr_data <- verify_ExprX_dataset(expr_data)

  # method
  if(missing(method)){
    stop("method not specified.")
  }
  if(!is.vector(method)){
    stop("method should be a vector of length 1.")
  }
  if(length(method) != 1){
    stop("method has length of ", length(method), ". Expect 1.")
  }

  methods.ok <- c("TMM","TMMwsp","RLE","upperquartile", "quantile", "none")
  if(!method %in% methods.ok){
    stop(method, "is not supported. It should be among:\n",
         paste(methods.ok, collapse = "\n")
    )
  }

  # combine ortholog expression matrix among species
  species      <- expr_data$Species
  mat_expr_raw <- NULL
  mat_rowname  <- NULL
  mat_genename <- list()
  mat_species  <- NULL
  for(i in seq_along(species)){
    this_expr <- expression_matrix_extract(
      expr_data = expr_data, species = species[i],
      geneset = "ortholog", norm = FALSE
      )
    mat_expr_raw      <- cbind(mat_expr_raw, this_expr)
    mat_rowname       <- paste0(mat_rowname, row.names(this_expr), sep = ":")
    mat_genename[[i]] <- attr(this_expr, "GeneName")
    mat_species       <- append(mat_species, rep(species[i], ncol(this_expr)))
  }
  row.names(mat_expr_raw)       <- sub(":$","",mat_rowname)
  attr(mat_expr_raw, "Species") <- mat_species
  expr_data$OrthologExprMatrix  <- mat_expr_raw

  # do normalization
  mat_expr_norm <- normalize_matrix(x = mat_expr_raw, method = method)
  attr(mat_expr_norm, "Species") <- mat_species
  expr_data$OrthologExprNorm     <- mat_expr_norm
  attr(expr_data$OrthologExprNorm, "GeneName") <- mat_genename

  return(expr_data)
}

#' Normalalize expression data matrix with specified method
#'
#' This function can normalize expression data matrix with specified approach.
#' Several approaches are implemented, including TMM, TMMwsp, RLE, upperquantile
#' and quantile. These approaches are implemented by invoking edgeR and limma.
#'
#' @param x
#' Expression data matrix, with rows for genes and columns for samples.
#' @param method
#' Method for normalization. It can be TMM, TMMwsp, RLE, upperquantile, quantile
#' or none.
#'
#' @return
#' Expression data matrix after normalization.
#'
#' @examples
#' mat_tpm.nrm  <- normalize_matrix(mat_tpm, "TMM")
#' mat_fpkm.nrm <- normalize_matrix(mat_fpkm, "quantile")
#'
#' @export
normalize_matrix <- function(x, method){

  # check expression data matrix
  if(missing(x)){
    stop("x not specified.")
  }
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(nrow(x)<2){
    stop("x has ", nrow(x), " rows.")
  }
  if(ncol(x)<2){
    stop("x has ", ncol(x), " columns.")
  }

  # check specifed method
  if(missing(method)){
    stop("method isn't specified.")
  }
  if(!is.vector(method)){
    stop("method should be a vector.")
  }
  if(length(method) != 1){
    stop("method has length of ", length(method), ". Expect 1.")
  }

  methods.ok <- c("TMM","TMMwsp","RLE","upperquartile","quantile", "none")
  if(!method %in% methods.ok){
    stop(
      method, "is not a supported method. Should be:\n",
      paste(methods.ok, collapse = "\n")
      )
  }

  x.nrm <- NULL

  ####################### normalize using edgeR ################################
  ##
  ## method used here:
  ## 1) first get CPM (thus sum of each column become 1M)
  ## 2) calculate norm.factor using calcNormFactor
  ## 3) normalize CPM (divide by norm.factor)
  ##
  ## alternative method:
  ## 1) convert the matrix to DGEList object
  ## 2) calculate both library.size and norm.factor using calcNormFactors
  ## 3) do normalizatin by dividing with both library.size and norm.factor
  ##
  ## Note/Warning:
  ##
  ## 1) the calculation of norm.factor is not affected by library size, thus
  ## the bias cannot be completely removed without also calculating CPM or
  ## considering library.size.
  ## 2) no matter what strategy (calculate CPM or use library.size), the
  ## differential expression analysis results should be the same.
  ##############################################################################
  if(method != 'quantile'){
    if("edgeR" %in% rownames(installed.packages())){
      library(edgeR)
    }
    else{
      stop("edgeR cannot be loaded. Install it first.")
    }
    # normalize
    x.cpm = edgeR::cpm(x)
    nrm.factor <- edgeR::calcNormFactors(object = x.cpm, method = method)
    x.nrm = t(t(x.cpm)/as.vector(nrm.factor))
    # indicate the method used in the attribute
    attr(x.nrm, "method") <- method
    attr(x.nrm, "norm.factor") <- nrm.factor
  }

  ######################### normalize using limma ##############################
  ##
  ## quantile is not the recommended approach for RNA-Seq data, but usually
  ## work nice for microarray data. By implement the quantile approach for
  ## normalization (and also the limma approach for differential expression
  ## analysis), ExprX will have the flexibility to be adopted for interspecies
  ## analysis of microarray data.
  ##
  ##############################################################################
  if(method == 'quantile'){
    if("limma" %in% rownames(installed.packages())){
      library(limma)
    }
    else{
      stop("limma cannot be loaded. Install it first.")
    }
    x.nrm <- normalizeQuantiles(x)
    row.names(x.nrm) <- row.names(x)
    colnames(x.nrm)  <- colnames(x)
  }

  return(x.nrm)
}
