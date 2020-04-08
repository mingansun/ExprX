################################################################################
# update ExprX object by integrating ortholog matching results
#------------------------------------------------------------------------------
# ortholog_expression_merge
# expression_matrix_extract
################################################################################


#' Integrate the expression data of orthologues among species
#'
#' Take the original ExprX object generated with \code{\link{make_ExprX_dataset}}
#' and the orthologue pairing data generated using \code{\link{ortholog_match}}
#' as input, this function integrates together the expression data for all 1:1
#' orthologs among compared species. The integrated data will be appended to
#' the original ExprX object and returned as an updated object. To be noted, for
#' orthologs that don't have matched expression data, they will be excluded from
#' the returned data.
#'
#' @param expr_data
#' ExprX object. Contains the original expression data for the samples of
#' compared species. It can be generated with \code{\link{make_ExprX_dataset}}.
#' @param orth_data
#' Ortholog data. Contains the ortholog matching results among compared
#' species. It can be generated with \code{\link{ortholog_match}}.
#'
#' @return
#' Updated ExprX object. It has the same structure to expr_data, but with
#' the expression data matrix of 1:1 orthologs for all samples from compared
#' species appended.
#'
#' @examples
#' # update the ExprX object to include expression data matrix for matched
#' # 1:1 orthologs. The updated object can then be used for data normalization
#' # and differential analysis
#' hs2mm.data <- ortholog_expression_merge(hs2mm.data, hs2mm.orth)
#'
#' @export
ortholog_expression_merge <- function(expr_data, orth_data){
  # expr_data
  if(missing(expr_data)){
    stop("expr_data is not given.")
  }
  expr_data <- verify_ExprX_dataset(expr_data)
  if("Ortholog" %in% names(expr_data)){
    warning("Ortholog already exists in expr_data. It will be replaced.")
  }
  if("OrthologExpr" %in% names(expr_data)){
    warning("OrthologExpr already exists in expr_data. It will be replaced.")
  }

  # orth_data
  if(missing(orth_data)){
    stop("orth_data is not given.")
  }
  orth_data <- verify_ortholog_data(orth_data)

  # check if species match
  species_expr <- expr_data$Species
  species_orth <- attr(orth_data, "Species")

  if(length(species_expr) != length(species_orth)) {
    stop(
      "Species numbers differ between expr_data and orth_data:\n",
      "expr_data: n = ", length(species_expr), "\n",
      "orth_data: n = ", length(species_orth)
    )
  }

  if(length(setdiff(species_expr,species_orth)) != 0){
    stop("Species lists differ between expr_data and orth_data:\n",
         "expr_data: ", paste(species_expr, collapse = ", "), "\n",
         "orth_data: ", paste(species_orth, collapse = ", ")
    )
  }

  if(sum(!sort(species_expr) == sort(species_orth)) > 0){
    warning("Species lists differ in order between expr_data and orth_data:\n",
             "expr_data: ", paste(species_expr, collapse = ", "), "\n",
             "orth_data: ", paste(species_orth, collapse = ", ")
    )
  }

  ## only keep ortholog pairs (from orth_data) that also present in expr_data
  # get filter rule by checking through all species
  orth_flt_rule <- rep(TRUE, nrow(orth_data[[1]]))
  message("Number of ortholog pairs absent from expr_data")
  for(i in 1:length(orth_data)){
    this_flt_rule <- orth_data[[i]]$GeneID %in% row.names(expr_data$RawExpr[[i]])
    orth_flt_rule <- orth_flt_rule & this_flt_rule
    message(species_orth[i], ":\t", sum(!this_flt_rule))
  }

  # filter ortholog data based on filtering rule
  orth_data_flt <- orth_data
  for(i in 1:length(orth_data)){
    orth_data_flt[[i]] <- orth_data[[i]][orth_flt_rule,]
  }

  message("")
  message("Original  ortholog number: ", length(orth_flt_rule) )
  message("Discarded ortholog number: ", sum(!orth_flt_rule) )
  message("Resulted  ortholog number: ", sum( orth_flt_rule) )

  expr_data$Ortholog     <- orth_data_flt
  expr_data$OrthologExpr <- list()

  # sometimes not all genes of the ortholog pairs can be found in the expression matrix
  # need to pay special attention to make sure to return correct data w/o NA rows
  for(i in seq_along(species_expr)){
    species_this <- species_expr[i]

    # get index of species_this in orth_data
    orth_idx <- which(species_this == species_orth)
    orth_inf <- orth_data_flt[[orth_idx]]

    # get raw expr data for species_this
    expr_raw <- expr_data$RawExpr[[i]]
    # get raw expr data for ortholog genes
    expr_orth <- expr_raw[match(orth_inf$GeneID, row.names(expr_raw)),]

    # add raw expr of ortholog genes to $OrthologExpr
    expr_data$OrthologExpr[[i]] <- expr_orth
    names(expr_data$OrthologExpr)[i] <- species_this
    # keep gene names as attr, probably be used later
    attr(expr_data$OrthologExpr[[i]], "GeneName") <- orth_inf$GeneName[orth_inf$GeneID %in% row.names(expr_raw)]

  }

  return(expr_data)
}

#' Extract expression data matrix from ExprX object
#'
#' This function parse ExprX object to get the expression data matrix for
#' specific species or groups of genes. The input ExprX object can be
#' originally generated with \code{\link{make_ExprX_dataset}}, or updated
#' with \code{\link{ortholog_expression_merge}} or
#' \code{\link{ortholog_expression_normalize}}. The returned expression
#' data matrix can be normalized or not, depending on the input ExprX object.
#'
#' @param expr_data
#' ExprX object. It can be the ExprX object originally generated with
#' \code{\link{make_ExprX_dataset}}, or the updated object generated with
#' \code{\link{ortholog_expression_merge}} or
#' \code{\link{ortholog_expression_normalize}}.
#' @param species
#' Species name (eg. human, mouse etc). If not specified, data for all species
#' will be returned.
#' @param geneset
#' all or ortholog. Indicating the set of genes to be extracted. If set as all,
#' expression data for all genes will be returned. If set as ortholog, only
#' expression data for 1:1 orthologs will be returned. Expression data for
#' matched orthologs are only available after the ExprX objected has been
#' updated with \code{\link{ortholog_expression_merge}}.
#' @param normalized
#' TRUE or FALSE. Indicating if extract the normalized expression data. If set
#' as TRUE, normalized expression data will be returned. If set as FALSE, raw
#' expression values will be returned. To be noted, normalized expression data
#' are available for ortholog genes if only the input expr_data has been
#' processed by \code{\link{ortholog_expression_normalize}}.
#'
#' @return
#' Expression data matrix, with rows for genes and columns for samples.
#'
#' @examples
#' hs2mm.expr      <- expression_matrix_extract(hs2mm.orth, "all", "ortholog", FALSE)
#' hs2mm.expr.norm <- expression_matrix_extract(hs2mm.orth, "all", "ortholog", TRUE)
#' hs2mm.expr.hs   <- expression_matrix_extract(hs2mm.orth, "human", "ortholog", FALSE)
#' hs2mm.expr.mm   <- expression_matrix_extract(hs2mm.orth, "human", "ortholog", FALSE)
#' @export
expression_matrix_extract <- function(expr_data, species, geneset, normalized) {

  # expr_data
  if(missing(expr_data)){
    stop("expr_data isn't specified.")
  }
  expr_data <- verify_ExprX_dataset(expr_data)

  # species
  if(missing(species)){
    stop("species isn't specified.")
  }
  if(!is.vector(species)){
    stop("species should be a vector. But detected as ", typeof(species))
  }
  if(length(species) != 1){
    stop("species has length of ", length(species), ". Expect length of 1.")
  }
  if(!species %in% expr_data$Species){
    stop(
      species, " doesn't present in expr_data. It should be among: ",
      paste(expr_data$Species, collapse = "\n")
      )
  }

  # geneset
  if(missing(geneset)){
    stop("geneset is not specified. Should be all or ortholog.")
  }
  if(!is.vector(geneset)){
    stop("geneset should be a vector. But detected as ", typeof(geneset))
  }
  if(length(geneset) != 1){
    stop("geneset has length of ", length(geneset), ". Expect length of 1.")
  }
  if(!geneset %in% c("all", "ortholog")){
    stop("geneset is specified as ", geneset, ".It should be among: all or ortholog")
  }
  if(geneset == "ortholog" & !"OrthologExpr" %in% names(expr_data)){
    stop("expr_data lacks OrthologExpr data. Use ortholog_expression_merge to generate it first.")
  }

  # normalize
  if(missing(normalized)){
    stop("normalized is not specified.")
  }
  if(!normalized %in% c(TRUE, FALSE)){
    stop("normalized is specified as ", normalized, ".It should be TRUE or FALSE")
  }
  if(normalized == TRUE & !"OrthologExprNorm" %in% names(expr_data)){
    stop("expr_data lacks OrthologExprNorm data. Use ortholog_expression_normalize to generate it first.")
  }

  # check parameter conflict
  if(geneset == "all" & normalized == TRUE){
    stop("geneset = all and normalized = TRUE conflict.")
  }

  ### Get wanted data
  mat <- NULL
  species_idx <- which(expr_data$Species == species)

  ## expression of all genes
  if(geneset == "all"){
    mat <- expr_data$RawExpr[[species_idx]]
  }

  ## expression of ortholog genes
  if(geneset == "ortholog"){
    if(normalized == TRUE){
      mat <- expr_data$OrthologExprNorm[[species_idx]]
    }
    else{
      mat <- expr_data$OrthologExpr[[species_idx]]
    }
  }

  return(as.matrix(mat))
}
