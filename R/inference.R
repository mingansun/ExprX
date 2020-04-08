################################################################################
# perform differential expression analysis
#------------------------------------------------------------------------------
# ortholog_expression_compare
# DE_rankprod
# DE_limma
# DE_ttest
# expression_statistics
################################################################################

#' Differential expression analysis of 1:1 orthologs between species
#'
#' This function performs differential expression analysis of 1:1 orthologs
#' between species. Expression data of 1:1 orthologs after normalization are
#' used for differential expression analysis.Statistics such as average
#' expression level, log2foldChange and p-values are calculated and returned
#' as a dataframe.
#'
#' @param expr_data
#' ExprX object which contains normalized expression data. It can be generated
#' with \code{\link{expression_matrix_normalize}}.
#'
#' @param method
#' Method for differential expression analysis. Can be RankProd, limma, t-test.
#'
#' @param p_adjust
#' Method for p-value adjustment. Can be "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none". P-value adjustment is
#' performed using the \code{p.adjust} function. Default: fdr.
#'
#' @return
#' Dataframe, with basic information (eg. GeneID, GeneName) of 1:1 orthologs and
#' statistics (eg. log2foldChange, adjusted p-value) from differential analysis.
#' The returned matrix can be filtered based on P-value and log2foldChange to
#' get differential genes.
#'
#' @examples
#' hs2mm.de <- ortholog_expression_compare(hs2mm.data, method = "RankProd", p_adjust = "fdr")
#'
#' @export
ortholog_expression_compare <- function(expr_data, method, p_adjust){

  # expr_data
  if(missing(expr_data)){
    stop("expr_data isn't given.")
  }
  expr_data <- verify_ExprX_dataset(expr_data)

  # method
  if(missing(method)){
    stop("method isn't specified.")
  }
  if(!is.vector(method)){
    stop("method should be a vector of length 1")
  }
  if(length(method) != 1){
    stop("method has length of ", length(method), ". Expect 1 element.")
  }

  methods.ok <- c("RankProd", "t-test", "limma")
  if(!method %in% methods.ok){
    stop(method, "is not supported. It should be among:\n",
         paste(methods.ok, collapse = "\n")
    )
  }

  # p_adjust: set as fdr if undefined
  if(missing(p_adjust)){
    message("p_adjust is not specified. Use default: fdr")
    p_adjust <- "fdr"
  }

  ##### examine and assign data

  # make sure OrthologExprNorm data exists in expr_data
  if(!"OrthologExprNorm" %in% names(expr_data)){
    stop("expr_data has no OrthologExprNorm. Generated it using expression_matrix_normalize first.")
  }
  expr_nrm <- expr_data$OrthologExprNorm

  # make sure expr_data has 2 species
  if(!"Species" %in% names(attributes(expr_nrm)) ) {
    stop("Species attribute is not available. Make sure expr_data is generated with expression_matrix_normalize.")
  }

  species_list <- attr(expr_nrm, "Species")
  species_uniq <- unique(species_list)
  verify_group_list(species_list)

  ##### perform differential analysis #####
  # multiple methods can be chosen for use

  deg_out <- NULL

  ###################### perform RankProduct test #############################
  # ------------------------------------------------------------------------- #
  # RankProduct test is species, because it returns two sets of p-values (ie.
  # p_value_up and p_value_down). Each set is helpful for determining up- and
  # down-regulated genes, respectively.
  #
  # For these reason, these two sets of p-values are combined for output:
  # (1) if a gene is up-regulated, keep the p_value_up.
  # (2) Otherwise, keep the p_value_down.
  #
  # The extraction of significant differential genes won't be affected by the
  # merging of p-values.
  #
  # Importantly, after merging the two sets of p-valuues, the output table of
  # RankProduct test can be of the same format of other tests. It will also
  # make downstream analysis much easier.
  # ------------------------------------------------------------------------- #
  #############################################################################
  if(method == "RankProd"){
    deg_out <- DE_rankprod(expr_nrm, species_list = species_list)
    # adjust p-value if required
    # RankProd returns two sets of p-values, for the testing of up- and down-
    # regulation, respectively. To be consistent with other methods, the two
    # sets of p-values are merged using a simple procedure:
    # 1) if log2fold>0, keep p-values for testing up-regulation
    # 2) otherwise, keep p-values for testing down-regulation
    # the calling of DEGs won't be affected by the merging of p-values.
    if(p_adjust != "none") {
      deg_out$pval <- ifelse(
        deg_out$log2fold>0,
        p.adjust(deg_out$pval_up,   method = p_adjust),
        p.adjust(deg_out$pval_down, method = p_adjust)
      )
    }
    else{
      deg_out$pval <- ifelse(
        deg_out$log2fold>0,
        deg_out$pval_up,
        deg_out$pval_down
      )
    }
  }

  ############## perform differential analysis using limma ####################
  # ------------------------------------------------------------------------- #
  # limma itself maybe not the optimal method for analyzing TPM/FPKM/RPKM
  # values, because the mean vs. variance relationship cannot be established
  # properly in this case. For RNA-Seq, the mean vs. variance relationship is
  # more evident for gene-level read counts, but not TPM/FPKM/RPKM values
  #
  # however, limma is still incorporated here for two reasons:
  # 1. provide a potentially helpful way to analyze microarray data for
  #    interspecies differential analysis
  # 2. based on read-counts for each gene, we may design a way to establish the
  #    mean vs. variance relationship by designing some approach similar to
  #    voom. Then, limma can be applied for interspecies analysis of RNA-Seq
  #    data. However, we decided not to do that now because the incorpotation
  #    of read count data will make the data processing methods (so far only
  #    need STAR or alternatives) way more complicated
  #
  # ------------------------------------------------------------------------- #
  #############################################################################
  if(method == "limma"){
    deg_out <- DE_limma(expr_nrm, species_list = species_list, log = TRUE)
    # adjust p-value if required
    if(p_adjust != "none") {
      deg_out$pval   <- p.adjust(deg_out$pval,   method = p_adjust)
    }
  }

  ############### perform t-test ##############################################
  # ------------------------------------------------------------------------- #
  # perform t-test through each ortholog genes between species, then adjust
  # the obtained p-values.
  #
  # Since t-test is performed on each gene - without borrowing information
  # from data of other genes, it is highly sensitive to the noise in the data.
  #
  # Comparison against RankProd results using real data between human and mouse
  # show that the p-values from these two methods can be highly correlated
  # (r ~ 0.6, p<2.2e-16), but the overlapping of DEGs (fdr<0.05, log2fold>1) is
  # moderate.
  #
  # Manual inspection of DEGs from both methods show that the obtained DEG list
  # from RankProd is more reasonable and fits better with previous knowledge.
  # ------------------------------------------------------------------------- #
  #############################################################################
  if(method == "t-test"){
    deg_out <- DE_ttest(expr_nrm, species_list = species_list)
    # adjust p-value if required
    if(p_adjust != "none") {
      deg_out$pval <- p.adjust(deg_out$pval, method = p_adjust)
    }
  }

  ##### reformat DE results and return as a dataframe
  deg.result <- data.frame(
    species1_geneid   = sub(":.*$", "", row.names(deg_out)),
    species2_geneid   = sub("^.*:", "", row.names(deg_out)),
    species1_genename = attr(expr_nrm, "GeneName")[[1]],
    species2_genename = attr(expr_nrm, "GeneName")[[2]],
    species1_value    = deg_out$value1,
    species2_value    = deg_out$value0,
    average_value     = deg_out$average,
    log2fold          = deg_out$log2fold,
    pvalue            = deg_out$pval
  )

  colnames(deg.result) <- c(
    paste0("GeneID_",     species_uniq[1]),
    paste0("GeneID_",     species_uniq[2]),
    paste0("GeneName_",   species_uniq[1]),
    paste0("GeneName_",   species_uniq[2]),
    paste0("Expression_", species_uniq[1]),
    paste0("Expression_", species_uniq[2]),
    "Expression_average",
    "log2foldChange",
    "P_value"
  )

  return(deg.result)
}

#' Differential expression analysis using RankProd approach
#'
#' @param expr_matrix
#' Matrix, with gene expression data. The rows are for genes, and columns for
#' samples. The samples should be from two conditions (ie species), as indicated
#' by species_list.
#' @param species_list
#' List of species, which should match the columns of expr_matrix. There are
#' should be two species, each with at least 2 replicates.
#'
#' @return
#' Data frame, with p-values and other statistics from differential analysis.
#'
#' @examples
#' hs2mm.RP_result <- DE_rankprod(hs2mm.expr.norm, hs2mm.species_list)
DE_rankprod <- function(expr_matrix, species_list){
  #
  if(missing(expr_matrix)){
    stop("expr_matrix is not given.")
  }
  verify_matrix(expr_matrix)

  #
  if(missing(species_list)){
    stop("species_list is not given.")
  }
  verify_group_list(species_list)

  if(length(species_list) != ncol(expr_matrix)){
    stop("Length of x doesn't equal ncol(expr_matrix): ", length(species_list), " != ", ncol(expr_matrix))
  }

  expr_grp <- rep(c(1,0), table(species_list))

  deg_out  <- RankProducts(log2(expr_matrix+0.5), expr_grp, logged = TRUE, plot=FALSE, fast = FALSE, rand = 229)
  deg_pval_up = deg_out$pval[,1]
  deg_pval_dn = deg_out$pval[,2]

  expr_stat <- expression_statistics(expr_matrix, expr_grp)
  expr_stat$pval_up   <- deg_pval_up
  expr_stat$pval_down <- deg_pval_dn

  return(expr_stat)
}

#' Differential expression analysis using limma approach
#'
#' @inheritParams DE_rankprod
#' @param log2_scale
#' TRUE/FALSE, indicate if the expression data should be log2-scaled before use.
#'
#' @return
#' Data frame, with p-values and other statistics from differential analysis.
#'
#' @examples
#' hs2mm.limma_result <- DE_limma(hs2mm.expr.norm, hs2mm.species_list)
DE_limma <- function(expr_matrix, species_list, log2_scale = TRUE){
  #
  if(missing(expr_matrix)){
    stop("expr_matrix is not given.")
  }
  verify_matrix(expr_matrix)

  #
  if(missing(species_list)){
    stop("species_list is not given.")
  }
  verify_group_list(species_list)

  if(length(species_list) != ncol(expr_matrix)){
    stop("Length of x doesn't equal ncol(expr_matrix): ", length(species_list), " != ", ncol(expr_matrix))
  }

  #
  if(missing(log2_scale)){
    log2_scale = TRUE
  }
  if(log2_scale != TRUE & log2_scale != FALSE){
    stop("log2_scale should be TRUE or FALSE.")
  }

  message("Perform differential analysis using limma")

  # design matrix
  expr_fac <- factor(species_list)
  design.mat <- model.matrix(~0 + ~expr_fac)
  colnames(design.mat) <- unique(species_list)

  # do log2-transformation or not based on the given parameter
  if(log == TRUE){
    eset <- log2(expr_matrix + 0.5)
  }
  else{
    eset <- expr_matrix
  }

  # linear fit
  fit <- lmFit(eset, design.mat)

  # DE analysis
  fit <- eBayes(fit)

  # p-value extraction
  expr_pval <- fit$F.p.value

  #
  expr_grp  <- rep(c(1,0), table(species_list))
  expr_stat <- expression_statistics(expr_matrix, expr_grp)
  expr_stat$pval   <- expr_pval

  return(expr_stat)
}

#' Differential expression analysis using t-test approach
#'
#' @inheritParams DE_rankprod
#'
#' @return
#' Data frame, with p-values and other statistics from differential analysis.
#'
#' @examples
#' hs2mm.tt_result <- DE_ttest(hs2mm.expr.norm, hs2mm.species_list)
DE_ttest <- function(expr_matrix, species_list){
  #
  if(missing(expr_matrix)){
    stop("expr_matrix is not given.")
  }
  verify_matrix(expr_matrix)

  #
  if(missing(species_list)){
    stop("species_list is not given.")
  }
  verify_group_list(species_list)

  if(length(species_list) != ncol(expr_matrix)){
    stop("Length of x doesn't equal ncol(expr_matrix): ", length(species_list), " != ", ncol(expr_matrix))
  }

  message("Perform t-test through all genes")

  expr_grp <- rep(c(1,0), table(species_list))
  expr_pval <- apply(expr_matrix, 1, function(v){x <- v[expr_grp==1]; y <- v[expr_grp!=1]; t.test(x,y)$p.value})
  # set as 1 if cannot calculated p-value (ie return NA)
  expr_pval[is.na(expr_pval)] <- 1

  expr_stat <- expression_statistics(expr_matrix, expr_grp)
  expr_stat$pval   <- expr_pval

  return(expr_stat)
}

#' Calculate statistics from expression matrix
#'
#' Based on the expression data matrix for differential expression analysis,
#' this function calculate basic statistics (eg. average expression level,
#' expression level for each species, log2foldChange etc.), which will be
#' returned together with differential analysis results (eg. P-values etc,
#' which are calculated else where).
#'
#' @inheritParams DE_rankprod
#'
#' @return
#' Data frame, with basic statistics (eg. average expression level,
#' expression level for each species, log2foldChange etc.) calculated from
#' expression data matrix.
#'
#' @examples
#' hs2mm.stat <- expression_statistics(hs2mm.expr.norm, hs2mm.species_list)
expression_statistics <- function(expr_matrix, group_list){
  #
  if(missing(expr_matrix)){
    stop("expr_matrix is not given.")
  }
  verify_matrix(expr_matrix)

  #
  if(missing(group_list)){
    stop("group_list is not given.")
  }
  verify_group_list(group_list)
  if(length(group_list) != ncol(expr_matrix)){
    stop(
      "Length of group_list doesn't equal ncol(expr_matrix): ",
      length(group_list), " != ", ncol(expr_matrix)
      )
  }

  # get expr for each group and their average
  expr_grp_1  <- expr_matrix[,group_list==1]
  expr_grp_0  <- expr_matrix[,group_list!=1]
  expr_mean_1 <- apply(expr_grp_1, 1, mean)
  expr_mean_0 <- apply(expr_grp_0, 1, mean)
  expr_mean   <- (expr_mean_1 + expr_mean_0) / 2

  # pseudocount of 0.5 is added to enable dividing when encounting 0
  log2fold <- log2( (expr_mean_1+0.5) / (expr_mean_0 + 0.5) )

  ##
  lst <- data.frame(value1 = expr_mean_1, value0 = expr_mean_0, average = expr_mean, log2fold = log2fold)
  return(lst)
}
