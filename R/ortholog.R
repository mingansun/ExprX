################################################################################
#    match 1-to-1 orthologs among different species
# ------------------------------------------------------------------------------
# ortholog_match
# ortholog_filter
# verify_ortholog_data
# summarize_ortholog_gene
################################################################################

#' Match 1-to-1 ortholog between two species
#'
#' Ortholog matching is performed based on ortholog information retrieved from
#' ENSEMBL. The R package biomaRt is used to access ENSEMBL.To work properly,
#' make sure biomaRt package is installed and the internet connection is ok.
#' It may take sometime to download data from ENSEMBL - depending on the
#' network speed. Only species with genome availale in ENSEMBL database
#' are supported. Use \code{\link{list_species}} to get the information for
#' all supported species.
#'
#' @param species_1
#' Name of species 1 (eg human, mouse etc). Use \code{\link{list_species}} to
#' get all supported species.
#' @param species_2
#' Name of species 2 (eg human, mouse etc). Use \code{\link{list_species}} to
#' get all supported species.
#'
#' @return
#' List with information for 1-to-1 ortholog pairs. The list has two data
#' frames (each for one species) with columns GeneID, GeneName, Chrom,
#' GeneType. The genes in these two frames are paired orthologs of the same
#' order. Additonal attributes including Species and SpeciesAbbr are also
#' attached to the return list.
#'
#' @examples
#' hs2mm.orth <- ortholog_match("human", "mouse")
#' hs2gg.orth <- ortholog_match("human", "chicken")
#' hs2dr.orth <- ortholog_match("human", "zebrafish")
#'
#' @export
ortholog_match <- function(species_1, species_2){

  # check species names
  if(missing(species_1) & missing(species_2)){
    stop("species_1 and species_2 not specified.")
  }
  else{
    if(missing(species_1)) {
      stop("species_1 not specified.")
    }
    if(missing(species_2)) {
      stop("species_2 not specified.")
    }
  }

  # check if both species are supported (by comparing against list_species()),
  # and get their abbreviations which will be used by biomaRt
  species.ok <- list_species(updated = TRUE)
  species    <- c(tolower(species_1), tolower(species_2))
  species_abbr  <- NULL
  #species_full  <- NULL

  message("Get species name abbreviation to be used by biomaRt ...")
  for(sp in seq_along(species)) {
    if(! species[sp] %in% tolower(species.ok$Species)) {
      stop(species[sp], " is not supported. Use list_species() to list all supported species.")
    }
    species_abbr[sp] <- species.ok$Dataset[which(tolower(species.ok$Species) == species[sp])]
    #species_full[sp] <- species.ok$FullName[which(species.ok$Species == species[sp])]
    message(species[sp], " => ", species_abbr[sp])
  }

  ## retrieve annotation from Ensembl using biomaRt package
  if("biomaRt" %in% rownames(installed.packages())){
    library(biomaRt)
  }
  else{
    stop("biomaRt cannot be loaded. Install it first.")
  }

  mart=useMart("ensembl", verbose = TRUE)

  # biomaRt connection
  conn     <- list()
  hom_tag  <- paste0(species_abbr, "_homolog_ensembl_gene")
  attr_lst <- list()

  message("\nConnect to biomaRt for each species ...")
  for(sp in seq_along(species)){
    # connect
    message(species[sp], " => ", species_abbr[sp])
    conn[[sp]] <- useDataset(paste0(species_abbr[[sp]], "_gene_ensembl"), mart, verbose = FALSE)
    # only keep homolog attribute
    attr_lst[[sp]] <- grep("_homolog_ensembl_gene", listAttributes(conn[[sp]])$name, value = TRUE)
  }

  # pairwise check if the homolog table available for each species
  # if yes, continue for reciprocal comparison to get 1:1 pairs
  # otherwise stop
  message("\nDetermine the homolog table to be retrieved ...")
  for(sp1 in seq_along(species)){
    for(sp2 in seq_along(species)){
      # skip self comparison
      if(sp1 == sp2){
        next
      }
      if(hom_tag[sp2] %in% attr_lst[[sp1]]) {
        message(species[sp1], " => ", species[sp2], ":\t", hom_tag[sp2])
      }
      else{
        stop(paste0(species[sp1], " has no attribute: ", hom_tag[sp2]))
      }
    }
  }

  ## get gene and homolog annotation for each species
  ann.inf <- list()
  hom.inf <- list()
  for(i in seq_along(species)){
    # i & j are index for each species
    j = ifelse(i==1,2,1)
    message(paste0("\nRetrieving homolog annotation for ", species[i], " ..."))
    ann.inf[[i]] = getBM(
      attributes=c("ensembl_gene_id","external_gene_name", "chromosome_name", "gene_biotype"),
      mart=conn[[i]], filters="ensembl_gene_id", values=""
    )
    hom.inf[[i]] = getBM(
      attributes=c("ensembl_gene_id", hom_tag[[j]]),
      mart=conn[[i]], filters="ensembl_gene_id",
      values=unique(ann.inf[[i]]$ensembl_gene_id)
    )
  }

  ## match ortholog between two species

  # remove duplicated lines
  hom.inf.rmd <- list()
  for(i in 1:2){
    hom.inf.rmd[[i]] <- hom.inf[[i]][!duplicated(apply(hom.inf[[i]], 1, function(v){paste0(v[1], ":", v[2])})),]
  }

  # get IDs unique in both tables (ie. genes with 1-to-1 orthologs between sp1 & sp2)
  # check both table to make sure the ID only match to 1 ortholog gene
  # eg. for human gene IDs, it should:
  # occur only 1 time in human table hom.inf.rmd[[1]][,1]: thus match only 1 mouse gene
  # occur only 1 time in mouse table hom.inf.rmd[[2]][,2]: thus matched only by 1 mouse gene
  id.1t1.lst <- list()

  message("\nNumber of genes with only 1 match in other species:")
  for(i in seq_along(species)){
    j = ifelse(i == 1, 2, 1)
    id.1t1.lst[[i]] <- get_uniq_id_pair(hom.inf.rmd[[i]][,1], hom.inf.rmd[[j]][,2])
    message(species[i], ":\t", length(id.1t1.lst[[i]]))
  }

  # get homolog annotation table for 1-to-1 ortholog genes based on id.lt1.lst
  hom.inf.flt <- list()

  for(i in seq_along(species)){
    j = ifelse(i == 1, 2, 1)
    hom.inf.flt[[i]] <- hom.inf.rmd[[i]][hom.inf.rmd[[i]][,1] %in% id.1t1.lst[[i]] & hom.inf.rmd[[i]][,2] %in% id.1t1.lst[[j]],]
  }
  hom.pair.num <- NULL
  if(nrow(hom.inf.flt[[1]]) == nrow(hom.inf.flt[[2]])) {
    hom.pair.num <- nrow(hom.inf.flt[[1]])
    message("\nNumber of 1-to-1 ortholog pairs before filtering:\t", nrow(hom.inf.flt[[i]]))
  }
  else{
    stop(
      "Number of obtained 1-to-1 ortholog genes differ among species:\n",
      species[1], ": ", nrow(hom.inf.flt[[1]]), "\n",
      species[2], ": ", nrow(hom.inf.flt[[2]])
      )
  }

  # further filter to make sure the matched pairs are the same in both tables
  # these gene pairs are considered as correct
  hom.pair <- list()
  for(i in seq_along(species)) {
    j <- ifelse(i == 1, 2, 1)
    hom.pair[[i]] <- as.vector(apply(hom.inf.flt[[i]], 1, function(v){paste0(v[i], ':', v[j])}))
  }

  hom.pair.ok <- hom.inf.flt[[1]][hom.pair[[1]] %in% hom.pair[[2]],]
  colnames(hom.pair.ok) <- species
  message("Number of 1-to-1 ortholog pairs after  filtering:\t", nrow(hom.pair.ok), "\n")

  # add chrom and genetype annotation
  ann.pair <- list()
  for(i in 1:2){
    ann.pair[[i]] <- ann.inf[[i]][match(hom.pair.ok[,i], ann.inf[[i]]$ensembl_gene_id),]
  }

  # check if hom.pair.ok (ortholog pairs) and ann.pair (annotation for ortholog pairs) are of equal rows and of the same order
  hom.pair.inf <- list()
  for(i in seq_along(species)) {
    if( nrow(ann.pair[[i]]) == nrow(hom.pair.ok) & sum(ann.pair[[i]][,1] != hom.pair.ok[,i]) == 0 ) {
      hom.pair.inf[[i]] <- cbind(hom.pair.ok[,i], ann.pair[[i]][,-1])
      colnames(hom.pair.inf[[i]]) <- c("GeneID", "GeneName", "Chrom", "GeneType")
    }
    else{
      if( nrow(ann.pair[[i]]) != nrow(hom.pair.ok) ) {
        stop("Rows of ortholog pairs and gene annotation don't match for ", species[i], ": ", nrow(hom.pair.ok[[i]]), " != ", nrow(ann.pair[[i]]))
      }
      if( sum(ann.pair[[i]][,1] != hom.pair.ok[,i]) >0 ) {
        stop("Rows of ortholog pairs and gene annotation are of different order for species: ", species[i])
      }
    }
  }

  # return the data.frame with ortholog annotation and additional information for genes
  names(hom.pair.inf) <- species
  attr(hom.pair.inf, "Species")     <- species
  attr(hom.pair.inf, "SpeciesAbbr") <- species_abbr
  #attr(hom.pair.inf, "SpeciesFull") <- species_full

  return(hom.pair.inf)
}

#' Filter 1-to-1 orthologous pairs based on genetype or chromosomes
#'
#' Take the result generated using ortholog_match() as input, and perform
#' filtering based on gene type (eg. protein_coding), chromosome name, or
#' gene list. It will return the filtered data with the same data structure.
#'
#' @param x
#' Data with information for paired 1-to-1 orthologous pairs. It should be the
#' list generated using \code{\link{ortholog_match}}.
#' @param genetype_include
#' Gene types to be included. Use \code{\link{summarize_ortholog_gene}} to get
#' full list. genetype_include and genetype_exclude are exclusive.
#' @param genetype_exclude
#' Gene types to be excluded. Use \code{\link{summarize_ortholog_gene}} to get
#' full list.  genetype_include and genetype_exclude are exclusive.
#' @param chrom_include
#' Chromosomes to be included. Use \code{\link{summarize_ortholog_gene}} to
#' get full list. chrom_include and chrom_exclude are exclusive.
#' @param chrom_exclude
#' Chromosomes to be excluded. Use \code{\link{summarize_ortholog_gene}} to
#' get full list. chrom_include and chrom_exclude are exclusive.
#' @param geneid_include
#' Ensembl gene IDs to be included. geneid_include and geneid_exclude are
#' exclusive.
#' @param geneid_exclude
#' Ensembl gene IDs to be excluded. geneid_include and geneid_exclude are
#' exclusive.
#'
#' @return
#' List, filtered data with the same structre as the input data d.
#'
#' @examples
#' hs2mm.orth <- ortholog_match("human", "mouse")
#'
#' hs2mm.orth.flt1 <- ortholog_filter(hs2mm.orth, genetype_include = "protein_coding")
#' hs2mm.orth.flt2 <- ortholog_filter(hsa2mmu.df, chrom_exclude = c("X", "Y", "M"))
#' hs2mm.orth.flt3 <- ortholog_filter(
#'   hs2mm.orth,
#'   genetype_include = "protein_coding",
#'   chrom_exclude    = c("X", "Y", "M"),
#'   geneid_exclude   = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419")
#' )
#'
#' @export
ortholog_filter <- function(x, genetype_include, genetype_exclude, chrom_include, chrom_exclude, geneid_include, geneid_exclude){

  # check x
  if(missing(x)){
    stop("x not specified.")
  }
  x <- verify_ortholog_data(x)

  # filters
  if(
    missing(genetype_include) & missing(genetype_exclude) &
    missing(chrom_include)    & missing(chrom_exclude   ) &
    missing(genetype_include) & missing(geneid_exclude  )
    ) {
    stop("No filtering parameter is given.")
  }

  ## check if some filters conflict
  if(!missing(genetype_include) & !missing(genetype_exclude)){
    stop("genetype_include and genetype_exclude are exclusive.")
  }
  if(!missing(chrom_include) & !missing(chrom_exclude)){
    stop("chrom_include and chrom_exclude are exclusive.")
  }
  if(!missing(geneid_include) & !missing(geneid_exclude)){
    stop("geneid_include and geneid_exclude are exclusive.")
  }

  # construct filter rules
  flt.rule <- rep(TRUE, nrow(x[[1]]))
  for(i in 1:length(x)){
    message("Check data for ", attr(x, "Species")[i])
    d <- x[[i]]

    # gene type
    if(!missing(genetype_include)){
      flt <- d$GeneType %in% genetype_include
      flt.rule <- flt.rule & flt
      message("genetype_include matches:\t", sum(flt))
    }
    if(!missing(genetype_exclude)){
      flt <- !d$GeneType %in% genetype_include
      flt.rule <- flt.rule & flt
      message("genetype_exclude matches:\t", sum(!flt))
    }

    # chrom
    if(!missing(chrom_include)){
      flt <- d$Chrom %in% chrom_include
      flt.rule <- flt.rule & flt
      message("chrom_include    matches:\t", sum(flt))
    }
    if(!missing(chrom_exclude)){
      flt <- !d$Chrom %in% chrom_exclude
      flt.rule <- flt.rule & flt
      message("chrom_exclude    matches:\t", sum(!flt))
    }

    # gene ID
    if(!missing(geneid_include)){
      flt <- d$GeneID %in% gendid_include
      flt.rule <- flt.rule & flt
      message("geneid_include   matches:\t", sum(flt))
    }
    if(!missing(geneid_exclude)){
      flt <- !d$GeneID %in% geneid_exclude
      flt.rule <- flt.rule & flt
      message("geneid_exclude   matches:\t", sum(!flt))
    }
    message("")
  }

  message("Original gene number: ", nrow(x[[1]]) )
  message("Filtered gene number: ", sum(flt.rule))

  x.flt <- x
  for(i in 1:length(x)){
    x.flt[[i]] <- x[[i]][flt.rule,]
    # attr(x.flt[[i]], "GeneName") <- attr(x[[i]], "GeneName")[flt.rule]
  }

  return(x.flt)
}

#' Verify if the orthologue matching data is normal
#'
#' Examine given data to see if its format is normal. Further, if no "Species"
#' attribute is present, try to set as Species_1, Species_2 etc, and return
#' the updated data. The species name is required for some other functions.
#'
#' @param x
#' Data with information for paired 1-to-1 orthologous pairs. It should be the
#' list generated using \code{\link{ortholog_match}}.
#'
#' @return
#' List, unchanged/updated data with the same structre as the input data x.
#'
#' @examples
#' hs2mm.orth <- ortholog_match("human", "mouse")
#' hs2mm.orth <- verify_ortholog_data(hs2mm.orth)
#'
#' @export
verify_ortholog_data <- function (x) {

  if(missing(x)) {
    stop("x not specified.")
  }
  if(!is.list(x)) {
    stop("x should be an object generated using ortholog_match.")
  }
  if(length(x) == 0) {
    stop("x is empty.")
  }
  if(length(x) != 2) {
    warning("x has ", length(x), " items. Expect 2.")
  }
  for(i in 1:2) {
    if(nrow(x[[i]]) < 1){
      stop("Element ", i, " of x has 0 row.")
    }
    if(ncol(x[[i]]) != 4) {
      stop("Element ", i, " of x has ", ncol(x[[i]]), " column. Expect 4.")
    }
  }
  # set Species as Species_n if no "Species" attribute
  if(length(attr(x, "Species")) == 0) {
    attr(x,"Species") <- paste0("Species_", 1:length(x))
    warning("x has no Species attribute. Set to Species_1, Species_2 etc.")
  }
  if(length(attr(x, "Species")) != length(x)) {
    attr(x,"Species") <- paste0("Species_", 1:length(x))
    warning("attr(x,\"Species\") has different length to x. Change to Species_1, Species_2 etc.")
  }

  return(x)
}

#' Summarize the orthologues for each genetype or chromosome
#'
#' @param x
#' Data frame with ortholog annotations for two species. It is generated
#' using ortholog_match(species1, species2).
#' @param group
#' genetype or chrom, the number of genes for each group based on specified
#' tags will be listed.
#'
#' @return
#' Data frame with the count of genes for each genetype or chrom in each species.
#'
#' @examples
#' hs2mm.orth <- ortholog_match("human", "mouse")
#' summarize_ortholog_gene(hs2mm.orth, "genetype")
#'
#' @export
summarize_ortholog_gene <- function(x, group = c("genetype", "chrom")){

  # check x
  if(missing(x)){
    stop("x not specified.")
  }
  x <- verify_ortholog_data(x)

  # check group
  if(missing(group)){
    stop("group not specified. Should be: genetype or chrom.")
  }
  if(length(group)>1){
    stop("group should be: genetype or chrom.")
  }
  group <- tolower(group)
  if(!group %in% c("genetype", "chrom")){
    stop("group should be: genetype or crom.")
  }

  # retrieve and summary information
  species   <- attr(x, "Species")
  gene_grp  <- NULL
  gene_cnt  <- list()

  for(i in 1:length(x)) {
    if(group == "genetype") {
      gene_grp      <- unique(c(gene_grp, x[[i]]$GeneType))
      gene_cnt[[i]] <- table(x[[i]]$GeneType)
    }
    if(group == "chrom") {
      gene_grp      <- unique(c(gene_grp, x[[i]]$Chrom))
      gene_cnt[[i]] <- table(x[[i]]$Chrom)
    }
  }

  gene_grp <- sort(gene_grp)

  mat <- matrix(0, nrow = length(gene_grp), ncol = length(species))
  row.names(mat) <- gene_grp
  colnames(mat)  <- species

  for(i in seq_along(species)) {
    mat[match(names(gene_cnt[[i]]), gene_grp),i] = as.vector(gene_cnt[[i]])
  }

  return(mat)
}
