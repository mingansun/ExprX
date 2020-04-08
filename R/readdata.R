################################################################################
#    read expression data files to create ExprX dataset
# ------------------------------------------------------------------------------
# make_ExprX_dataset
# make_meta_table
# save_meta_table
# read_meta_table
# read_geneid_from_file
# read_data_matrix_from_file
# verify_ExprX_dataset
################################################################################

#' Create ExprX dataset
#'
#' By parsing the meta table (as a data frame or CSV file) which contains
#' information about expression data files (usually contain TPM, FPKM or RPKM
#' values) for different species, this function will read these data files to
#' create an object which contains the expression levels of the replicates of
#' different species. The created ExprX object can also contain additional data
#' such as orthologue pairs, normalized expression etc, and will be used by
#' most of the subsequent analysis.
#'
#' @param x
#' Data frame or CSV file, which contains the information for each expression
#' data file for different species.The function can automatically determine if
#' the given x is a data frame or file name. If it's a data frame, then use
#' directly. Otherwise, it will read the information by reading the CSV file.
#' @param data_dir
#' The directory to find expression data files. If the absolute directory of
#' expression data files are specified, then this parameter won't be matter.
#' But if what is specified is relative directory (or no directory), then this
#' parameter will determine how to find the expression data files. Default: ./
#' (ie current working directory)
#'
#' @return
#' ExprX object, with information for expression values for each sample of
#' each species. The created object will be used by most of the analysis,
#' including orthologue matching, data normalization and differential
#' expression analysis. After each operation, the results can be appended to
#' the original ExprX object.
#'
#' @examples
#' # Make an ExprX object based on manually created data frame which contains
#' the meta information
#' hs2mm.data <- make_ExprX_dataset(hsa2mmu_meta_info)
#'
#' # make an ExprX object based on the CSV file with meta information
#' hs2mm.data <- make_ExprX_dataset("hsa2mmu_metatable.csv")
#'
#' @export
make_ExprX_dataset <- function(x, data_dir){

  opts <- options(stringsAsFactors = FALSE)

  # check x
  if(missing(x)){
    stop("x not given.")
  }

  # check and convert data_dir
  if(missing(data_dir)){
    data_dir = "."
  }
  if(!dir.exists(data_dir)){
    stop("data_dir unexist: ", data_dir)
  }
  data_dir <- sub("/$", "", data_dir)

  ## check if x is a data frame or file name
  meta.df <- NULL
  # if it's a data matrix, use directly
  if(is.data.frame(x)){
    message("x is a dataframe. Use directly.\n")
    meta.df <- x
  }
  # if it's a file name, read it to a dataframe for use
  else{
    if(is.character(x) & length(x) == 1){
      if(file.exists(x)){
        message("x is detected as a file name. Read to a data frame:\n")
        meta.df <- read_meta_table(x)
      }
      else{
        stop("x is detected as a file name, but doesn't exist.")
      }
    }
    else{
      stop("x should be a data frame or a file containing meta information.")
    }
  }

  ## display the meta data after formatting and sorting
  meta.df  <- meta.df[order(meta.df$Species, meta.df$RepIndex),]
  meta.fmt <- meta.df[,c(2:ncol(meta.df),1)]
  message(
    paste(colnames(meta.fmt), collapse = "\t"), "\n",
    paste(apply(meta.fmt, 1, function(x){paste(x, collapse = "\t")}), collapse = "\n")
  )

  ## examine if miss some required columns
  need.column <- c("File", "Species", "IdColumn", "ExprColumn", "ExprType", "RepIndex")
  meta.column <- names(meta.df)
  miss.column <- NULL
  # check through all required columns to find the missing ones
  for(i in need.column) {
    if(tolower(i) %in% tolower(meta.column)){
      names(meta.df)[which(tolower(i) == tolower(meta.column))] = i
    }
    else{
      miss.column <- append(miss.column, i)
    }
  }
  # alert and exit if lack some columns
  if(length(miss.column)>0){
    stop(
      "Lack required information:\n",
      paste(miss.column, collapse = "\n")
    )
  }

  ## examine species list
  if(length(unique(meta.df$Species)) != 2 | min(table(meta.df$Species)) < 2){
    stop("Require 2 species, each with >=2 replicates.")
  }

  ## read expression data files
  # check if all files exist
  # get their exact location by consider data_dir
  miss.file <- NULL
  for(i in seq_along(meta.df$File)){
    f <- meta.df$File[i]
    # adjust if relative directory is specified for the data files (ie. not started with /)
    f.abs <- grep("^/",f)
    if(length(f.abs) == 0){
      f <- paste0(data_dir, '/', f)
    }
    else{
      if(length(f.abs) == 1 & f.abs != 1){
        f <- paste0(data_dir, '/', f)
      }
    }
    if(!file.exists(f)){
      miss.file <- append(miss.file, f)
    }
    meta.df$File[i] <- f
  }
  if(length(miss.file) > 0){
    stop(
      "Some files unexist:\n",
      paste(miss.file, collapse = "\n")
    )
  }

  ## examine the type of expression data (tpm, fpkm, rpkm, etc)
  type.ok     <- c("tpm", "rpkm", "fpkm")
  type.unique <- unique(meta.df$ExprType)
  if(length(type.unique) == 1){
    if(!type.unique %in% type.ok){
      stop("Type of expression data should be among: tpm, rpkm, fpkm.")
    }
  }
  else{
    stop(
      "Type of expression data are inconsistent among the files:\n",
      paste(expr.type.unique, collapse = "\n")
      )
  }

  ## examine replicate and their indice
  for(sp in unique(meta.df$Species)){
    rep.id <- meta.df$RepIndex[meta.df$Species == sp]

    # replicate index for each species should be non-duplicated
    if(max(table(rep.id)) > 1){
      stop(
        "Replicate indexs are duplicated for: ", sp, "\n",
        paste(rep.id, "\n")
        )
    }

    # replicate indexes should start from 1 and be continuous
    if(min(rep.id) != 1 | max(rep.id) != length(rep.id)) {
      stop(
        "The replicated indice should start from 1 and continuous. But for ", sp, ":\n",
        paste(sort(rep.id), collapse = ", ")
        )
    }
  }

  ## make an object (list) to store all data, which will be used by most of the
  ## subsequent analyses. Hereby it is mentioned as an ExprX object.

  lst <- list()

  # description
  lst$Desc <- meta.df

  # species, replicates, expr type, files
  lst$Species   <- unique(meta.df$Species)
  lst$FullName  <- unique(meta.df$FullName)
  lst$AbbrName  <- unique(meta.df$AbbrName)
  lst$Replicate <- as.vector(table(meta.df$Species))
  lst$ExprType  <- unique(meta.df$Expr_type)
  lst$FileList  <- list(
    meta.df$File[meta.df$Species == lst$Species[1]],
    meta.df$File[meta.df$Species == lst$Species[2]]
  )
  names(lst$FileList) <- lst$Species

  # get gene list for each species
  message("\nRead gene IDs for each species ...")

  lst$FullGeneList <- list()
  lst$FullGeneNum  <- NULL

  for(i in seq_along(lst$Species)) {
    lst$FullGeneList[[i]] <- read_geneid_from_file(meta.df$File[meta.df$Species == lst$Species[i]], col_geneid = 1, header = TRUE)
    lst$FullGeneNum[i]    <- length(lst$FullGeneList[[i]])
    message(lst$Species[i], ": ", lst$FullGeneNum[i], " genes.")
  }
  names(lst$FullGeneList) <- lst$Species
  names(lst$FileList)     <- lst$Species

  ## get raw expression data for each species
  message("\nLoading expression data for each species ...")

  lst$RawExpr <- list()
  for(i in seq_along(lst$Species)){
    lst$RawExpr[[i]] <- read_data_matrix_from_file(
      lst$FileList[[i]], lst$Desc$IdColumn[1], lst$Desc$ExprColumn[1], header = TRUE, geneid_list = lst$FullGeneList[[i]]
      )
    colnames(lst$RawExpr[[i]]) <- paste0(lst$Species[i], "_", 1:lst$Replicate[i])
    message(lst$Species[i], ": ", lst$FullGeneNum[i], " genes from ", length(lst$FileList[[i]]), " files.")
  }

  names(lst$RawExpr) <- lst$Species
  return(lst)
  on.exit(options(opts))
}

#' Create meta table with information for expression data files of each species
#'
#' Based on given expression data file list, species list and information about
#' which columns have gene IDs and expression values, this function can create
#' a data frame with meta informaton for the expression data files for each
#' species.
#'
#' @param file_list
#' Vector, list of expression data files.
#' @param species_list
#' Vector, list of the species names (eg. human, mouse) which are in
#' corresponding to the expression data files. The list of all supported species
#' by ExprX can be obtained with the \code{\link{list_species}} function.
#' @param col_id
#' Integer, the column with gene ID in the expression data files. For most cases
#' (eg. RSEM results), the Gene IDs are in the 1st column.
#' @param col_expr
#' Integer, the column with expression values in the expression data files. For
#' example, the TPM values are in the 6th column of the RSEM result.
#' @param data_type
#' Type of gene expression data. Should be: tpm, fpkm or rpkm.
#'
#' @return
#' Data frame, with information for each expression data file for each species.
#'
#' @examples
#' hs2mm.meta <- make_meta_table(
#'   file_list = c(
#'   "tpm_dir/hsa_tpm_1.txt", "tpm_dir/hsa_tpm_2.txt",
#'   "tpm_dir/mmu_tpm_1.txt", "tpm_dir/mmu_tpm_2.txt"
#'   ),
#'   species_list = rep(c("human", "mouse"), each = 2),
#'   col_id = 1, col_expr = 6, data_type = "tpm"
#' )
#'
#' @export
make_meta_table <- function(file_list, species_list, col_id, col_expr, expr_type){

  # check list of file names
  if(missing(file_list)) {
    stop("file_list not specified.")
  }

  # check species list
  if(missing(species_list)) {
    stop("species_list not specified. Use list_species() to get all supported species.")
  }

  # make sure file list is a vector
  if(!is.vector(file_list)) {
    if(is.factor(file_list)) {
      file_list <- as.vector(as.character(file_list))
    }
    else{
      stop("file_list should be a vector.")
    }
  }

  # make sure species list is a vector
  if(!is.vector(species_list)){
    if(is.factor(species_list)){
      species_list <- as.vector(as.character(species_list))
    }
    else{
      stop("species_list should be a vector.")
    }
  }

  # make sure expression data file list and species list are of the same length and each
  # species has at least 2 replicates
  num.file    <- length(file_list)
  num.species <- length(species_list)
  cnt.species <- sort(table(species_list))

  if(num.file != num.species){
    stop("file_list and species_list are of different length: ", num.file, " != ", num.species)
  }

  if(min(table(species_list))<2){
    stop(
      "Require at least 2 replicates per species.\n",
      "But ", names(cnt.species)[1], "only has ", cnt.species[1]
      )
  }

  # check if some files unexist
  file.unexist <- NULL
  for(f in file_list){
    if(!file.exists(f)){
      file.unexist <- append(file.unexist, f)
    }
  }
  if(length(file.unexist) > 0){
    stop("These files unexist:\n", paste(file.unexist, collapse = "\n"))
  }

  # check if the specified data type is supported
  data_type <- tolower(data_type)
  data_type.ok <- c("tpm", "rpkm", "fpkm")
  if(!data_type %in% data_type.ok){
    stop(
      "The data type ", data_type, " invalid. ",
      "Should be among: ", paste(data_type.ok, collapse = ", ")
      )
  }

  # check if two species are given and both are supported
  species.uniq <- sort(unique(species_list))
  if(length(species.uniq) != 2) {
    stop(
      "Data should be for 2 species, but got ", species.num, " here:\n",
      paste(species.uniq, collapse = "\n")
      )
  }

  species.ok <- list_species()
  species.no <- NULL
  for(sp in species.uniq){
    if(! sp %in% species.ok$Species){
      species.no <- append(species.no, sp)
    }
  }
  if(length(species.no)>0){
    stop(
      "Some given species are unsupported:\n", paste(species.no, collapse = "\n"), "\n",
      "Use list_species() to get all supported species."
      )
  }

  # get scientific name & abbreviation
  species.full  <- species.ok$FullName[match(species, species.ok$Species)]
  species.abbr  <- species.ok$AbbrName[match(species, species.ok$Species)]

  # make a dataframe, then sort and add replicate index
  df <- data.frame(
    File           = file_list,
    Species        = species_list,
    FullName       = species.full,
    AbbrName       = species.abbr,
    IdColumn       = rep(col_id,    num.file),
    ExprColumn     = rep(col_expr,  num.file),
    ExprType       = rep(data_type, num.file)
  )
  df <- df[order(df$Species, df$File),]
  df$RepIndex <- c(1:table(df$Species)[1], 1:table(df$Species)[2])

  return(df)
}

#' Write data frame with meta information as a CSV file
#'
#' @param x
#' Data frame, with information about expression data files for each
#' species. It can be created with \code{\link{make_meta_table}} function.
#' @param filename
#' The name of the file to store the meta data.
#' @param overwrite
#' TRUE/FALSE, to overwrite an existed file or not.
#' Default: FALSE
#'
#' @examples
#' save_meta_table(hs2mm.meta, "hs2mm.metafile.csv")
#'
#' @export
save_meta_table <- function(x, filename, overwrite = FALSE){
  if(missing(x)){
    stop("x not specified.")
  }
  if(missing(filename)){
    stop("filename not specified.")
  }
  if(!is.data.frame(x)){
    stop("x should be a data frame.")
  }
  if(file.exists(filename)){
    if(overwrite == FALSE | missing(overwrite)){
      stop(filename, " already exists. Set overwrite=TRUE to overwrite it.")
    }
    else{
      warning(filename, " already exists. It is overwritten now.")
    }
  }
  write.csv(x, filename, row.names = FALSE)
}

#' Read meta information from file
#'
#' It reads the meta information from the file (CSV format) created manually or
#' with the \code{\link{save_meta_table}} function. Alternatively, the meta
#' information can be created directly with the \code{\link{make_meta_table}}
#' function.
#'
#' @param filename
#' File name, for the CSV file created manually or by using the
#' \code{\link{save_meta_table}} function.
#'
#' @return
#' Data frame, with information about expression data files for each species.
#'
#' @examples
#' hs2mm.meta <- read_meta_table("hs2mm.meta_table.csv")
#'
#' @export
read_meta_table <- function(filename){
  if(missing(filename)){
    stop("filename not specified.")
  }
  if(!is.vector(filename) | length(filename) != 1 | !is.character(filename)){
    stop("filename should be a file name.")
  }
  if(!file.exists(filename)){
    stop(filename, "doesn't exist.")
  }

  return(read.csv(filename))
}

#' Read expression data files to get gene IDs
#'
#' Read the given expression data files to get the gene IDs. If the gene IDs
#' are consistent among the files, return the them. Otherwise alert and exit.
#'
#' @inheritParams read_data_matrix_from_file
#'
#' @return
#' Vector with the full list of gene names.
#'
#' @examples
#' filelist <- c("file1.tpm.txt", "file2.tpm.txt")
#' geneid_list <- read_geneid_from_file(filelist, col_geneid = 1)
read_geneid_from_file <- function(file_list, col_geneid, header = TRUE){

  opts <- options(stringsAsFactors = FALSE)

  # check file_list
  if(missing(file_list)){
    stop("file_list not specified.")
  }

  if(is.vector(file_list)){
    if(length(file_list) == 0){
      stop("file_list is empty.")
    }
  }
  else{
    stop("file_list is not a vector.")
  }

  for(f in file_list){
    if(!file.exists(f)){
      message("File unexist: ", f)
    }
  }

  # check col_geneid
  if(!is.numeric(col_geneid)){
    stop("col_geneid should be an integer.")
  }

  # check header
  if(!header %in% c(TRUE, FALSE)){
    stop("header should be: TRUE or FALSE.")
  }

  # read and check gene IDs
  d1 <- read.table(file_list[1], header = header)
  id.lst <- d1[,col_geneid]

  if(length(file_list) >= 2){
    for(i in 2:length(file_list)){
      d.tmp  <- read.table(file_list[i], header = header)
      id.tmp <- d.tmp[,col_geneid]
      if(length(id.tmp) != length(id.lst) | sum(id.tmp %in% id.lst) != length(id.lst)){
        stop("Gene IDs inconsistent between the files:\n", file_list[i], "\n", file_list[i-1])
      }
    }
  }

  return(id.lst)
  on.exit(options(opts))
}

#' Read data files to create an expression matrix
#'
#' Read all given expression data files to get gene IDs and expression values,
#' then merge them to create a expression matrix.
#'
#' @param file_list
#' Vector, a list of expression data files
#' @param col_geneid
#' Integer, indicate which column has gene IDs
#' @param col_value
#' Integer, indicate which column has expression values
#' @param header
#' TRUE/FALSE, indicate if the files have header line. Default: TRUE.
#' @param geneid_list
#' Vector, a list of gene IDs. It is optional. If given, genes absent from the
#' list will be excluded. Otherwise, all genes read from the expression data
#' files will be kept.
#'
#' @return
#' Matrix, with expression values for each gene in each expression data file.
#' For the created matrix, each row is for a gene and each column for a sample.
#'
#' @examples
#' filelist <- c("file1.tpm.txt", "file2.tpm.txt")
#'
#' genelist <- c(
#'   "ENSG00000000003", "ENSG00000000005", "ENSG00000000419",
#'   "ENSG00000000457", "ENSG00000000460")
#'
#' expr.all  <- read_data_matrix_from_file(filelist, 1, 6)
#' expr.part <- read_data_matrix_from_file(filelist, 1, 6, genelist)
read_data_matrix_from_file <- function(file_list, col_geneid, col_value, header = TRUE, geneid_list) {

  opts <- options(stringsAsFactors = FALSE)

  # file_list: file names
  if(missing(file_list)){
    stop("file_list not specified.")
  }
  if(is.factor(file_list)){
    file_list <- is.vector(is.character(file_list))
  }
  if(length(file_list) == 0){
    stop("x is empty")
  }
  file.miss <- NULL
  for(x in file_list){
    if(!file.exists(x)){
      file.miss <- append(file.miss, x)
    }
  }
  if(length(file.miss) > 0){
    stop(
      "Some files cannot be found:\n",
      paste(file.miss, collapse = "\n")
    )
  }

  # col_geneid
  if(missing(col_geneid)){
    stop("col_geneid unspecified.")
  }
  if(!is.numeric(col_geneid)){
    stop("col_geneid should be an integer.")
  }

  # col_value
  if(missing(col_value)){
    stop("col_value unspecified.")
  }
  if(!is.numeric(col_value)){
    stop("col_value should be an integer.")
  }

  # header
  if(!header %in% c(TRUE, FALSE)){
    stop("header should be: TRUE or FALSE.")
  }

  # geneid
  id.lst <- NULL
  if(missing(geneid_list)) {
    d <- read.table(file[1], header = TRUE)
    id.lst <- d[,col_geneid]
  }
  else{
    if(is.factor(geneid_list)){
      geneid_list <- as.vector(as.character(geneid_list))
    }
    if(is.vector(geneid_list)){
      if(length(geneid_list) == 0){
        stop("geneid_list is empty.")
      }
      else{
        id.lst <- geneid_list
      }
    }
    else{
      stop("geneid_list should be a vector")
    }
  }

  ##
  n_row <- length(id.lst)
  n_col <- length(file_list)

  mat <- matrix(0, n_row, n_col)
  row.names(mat) <- id.lst
  colnames(mat)  <- file_list
  for(i in seq_along(file_list)) {
    d <- read.table(file_list[i], header = TRUE)
    id.tmp  <- d[,col_geneid]
    val.tmp <- d[,col_value]
    mat[,i] <- val.tmp[match(id.lst, id.tmp)]
  }

  return(mat)
  on.exit(options(opts, add = TRUE))
}

#' Examine if the ExprX dataset is valid
#'
#' Examine if the given ExprX object is valid, and print warnings if missing
#' some required content.
#'
#' @param x
#' ExprX object generated using \code{\link{make_ExprX_dataset}}.
#'
#' @return
#' ExprX object with the same structure of input x
#'
#' @examples
#' hs2mm.data <- verify_ExprX_dataset(hs2mm.data)
#' @export
verify_ExprX_dataset <- function (x) {
  if(missing(x)) {
    stop("x is not specified.")
  }
  if(!is.list(x)) {
    stop("x should be the ExprX object generated with make_ExprX_dataset.")
  }
  item.need <- c(
    "Desc", "Species", "FullName", "AbbrName", "Replicate",
    "FileList", "FullGeneList", "FullGeneNum", "RawExpr"
    )
  item.miss <- NULL
  for(i in item.need){
    if(! i %in% names(x)){
      item.miss <- append(item.miss, i)
      warning("x lacks the item: ", i)
    }
  }
  if(length(item.miss)>0){
    stop("x lacks some required items.")
  }
  return(x)
}
