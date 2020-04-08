################################################################################
#                 generic functions
################################################################################

#' List all supported species
#'
#' Only species with genome available in ENSEMBL database are supported by
#' ExprX. This function will show the information for all these species.
#'
#' @return
#' Data frame, with the information for supported species
#'
#' @examples
#' list_species()
#'
#' @export
list_species <- function(){

  opts <- options(stringsAsFactors = FALSE)

  # read species list file from ExprX/extdata folder. If unexist, then try
  # ExprX/inst/extdata folder. Exit with error if fail again
  species_file <- paste0(path.package("ExprX"), "/extdata/species_list.csv")
  if(!file.exists(species_file)){
    species_file <- paste0(path.package("ExprX"), "/inst/extdata/species_list.csv")
  }
  if(!file.exists(species_file)){
    stop(getwd(), "\n", "Species list file doesn't exist: ", species_file)
  }

  df_species <- read.csv(species_file)
  return(df_species)
  on.exit(options(opts, add = TRUE))
}

#' Check if the Ensembl gene IDs are valid
#'
#' @param id
#' Vector. List of ENSEMBL gene IDs
#' @param species
#' Species name. Use list_species() to list all species.
#'
#' @return
#' Logical, return TRUE if all valid, otherwise alert and return FALSE
#'
#' @examples
#' valid_ensembl_geneid(hsa.id, "human")
valid_ensembl_geneid <- function(id, species){

  if (missing(id)){
    stop("id is not specified.")
  }
  if (missing(species)){
    stop("species is not specified.")
  }

  res <- grep("^ENS\\w+\\d+$", id)
  if (length(res)==0) {
    warning(paste0(id, ": the ID looks abnormal."))
    return(FALSE)
  }
  else{
    if(re != 1){
      warning(paste0(id, ": the ID looks abnormal."))
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  }
}

#' Extract IDs that are unique in both lists
#'
#' Based on two lists of IDs, this function check both of them and
#' return IDs that are unique in both lists. By applying to the
#' homolog annotation tabels of two species, it can can get
#' candidate 1-to-1 orthologue pairs which will be further
#' compared.
#'
#' @param x
#' Vector, ID list 1
#' @param y
#' Vector, ID list 2
#'
#' @examples
#' get_uniq_id_pair(human.ortholog.ann$ensembl_gene_id, mouse.ortholog.ann$hsapiens_homolog_ensembl_gene)
#' get_uniq_id_pair(human.ortholog.ann[,1], mouse.ortholog.ann[,2])
get_uniq_id_pair <- function(x,y){

  if(missing(x)){
    stop("x is not given.")
  }
  if(!is.vector(x)){
    stop("x should be a vector.")
  }
  if(missing(y)){
    stop("y is not given.")
  }
  if(!is.vector(y)){
    stop("y should be a vector.")
  }

  x.cnt <- table(x)
  y.cnt <- table(y)

  return(
    intersect(names(x.cnt)[x.cnt==1], names(y.cnt)[y.cnt==1])
  )
}

#' Verify if the object is a valid matrix
#'
#' @param x
#' Data Matrixm
#'
#' @return
#' No return values. Stop with error if failed.
#'
#' @examples
#' verify_matrix(mat)
verify_matrix <- function(x){
  if(missing(x)){
    stop("x not given.")
  }
  if(!is.matrix(x)){
    stop("x is not a matrix.")
  }
  if(nrow(x)<2){
    stop("x should have >=2 rows. But got ", nrow(x), ".")
  }
  if(ncol(x)<4){
    stop("x should have >=4 columns. But got ", ncol(x), ".")
  }
}

#' Check matrix to get the rows with NAs
#'
#' @inheritParams  verify_matrix
#'
#' @return
#' Vector, row numbers that contain NAs
#'
#' @examples
#' na_row_idx <- check_na_rows(x)
check_na_rows <- function(x){
  if(missing(x)){
    stop("x is not given.")
  }
  return(
    which(apply(x, 1, function(v){sum(is.na(v))} ) > 0)
  )
}

#' Verify if the group list is normal
#'
#' The input can be things like group classes as numbers, or species classes
#' as characters and so on. This function will check if there are multiple
#' groups, and each group should have >=2 replicates. If failed, it will exit
#' with errors.
#'
#' @param x
#' Vector, with class labels.
#'
#' @return
#' No return values. Stop with error if failed.
#'
#' @examples
#' verify_group_list(c(1,1,1,0,0))
#' verify_group_list(c("human", "human", "human", "mouse", "mouse"))
verify_group_list <- function(group_list){
  # check missing data
  if(missing(group_list)){
    stop("group_list is not given.")
  }
  if(!is.vector(group_list)){
    stop("group_list should be a vector.")
  }
  # require at least 2 groups
  if(length(unique(group_list)) != 2){
    stop(
      "unique(group_list) should have length of 2, but got ", length(unique(group_list)), ":\n",
      paste(unique(group_list), collapse = "\n")
      )
  }
  # require each group with >=2 replicates
  group_cnt <- table(group_list)
  if(min(table(group_list))<2){
    stop(
      "Each group in group_list should has >=2 elements, but:\n",
      paste(paste0(names(group_cnt), ":\t", group_cnt), collapse = "\n")
    )
  }
}
