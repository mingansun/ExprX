################################################################################
# Visualization of differential analysis results
# ------------------------------------------------------------------------------
# ortholog_expression_plot
# plot_MA
# plot_volcano
################################################################################

#' Visualize differential expression analysis result
#'
#' Visualize differential analysis results of the 1:1 orthologous genes
#' between two species. The differential genes are defined based on the
#' given cutoffs of p-value and log2fold. The differential expression can be
#' visualized as MA-plot or Volcano-plot, with differential genes
#' highlighted by color.
#'
#' @param d
#' Dataframe for differential analysis results. It can be generated with the
#' function \code{\link{ortholog_expression_compare}}.
#' @param type
#' The type of visualization. MA or volcano. Default: MA
#' @param cutoff_pvalue
#' The cutoff of p-value for differential genes. Set as 1 to keep all genes.
#' Default: 0.05.
#' @param cutoff_log2fold
#' The cutoff of log2fold for determine differential genes. Set as 0 to
#' keep all genes. Default: 1.
#' @param main
#' Main title, default: NULL
#' @param xlab
#' xlabel, default: log2fold(species1/species2)
#' @param ylab
#' ylabel, default: "Average expression" for MA-plot, and "-log10(p-value)"
#' for Volcano-plot
#' @param xlim
#' Range of x-axis, default: NULL
#' @param ylim
#' Range of y-axis, default: NULL
#'
#' @examples
#'
#' # Visualize differential analysis result as MA-plot, with genes with
#' # pvalue<0.05 highlighted
#' ortholog_expression_plot(hs2mm.de, type = "MA", cutoff_pvalue = 0.05)
#'
#' # Visualize differential analysis result as Valcano-plot, with genes with
#' # pvalue<0.05 and abs(log2fold)>1 highlighted
#' ortholog_expression_plot(
#'   hs2mm.de, type = "volcano",
#'   cutoff_pvalue = 0.05, cutoff_log2fold = 1
#' )
#'
#' @export
ortholog_expression_plot <- function( d, type, cutoff_pvalue = 0.05, cutoff_log2fold = 1, main, xlab, ylab, xlim, ylim ){

  ##### check parameters
  # d
  if(missing(d)){
    stop("d not given.")
  }
  if(!is.data.frame(d)){
    stop("d should be generated with ortholog_expression_compare.")
  }
  if(ncol(d) != 9){
    stop("d should have 9 columns, but got ", ncol(d))
  }
  if(nrow(d) < 2){
    stop("d should have at least 2 rows, but got ", nrow(d))
  }
  if(!"log2foldChange" %in% colnames(d)){
    stop("d doesn't have the column for log2foldChange")
  }
  if(!"P_value" %in% colnames(d)){
    stop("d doesn't have the column for P_value")
  }

  # type
  if(missing(type)){
    stop("type not defined")
  }
  type.ok <- c("MA", "volcano")
  if(!type %in% type.ok){
    stop("The plot type should be among: ", paste0(type.ok, collapse = ", "))
  }

  # pvalue cutoff
  if(missing(cutoff_pvalue)){
    cutoff_pvalue <- 0.05
  }

  # log2fold cutoff
  if(missing(cutoff_log2fold)){
    cutoff_log2fold <- 1
  }

  # main
  if(missing(main)){
    main <- paste0(sub("^.*_", "", colnames(d)[1]), "/", sub("^.*_", "", colnames(d)[2]))
  }
  if(missing(xlim)){
    xlim <- NULL
  }
  if(missing(ylim)){
    ylim <- NULL
  }

  ##### plot figure
  if(type == "MA"){
    if(missing(xlab)){
      xlab <- "Average expression"
    }
    if(missing(ylab)){
      ylab <- paste0("log2fold(", sub("^.*_", "", colnames(d)[1]), "/", sub("^.*_", "", colnames(d)[2]), ")")
    }
    plot_MA(d, cutoff_pvalue, cutoff_log2fold, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
  }

  if(type == "volcano"){
    if(missing(xlab)){
      xlab <- paste0("log2fold(", sub("^.*_", "", colnames(d)[1]), "/", sub("^.*_", "", colnames(d)[2]), ")")
    }
    if(missing(ylab)){
      ylab <- "-log10(p-value)"
    }
    plot_volcano(d, cutoff_pvalue, cutoff_log2fold, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
  }
}

#' Visualize differential analysis result as MA-plot
#'
#' @inheritParams ortholog_expression_plot
#' @param xlab
#' xlabel, default: Average expression
#' @param ylab
#' ylabel, default: log2foldChange
#'
#' @examples
#' plot_MA(deg, cutoff_pvalue = 0.05, cutoff_log2fold = 1)
plot_MA <- function(d, cutoff_pvalue, cutoff_log2fold, main, xlab, ylab, xlim, ylim){
  ## check parameters
  if(missing(d)){
    stop("d not given.")
  }
  if(missing(cutoff_pvalue)){
    warning("cutoff_pvalue not defined. Set as 0.05 by default.")
    cutoff_pvalue <- 0.05
  }
  if(missing(cutoff_log2fold)){
    warning("cutoff_log2fold not defined. Set as 1 by default.")
    cutoff_log2fold <- 1
  }
  if(missing(main)){
    main <- ""
  }
  if(missing(xlab)){
    xlab <- "Average expression"
  }
  if(missing(ylab)){
    ylab <- "log2foldChange"
  }
  if(missing(xlim)){
    xlim <- NULL
  }
  if(missing(ylim)){
    ylim <- NULL
  }

  ## MA-plot
  sig.status <- d$P_value <= cutoff_pvalue & abs(d$log2foldChange)>=cutoff_log2fold

  plot(
    d[!sig.status,7], d[!sig.status,8], pch = 19,
    main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
    cex = 0.1, col = "grey"
  )
  abline(h = 0, lwd = 2, col = "pink")
  points(
    d[sig.status,7], d[sig.status,8], pch = 19,
    cex = 0.2, col = "red"
  )
}

#' Visualize differential analysis result as volcano plot
#'
#' @inheritParams ortholog_expression_plot
#' @param xlab
#' xlabel, default: log2foldChange
#' @param ylab
#' ylabel, default: -log10(p-value)
#'
#' @examples
#' plot_volcano(deg, cutoff_pvalue = 0.05, cutoff_log2fold = 1)
plot_volcano <- function(d, cutoff_pvalue, cutoff_log2fold, main, xlab, ylab, xlim, ylim){
  ## check parameters
  if(missing(d)){
    stop("d not given.")
  }
  if(missing(cutoff_pvalue)){
    warning("cutoff_pvalue not defined. Set as 0.05 by default.")
    cutoff_pvalue <- 0.05
  }
  if(missing(cutoff_log2fold)){
    warning("cutoff_log2fold not defined. Set as 1 by default.")
    cutoff_log2fold <- 1
  }
  if(missing(main)){
    main = ""
  }
  if(missing(xlab)){
    xlab = "log2foldChange"
  }
  if(missing(ylab)){
    ylab = "-log10(p-value)"
  }
  if(missing(xlim)){
    xlim <- NULL
  }
  if(missing(ylim)){
    ylim <- NULL
  }

  ## volcano plot
  sig.status <- d$P_value <= cutoff_pvalue & abs(d$log2foldChange)>=cutoff_log2fold

  plot(
    d[!sig.status,8], -log10(d[!sig.status,9]), pch = 19,
    main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
    cex = 0.1, col = "grey"
  )
  abline(v = 0, lwd = 2, col = "pink")
  points(
    d[sig.status,8], -log10(d[sig.status,9]), pch = 19,
    cex = 0.2, col = "red"
  )
}
