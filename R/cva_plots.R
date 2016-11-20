#' Rank Trace for Reduced-Rank CVA
#'
#' \code{cva_rank_trace}
#' 
#' @inheritParams rrcva
#'
#' @export

cva_rank_trace <- function(x, y, type = "cov", k = 0){
	gamma <- cov(y)
	rank_trace(x, y, gamma, type, k)
}

#' Rank Trace Plot for Reduced-Rank CVA
#'
#' \code{cva_rank_trace_plot}
#'
#' @inheritParams cva_rank_trace
#'
#' @export

cva_rank_trace_plot <- function(x, y, type = "cov", k = 0){
	gamma <- cov(y)
	rank_trace_plot(x, y, gamma, type, k) + 
	ggtitle(paste("CVA Rank Trace Plot, k = ", k, sep = ""))
}

#cv_pairwise <- function(x, y, cv_pair, rank = "full", type = "cov")

