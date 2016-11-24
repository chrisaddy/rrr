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
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	rank_trace_plot(x, y, gamma, type, k) + 
	ggtitle(paste("CVA Rank Trace Plot, k = ", k, sep = ""))
}

#cv_pairwise <- function(x, y, cv_pair, rank = "full", type = "cov")

#' Residual Plots for Reduced-Rank CVA
#'
#' \code{cv_residual_plot}
#' 
#' @inheritParams rrcva
#'
#' @export

cva_residual_plot <- function(x, y, rank = "full", type = "cov", k = 0){
	residuals <- cva_residuals(x, y, rank, type, k)
	GGally::ggpairs(residuals)
}

#' Pairwise Canonical Variates Plot
#'
#' @inheritParams rrcva
#' @param cv_pair integer. canonical variate pair to be plotted
#'
#' @export 

cv_pairwise_plot <- function(x, y, cv_pair = 1, type = "cov", k = 0){
	scores_object <- cv_scores(x, y, cv_pair, type, k)
	corr <- scores_object[["canonical_corr"]][cv_pair]
	x_axis <- scores_object$xi[,cv_pair]
	y_axis <- scores_object$omega[,cv_pair]
	df <- bind_cols(x_axis, y_axis)
	ggplot(df, aes_q(x = as.name(names(df)[1]), y = as.name(names(df)[2]))) + 
		geom_point() + 
		geom_smooth(method = "lm") + 
		labs(title = paste("CV", cv_pair, " Pairwise Plot, Canonical Correlation = ", round(corr, 4), sep = ""))
}


cv_allpairs_plot <- function(x, y, rank, type = "cov", k = 0){
	scores_object <- cv_scores(x, y, rank, type, k)
	all_pairs <- dplyr::bind_cols(scores_object[["xi"]], scores_object[["omega"]])
	GGally::ggpairs(all_pairs)
}