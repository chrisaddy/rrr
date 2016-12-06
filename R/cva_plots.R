#' Rank Trace for Reduced-Rank CVA
#'
#' \code{cva_rank_trace}
#' 
#' @inheritParams cva
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_rank_trace(galaxy_x, galaxy_y)
#'
#' @export

cva_rank_trace <- function(x, y, type = "cov", k = 0){
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	rank_trace(x, y, gamma, type, k)
}

#' Rank Trace Plot for Reduced-Rank CVA
#'
#' \code{cva_rank_trace_plot}
#'
#' @inheritParams cva_rank_trace
#' @inheritParams rank_trace_plot
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_rank_trace_plot(galaxy_x, galaxy_y)
#'
#' @export

cva_rank_trace_plot <- function(x, y, type = "cov", k = 0, interactive = TRUE){
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	static_plot <- rank_trace_plot(x, y, gamma, type, k) + 
	ggtitle(paste("CVA Rank Trace Plot, k = ", k, sep = ""))
	if(interactive == TRUE){
		plotly::ggplotly(static_plot)
	} else {
		static_plot
	}
}

#' Residual Plots for Reduced-Rank CVA
#'
#' \code{cva_residual_plot}
#' 
#' @inheritParams cva_rank_trace_plot
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_residual_plot(galaxy_x, galaxy_y, rank = 3)
#'
#' @export

cva_residual_plot <- function(x, y, rank = "full", type = "cov", k = 0, interactive = FALSE){
	residuals <- cva_residuals(x, y, rank, type, k)
	static_plot <- GGally::ggpairs(residuals) + 
		labs(title = "CVA Residuals")
	static_plot[1,1] <- ggplot2::ggplot()
	if(interactive == TRUE){
		plotly::ggplotly(static_plot)
	} else {
		static_plot
	}
}

#' Pairwise Canonical Variates Plot
#'
#' @inheritParams cva
#' @inheritParams cva_rank_trace_plot
#' @param cva_pair integer. canonical variate pair to be plotted
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_pairwise_plot(galaxy_x, galaxy_y)
#'
#' @export 

cva_pairwise_plot <- function(x, y, cva_pair = 1, type = "cov", k = 0, interactive = TRUE){
	scores_object <- cva_scores(x, y, cva_pair, type, k)
	corr <- scores_object[["canonical_corr"]][cva_pair]
	x_axis <- scores_object[["xi"]][,cva_pair]
	y_axis <- scores_object[["omega"]][,cva_pair]
	df <- bind_cols(x_axis, y_axis)
	ggplot(df, aes_q(x = as.name(names(df)[1]), y = as.name(names(df)[2]))) + 
		geom_point() + 
		geom_smooth(method = "lm") + 
		labs(title = paste("CV", cva_pair, " Pairwise Plot, Canonical Correlation = ", round(corr, 4), sep = ""))
}


cva_allpairs_plot <- function(x, y, rank, type = "cov", k = 0){
	scores_object <- cva_scores(x, y, rank, type, k)
	all_pairs <- dplyr::bind_cols(scores_object[["xi"]], scores_object[["omega"]])
	GGally::ggpairs(all_pairs)
}