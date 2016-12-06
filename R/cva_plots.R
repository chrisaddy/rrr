#' Rank Trace Plot for Reduced-Rank CVA
#'
#' \code{cva_rank_trace}
#'
#' @inheritParams rank_trace
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
#' cva_rank_trace(galaxy_x, galaxy_y, plot = FALSE)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

cva_rank_trace <- function(x, y, type = "cov", k = 0, plot = TRUE, interactive = FALSE){
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	if(plot == FALSE){
		rank_trace(x, y, gamma, type, k, plot = FALSE)
	}
	static_plot <- rank_trace(x, y, gamma, type, k) + 
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
#' @inheritParams cva_rank_trace
#' @inheritParams cva
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
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
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

#' 3D Plot of Residuals for Reduced-Rank CVA
#'
#'
#'
#' @inheritParams cva_residual_plot
#' @param cva_x integer number of the canonical variate used for the x-axis.
#' @param cva_y integer number of the canonical variate used for the y-axis.
#' @param cva_z integer number of the canonical variate used for the z-axis.
#' @param point_size size of points in scatter.
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

cva_residual_3D_plot <- function(x, y, cva_x = 1, cva_y = 2, cva_z = 3, rank = "full", type = "cov", k = 0, point_size = 3){
	residuals <- cva_residuals(x, y, rank, type, k) %>% dplyr::select(-index)
	resid_x <- residuals[, cva_x]
	resid_y <- residuals[, cva_y]
	resid_z <- residuals[, cva_z]
	resid_tbl <- dplyr::bind_cols(resid_x, resid_y, resid_z)
	names(resid_tbl) <- c("resid_x", "resid_y", "resid_z")
	ptsize <- point_size
	plotly::plot_ly(resid_tbl,
		x = ~resid_x,
		y = ~resid_y,
		z = ~resid_z,
		type = "scatter3d",
		mode = "marker",
		marker = list(size = ptsize),
		name = "CVA 3D Residual Scatter Plot")
	}

#' Pairwise Canonical Variates Plot
#'
#' @inheritParams cva
#' @inheritParams cva_rank_trace
#' @param cv_pair integer. canonical variate pair to be plotted
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
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export 

cva_pairwise_plot <- function(x, y, cv_pair = 1, type = "cov", k = 0, interactive = TRUE){
	scores_object <- cva_scores(x, y, cv_pair, type, k)
	corr <- scores_object[["canonical_corr"]][cv_pair]
	x_axis <- scores_object[["xi"]][,cv_pair]
	y_axis <- scores_object[["omega"]][,cv_pair]
	df <- bind_cols(x_axis, y_axis)
	ggplot(df, aes_q(x = as.name(names(df)[1]), y = as.name(names(df)[2]))) + 
		geom_point() + 
		geom_smooth(method = "lm") + 
		labs(title = paste("CV", cv_pair, " Pairwise Plot, Canonical Correlation = ", round(corr, 4), sep = ""))
}

cva_allpairs_plot <- function(x, y, rank, type = "cov", k = 0){
	scores_object <- cva_scores(x, y, rank, type, k)
	all_pairs <- dplyr::bind_cols(scores_object[["xi"]], scores_object[["omega"]])
	GGally::ggpairs(all_pairs)
}