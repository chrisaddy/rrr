#' Rank Trace Plot for Reduced-Rank CVA
#'
#' \code{cva_rank_trace} is a plot used to determine the effective dimensionality, i.e., \eqn{t = \mathrm{rank}\left(\mathbf{C}\right)},
#' of the reduced-rank regression equation, or the number of canonical variates to be used in the model.
#'
#' @inheritParams rank_trace
#'
#' @return ggplot object if plot is \code{TRUE}, a data frame of rank trace coordinates if plot = \code{FALSE},
#' or an interactive plotly object if interactive = \code{TRUE}.
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
#'
#' @seealso \code{\link{rank_trace}}
#'
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
#' \code{cva_residual_plot} is a scatter plot matrix used for diagnostics of the reduced-rank canonical variate analysis.
#' 
#' @inheritParams cva_rank_trace
#' @inheritParams cva
#'
#' @return ggplot object, scatterplot matrix.
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
#' @seealso \code{\link{rrr_residual_plot}}
#'
#' @export

cva_residual_plot <- function(x, y, rank = "full", type = "cov", k = 0){
	residuals <- cva_residual(x, y, rank, type, k)
	static_plot <- GGally::ggpairs(residuals) + 
		labs(title = "CVA Residuals")
	static_plot[1,1] <- ggplot2::ggplot()
}

#' 3D Plot of Residuals for Reduced-Rank CVA
#' 
#' \code{cva_residual_3D_plot} creates an interactive, html plotly plot that can be manipulated by the viewer.
#'
#' @return plotly object.
#'
#' @inheritParams cva_residual_plot
#' @param cva_x integer number of the canonical variate used for the x-axis.
#' @param cva_y integer number of the canonical variate used for the y-axis.
#' @param cva_z integer number of the canonical variate used for the z-axis.
#' @param point_size size of points in scatter.
#'
#' @export

cva_residual_3D_plot <- function(x, y, cva_x = 1, cva_y = 2, cva_z = 3, rank = "full", type = "cov", k = 0, point_size = 3){
	residuals <- cva_residual(x, y, rank, type, k)
	residuals[, 2:dim(residuals)[2]]
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
#' For a given canonical-variate pair, \code{cva_pairwise_plot} produces a scatter plot with the canonical variate scores of the predictor
#' variables along the \eqn{X}-axis against the canonical-variate scores of the response variables along the \code{Y}-axis.
#'
#' @inheritParams cva
#' @inheritParams cva_rank_trace
#' @param cv_pair integer. canonical variate pair to be plotted
#'
#' @return ggplot object if \code{interactive = FALSE}, the default, or an interactive plotly plot if \code{interactive = TRUE}.
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
#'
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

#cva_allpairs_plot <- function(x, y, rank, type = "cov", k = 0){
#	scores_object <- cva_scores(x, y, rank, type, k)
#	all_pairs <- dplyr::bind_cols(scores_object[["xi"]], scores_object[["omega"]])
#	GGally::ggpairs(all_pairs)
#}