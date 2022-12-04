#' Fit Reduced-Rank Regression Model
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x data frame or matrix of predictor variables
#' @param y data frame or matrix of response variables
#' @param type type of reduced-rank regression model to fit. \code{type = "identity"}, the default, uses \eqn{\mathbf{\Gamma} = \mathbf{I}} to fit a reduced-rank regression. \code{type = "pca"} fits a principal component analysis model as a special case of reduced-rank regression. \code{type = "cva"} fits a canonical variate analysis model as a special case of reduced-rank regression. \code{type = "lda"} fits a linear discriminant analysis model as a special case of reduced-rank regression.
#' @param rank rank of coefficient matrix.
#' @param k small constant added to diagonal of covariance matrices to make inversion easier.
#'
#' @return list containing estimates of coefficients and means, and eigenvalue-based diagnostics.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' rrr(tobacco_x, tobacco_y, rank = 1)
#'
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' rrr(digits_features, digits_features, type = "pca", rank = 3)
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' rrr(galaxy_x, galaxy_y, type = "cva", rank = 2)
#'
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' rrr(iris_x, iris_y, type = "lda")
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrr <- function(x, y, type = "identity", rank = "full", k = 0){
	if(is.matrix(type)){
		reduce_rank_regression(x, y, type, rank, k)
	} else {
		if(type == "lda"){
			full_rank <- dim(distinct(as_data_frame(y)))[1] - 1
		} else {
			full_rank <- dim(y)[2]
		}
		if(rank == "full"){
			reduced_rank <- full_rank
		} else {
			reduced_rank <- rank
		}
		ident <- diag(1, full_rank)
		switch(type,
			   identity = reduce_rank_regression(x, y, ident, reduced_rank, k),
			   pca = pca(x, reduced_rank, k),
			   cva = cva(x, y, reduced_rank, k),
			   lda = lda(x, y, reduced_rank, k))
	}
}

#' Rank Trace Plot
#'
#' \code{rank_trace} is a plot used to determine the effective dimensionality, i.e., \eqn{t = \mathrm{rank}\left(\mathbf{C}\right)},
#' of the reduced-rank regression equation.
#'
#' @inheritParams rrr
#' @param plot if FALSE, returns data frame of rank trace coordinates.
#' @param interactive if TRUE, creates an interactive plotly graphic.
#'
#' @return plot of rank trace coordinates if \code{plot = TRUE}, the default, or data frame of rank trace coordinates if \code{plot = FALSE}.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' rank_trace(tobacco_x, tobacco_y)
#' rank_trace(tobacco_x, tobacco_y, plot = FALSE)
#' rank_trace(tobacco_x, tobacco_y, type = "cva")
#'
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' rank_trace(digits_features, digits_features, type = "pca")
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' rank_trace(galaxy_x, galaxy_y, type = "cva")
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rank_trace <- function(x, y, type = "identity", k = 0, plot = TRUE, interactive = FALSE){
	ident <- diag(1, dim(y)[2])
	switch(type,
		   identity = rrr_rank_trace(x, y, ident, k, plot, interactive),
		   cva = cva_rank_trace(x, y, k, plot, interactive),
		   pca = pca_rank_trace(x, k, plot, interactive),
		   "type not recognized for this function")
}

#' Reduced-Rank Regression Residuals
#'
#' \code{residuals} calculates the regression residuals for reduced-rank regression and canonical variate analysis.
#'
#' @inheritParams rrr
#' @inheritParams rank_trace
#'
#' @return scatterplot matrix of residuals if \code{plot = TRUE}, the default, or a data frame of residuals if \code{plot = FALSE}.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' tobacco_rrr <- rrr(tobacco_x, tobacco_y, rank = 1)
#' residuals(tobacco_x, tobacco_y, rank = 1, plot = FALSE)
#' residuals(tobacco_x, tobacco_y, rank = 1)
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' residuals(galaxy_x, galaxy_y, type = "cva", rank = 2, k = 0.001)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

residuals <- function(x, y, type = "identity", rank = "full", k = 0, plot = TRUE){
	ident <- diag(1, dim(y)[2])
	switch(type,
		   identity = rrr_residual_plot(x, y, ident, rank, k, plot),
		   cva = cva_residual_plot(x, y, rank, k, plot),
		   "type not recognized for this function")
}

#' Pairwise Plots
#'
#' @inheritParams rrr
#' @param pair_x variable to be plotted on the \eqn{X}-axis
#' @param pair_y variable to be plotted on the \eqn{Y}-axis
#' @param interactive logical. If \code{interactive = FALSE}, the default, plots a static pairwise plot. If \code{interactive = FALSE} plots an interactive pairwise plot.
#' @param point_size size of points in scatter plot.
#'
#' @return ggplot2 object if \code{interactive = FALSE}; plotly object if \code{interactive = TRUE}.
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[,1:34]
#' digits_class <- pendigits[,35]
#' pairwise_plot(digits_features, digits_class, type = "pca", pair_x = 1, pair_y = 3)
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' pairwise_plot(galaxy_x, galaxy_y, type = "cva")
#'
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' pairwise_plot(iris_x, iris_y, type = "lda")
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

pairwise_plot <- function(x, y, type = "pca", pair_x = 1, pair_y = 2, rank = "full", k = 0, interactive = FALSE, point_size = 2.5){
	switch(type,
		   pca = pca_pairwise_plot(x, pair_x, pair_y, y, rank, k, interactive),
		   cva = cva_pairwise_plot(x, y, pair_x, k, interactive),
		   lda = lda_pairwise_plot(x, y, pair_x, pair_y, rank, k, interactive),
		   "type not recognized for this function")
}

#' Compute Latent Variable Scores
#'
#' @inheritParams rrr
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' scores(digits_features, digits_features, type = "pca", rank = 3)
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' scores(galaxy_x, galaxy_y, type = "cva", rank = 4)
#'
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' scores(iris_x, iris_y, type = "lda")
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

scores <- function(x, y, type = "pca", rank = "full", k = 0){
	switch(type,
		   pca = pca_scores(x, rank, k),
		   cva = cva_scores(x, y, rank, k),
		   lda = lda_scores(x, y, rank, k),
		   "type not recognized for this function")
}

#' All Pairs Plots
#'
#'
#'
#' @inheritParams pairwise_plot
#'
#' @return scatterplot matrix.
#'
#'

# plots all pairs of principal components
#
# @inheritParams scores
# @inheritParams pairwise_plot
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' digits_class <- pendigits[,35]
#' allpairs_plot(digits_features, digits_class, type = "pca", rank = 3)
#'
#' @export

allpairs_plot <- function(x, y, type = "pca", rank, k = 0){
	if(rank < 2){
		stop("Too few pairs for allpairs_plot. Choose rank >= 2")
	} else {
	switch(type,
		pca = pca_allpairs_plot(x, rank, k, y),
		lda = lda_pairwise_plot(x, y, rank, k))
	}
}


#' 3-D Reduced Rank Regression Plots
#'
#' Create three-dimensional, interactive plotly graphics for exploration and diagnostics.
#'
#' @inheritParams pairwise_plot
#' @param pair_z variable to be plotted on the \eqn{Z}-axis
#'
#' @return three-dimensional plot. If \code{type = "pca"} returns three principal components scores - defaulted to the first three - against each other.
#' If \code{type = "cva"} returns three-dimensional plot of residuals. If \code{type = "lda"} returns three-dimensional plot of three linear discriminant scores plotted against each other.
#'
#' @examples
#' \dontrun{
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' threewise_plot(digits_features, digits_class, type = "pca", k = 0.0001)
#'
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' threewise_plot(galaxy_x, galaxy_y, type = "cva")
#'
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' threewise_plot(iris_x, iris_y, type = "lda")
#' }
#'
#' @export

threewise_plot <- function(x, y, type = "pca", pair_x = 1, pair_y = 2, pair_z = 3, rank = "full", k = 0, point_size = 2.5){
	switch(type,
		   pca = pca_3D_plot(x, y, pair_x, pair_y, pair_z, rank, point_size),
		   cva = cva_residual_3D_plot(x, y, pair_x, pair_y, pair_z, rank, k, point_size),
		   lda = lda_3D_plot(x, y, pair_x, pair_y, pair_z, rank, k, point_size))
}
