# \code{cva} fits a reduced-rank canonical variate/correlation model. This is a special case of reduced-rank regression
# with the weight matrix set to the covariance matrix of \eqn{Y}, i.e., \eqn{\mathbf{\Gamma} = \mathbf{\Sigma}_{YY}}.
# Canonical variate analysis creates a set of new predictor variables that are linear combinations of the orignal predictors,
# and a set of new resonse variables that are linear combinations of the original responses, such that each pair of new predictor and new
# response maximizes correlation between the two and all pairs of canonical variates are independent of each other.
#
# @return list containing: matrix of means; matrix of canonical variate coefficients of predictor variables; matrix of canonical variate coefficients of response variables; vector of canonical correlations. 

cva <- function(x, y, rank = "full", k = 0) {
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	rrr_object <- reduce_rank_regression(x, y, gamma, rank, k)
	H <- ginv(rrr_object[["A"]])
	colnames(H) <- names(y)
    list(mean = rrr_object[["mean"]], G = rrr_object[["B"]], H = H, canonical_corr = rrr_object[["eigen_values"]])
}

# Canonical Variate Scores
#
# \code{cva_scores} creates linear combinations of predictor and response variables from the coefficients of reduced-rank canonical variate analysis.
# 
# @inheritParams cva
#
# @return list containing: data frame of canonical variate scores of predictor variables; data frame of canonical variate scores of response variables; vector of canonical correlations.
#
# @examples
# library(dplyr)
# data(COMBO17)
# galaxy <- as_data_frame(COMBO17)
# galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
# galaxy <- na.omit(galaxy)
# galaxy_x <- select(galaxy, -Rmag:-chi2red)
# galaxy_y <- select(galaxy, Rmag:chi2red)
# cva_scores(galaxy_x, galaxy_y, rank = 2)

cva_scores <- function(x, y, rank = "full", k = 0){
	cva_object <- cva(x, y, rank, k)
	correlation <- cva_object[["canonical_corr"]]
	xi <- as_data_frame(t(cva_object[["G"]] %*% organize(x)))
	names(xi) <- paste("xi", 1:dim(xi)[2], sep = "")
	omega <- as_data_frame(t(cva_object[["H"]] %*% organize(y)))
	names(omega) <- paste("omega", 1:dim(omega)[2], sep = "")
	list(xi = xi, omega = omega, canonical_corr = correlation)
}

# Error of Reduced-Rank Canonical Variate Analysis
#
# \code{cva_error} calculates the errors of a canonical-variate regression model built on a training set and applied to a test set.
# @return data frame of error.
#
# @examples
# set.seed(12345)
# library(dplyr)
# data(COMBO17)
# galaxy <- as_data_frame(COMBO17)
# galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
# galaxy <- na.omit(galaxy)
# galaxy_train <- 
# sample_size <- floor(0.66 * nrow(galaxy))
# train_ind <- sample(seq_len(nrow(galaxy)), size = sample_size)
# train <- galaxy[train_ind, ]
# test <- galaxy[-train_ind, ]
# train_x <- select(train, -Rmag:-chi2red)
# train_y <- select(train, Rmag:chi2red)
# test_x <- select(test, -Rmag:-chi2red)
# test_y <- select(test, Rmag:chi2red)
# cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)

cva_error <- function(x, y, x_new, y_new, rank = "full", k = 0){
	cva_object <- cva(x, y, rank, k)
	index <- data_frame(index = 1:dim(y_new)[1])
	error <- as_data_frame(t(cva_object[["H"]] %*% organize(y_new) - cva_object[["G"]] %*% organize(x_new)))
	names(error) <- paste("CV", 1:dim(error)[2], sep = "") 
	dplyr::bind_cols(index, error)
}

# Residuals of Reduced-Rank Canonical Variate Analysis
#
# \code{cva_residual} returns the multivariate residuals of the reduced-rank CVA regression.
#
# @return data frame of residuals.

cva_residual <- function(x, y, rank = "full", k = 0){
	cva_error(x, y, x, y, rank, k)	
}

# Residual Plots for Reduced-Rank CVA
#
# \code{cva_residual_plot} is a scatter plot matrix used for diagnostics of the reduced-rank canonical variate analysis.
# 
# @inheritParams cva_rank_trace
# @inheritParams cva
#
# @return ggplot object, scatterplot matrix.
#
# @examples
# library(dplyr)
# data(COMBO17)
# galaxy <- as_data_frame(COMBO17)
# galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
# galaxy <- na.omit(galaxy)
# galaxy_x <- select(galaxy, -Rmag:-chi2red)
# galaxy_y <- select(galaxy, Rmag:chi2red)
# cva_residual_plot(galaxy_x, galaxy_y, rank = 3)


cva_residual_plot <- function(x, y, rank = "full", k = 0, plot = TRUE){
	residuals <- cva_residual(x, y, rank, k)
	static_plot <- GGally::ggpairs(residuals) + 
		labs(title = "CVA Residuals")
	static_plot[1,1] <- ggplot2::ggplot()
	if(plot == FALSE){
		residuals
	} else {
		static_plot
	}
}

# \code{cva_rank_trace} is a plot used to determine the effective dimensionality, i.e., \eqn{t = \mathrm{rank}\left(\mathbf{C}\right)},
# of the reduced-rank regression equation, or the number of canonical variates to be used in the model.
#
# @return ggplot object if plot is \code{TRUE}, a data frame of rank trace coordinates if plot = \code{FALSE},
# or an interactive plotly object if interactive = \code{TRUE}.

cva_rank_trace <- function(x, y, k, plot, interactive){
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	if(plot == TRUE & interactive == TRUE){
		rt <- rrr_rank_trace(x, y, gamma, plot = TRUE, interactive = FALSE) +
			ggtitle(paste("CVA Rank Trace Plot, k = ", k, sep = ""))
		plotly::ggplotly(rt)
	} else if(plot == TRUE & interactive == FALSE){
		rrr_rank_trace(x, y, gamma, plot = TRUE, interactive = FALSE) +
			ggtitle(paste("CVA Rank Trace Plot, k = ", k, sep = ""))
		} else {
			rrr_rank_trace(x, y, gamma, plot = FALSE, interactive = FALSE)
	}
}


 # 3D Plot of Residuals for Reduced-Rank CVA
 
 # \code{cva_residual_3D_plot} creates an interactive, html plotly plot that can be manipulated by the viewer.

 # @return plotly object.

 # @param point_size size of points in scatter.


cva_residual_3D_plot <- function(x, y, cva_x = 1, cva_y = 2, cva_z = 3, rank = "full", k = 0, point_size = 3){
	residuals <- cva_residual(x, y, rank, k)
	resids <- residuals[, 2:dim(residuals)[2]]
	resid_x <- resids[, cva_x]
	resid_y <- resids[, cva_y]
	resid_z <- resids[, cva_z]
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


 # For a given canonical-variate pair, \code{cva_pairwise_plot} produces a scatter plot with the canonical variate scores of the predictor
 # variables along the \eqn{X}-axis against the canonical-variate scores of the response variables along the \code{Y}-axis.

 # @inheritParams cva
 # @inheritParams cva_rank_trace
 # @param cv_pair integer. canonical variate pair to be plotted

 # @return ggplot object if \code{interactive = FALSE}, the default, or an interactive plotly plot if \code{interactive = TRUE}.

 # @examples
 # library(dplyr)
 # data(COMBO17)
 # data(COMBO17)
 # galaxy <- as_data_frame(COMBO17)
 # galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
 # galaxy <- na.omit(galaxy)
 # galaxy_x <- select(galaxy, -Rmag:-chi2red)
 # galaxy_y <- select(galaxy, Rmag:chi2red)
 # cva_pairwise_plot(galaxy_x, galaxy_y)


cva_pairwise_plot <- function(x, y, cv_pair = 1, k = 0, interactive = FALSE){
	scores_object <- cva_scores(x, y, cv_pair, k)
	corr <- scores_object[["canonical_corr"]][cv_pair]
	x_axis <- scores_object[["xi"]][,cv_pair]
	y_axis <- scores_object[["omega"]][,cv_pair]
	df <- bind_cols(x_axis, y_axis)
	static_plot <- ggplot(df, aes_q(x = as.name(names(df)[1]), y = as.name(names(df)[2]))) + 
		geom_point() + 
		geom_smooth(method = "lm") + 
		labs(title = paste("CV", cv_pair, " Pairwise Plot, Canonical Correlation = ", round(corr, 4), sep = ""))
	if(interactive == TRUE){
		plotly::ggplotly(static_plot)
	} else {
		static_plot
	}
}