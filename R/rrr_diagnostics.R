#' Predict Multivariate Responses via Reduced-Rank Regression
#'
#' \code{rrr_predict} predicts a matrix of responses from the coefficients of a \code{rrr} object.
#'
#' @inheritParams rrr
#' @param rrr_object an object of type `list` that contains the means and coefficients of the reduced-rank regression.
#' @param x_new data frame or matrix of input variables used to predict the response matrix.
#'
#' @examples
#' data(tobacco)
#' set.seed(123)
#' train_index <- sample(dim(tobacco)[1], 18)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' tobacco_train_x <- tobacco_x[train_index, ]
#' tobacco_train_y <- tobacco_y[train_index, ]
#' tobacco_test_x <- tobacco_x[-train_index, ]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' tobacco_rrr <- rrr(tobacco_train_x, tobacco_train_y, gamma, rank = 1)
#' rrr_predict(tobacco_rrr, tobacco_test_x)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

rrr_predict <- function(rrr_object, x_new){
	num_obs <- dim(x_new)[1]
	coeffs <- rrr_object[["C"]]
	means <- matrix(rep(rrr_object[["mean"]], num_obs), ncol = num_obs)
	df <- t(means + coeffs %*% organize(x_new)) %>% 
		as_data_frame()
}

#' Reduced-Rank Regression Error
#'
#' \code{rrr_error} calculates the error from predicting response variables from a test set using coefficients calculated from training data.
#'
#' @inheritParams rrr_predict
#' @param y_new data frame or matrix of observed response variables.
#'
#' @examples
#' data(tobacco)
#' set.seed(123)
#' train_index <- sample(dim(tobacco)[1], 18)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' tobacco_train_x <- tobacco_x[train_index, ]
#' tobacco_train_y <- tobacco_y[train_index, ]
#' tobacco_test_x <- tobacco_x[-train_index, ]
#' tobacco_test_y <- tobacco_y[-train_index, ]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' tobacco_rrr <- rrr(tobacco_train_x, tobacco_train_y, gamma, rank = 1)
#' rrr_error(tobacco_rrr, tobacco_test_x, tobacco_test_y)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

rrr_error <- function(rrr_object, x_new, y_new){
	as_data_frame(y_new - rrr_predict(rrr_object, x_new))
}

#' RRR Residuals
#' 
#' \code{rrr_residuals}
#'
#' @inheritParams rrr
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' tobacco_rrr <- rrr(tobacco_x, tobacco_y, gamma, rank = 1)
#' rrr_residuals(tobacco_x, tobacco_y, gamma, rank = 1)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

rrr_residuals <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	object <- rrr(x, y, gamma_matrix, rank, type, k)
	rrr_error(object, x, y)
}