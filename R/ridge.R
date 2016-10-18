#' Ridge Regression
#'
#' \code{ridge} fits a generalized ridge regression model by minimizing sum of squared error subject to an elliptical restriction on model parameters.
#'
#' @param x data frame of input variables
#' @param y vector of response variables
#' @param k small constant to augment the diagonal of the covariance matrix \eqn{\Sigma_{XX}}
#' @radius radius of 

ridge <- function(x, y, k = 0){
			x_mat <- as.matrix(x)
			y_mat <- as.matrix(y)
			solve((t(x_mat) %*% x + k * diag(1, dim(x_mat)[2]))) %*% t(x_mat) %*% y_mat
}