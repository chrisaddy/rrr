#' Predict  RRR 
#'
#' \code{predict_rrr} predicts a matrix of responses from the coefficients of a \code{rrr} object.
#'
#' @inheritParams rrr
#' @param rrr_object an object of type `list` that contains the means and coefficients of the reduced-rank regression.
#' @param x_new data frame or matrix of input variables used to predict the response matrix.
#'
#' @export

rrr_predict <- function(rrr_object, x_new){
	num_obs <- dim(x_new)[1]
	coeffs <- rrr_object$C
	means <- matrix(rep(rrr_object$mean, num_obs), ncol = num_obs)
	t(means + coeffs %*% organize(x_new)) %>% as_data_frame
}

#' RRR Error 
#'
#' \code{predict_rrr} predicts a matrix of responses from the coefficients of a \code{rrr} object.
#'
#' @inheritParams rrr_predict
#' @param y_new data frame or matrix of observed response variables.
#'
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
#' @export

rrr_residuals <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	object <- rrr(x, y, gamma_matrix, rank, type, k)
	rrr_error(object, x, y)
}