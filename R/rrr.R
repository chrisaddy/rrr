#' Fit Reduced-Rank Regression Model
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x data frame or matrix of input variables
#' @param y data frame or matrix of response variables
#' @param gamma_matrix weight matrix
#' @param rank of the coefficient matrix to estimate. Default \code{rank = full} prodcues the standard multivariate regression technique.
#' @param type the format of the covariance matrix. \code{type = "cov"} runs the regression using the mean-centered covariance matrix. \code{type = "cor"} runs the regression using the mean-centered, standard-deviation-scaled correlation matrix.
#' @param k small number to add to the ridge.
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

rrr <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	if(rank == "full"){
		reduce_rank <- full_rank
	} else if(rank <= full_rank){
		reduce_rank <- rank
	} else {
		stop("rank out of bounds")
	}
	cov_x <- cov(x) + k * diag(1, dim(x)[2])
	cov_yx <- cov(y, x)
    cov_y <- cov(y) + k * diag(1, dim(y)[2])
    cov_xy <- t(cov_yx)
	sqrtm <- sqrt_matrix(gamma_matrix)
	weighted_matrix <- sqrtm %*%
		cov_yx %*%
		solve(cov_x) %*%
		cov_xy %*%
		sqrtm
	V_t <- eigen(weighted_matrix)$vectors[,1:reduce_rank] %>%
		as.matrix(ncol = reduce_rank)
	A_t <- solve(sqrtm) %*% V_t
	B_t <- t(V_t) %*%
		sqrtm %*%
		cov_yx %*%
		solve(cov_x)
	C_t <- A_t %*% B_t
	mu_y <- colMeans(y)
	mu_x <- colMeans(x)
	mu_t <- mu_y - C_t %*% mu_x
	list(mean = mu_t, A = A_t, B = B_t, C = C_t)
}

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
	y_new - rrr_predict(rrr_object, x_new)
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


