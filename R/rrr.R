#' Reduced-Rank Regression
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x data frame of input variables
#' @param y data frame of response variables
#' @param gamma_mat weight matrix
#' @param rank of the coefficient matrix to estimate. \code{rank = full} is standard multivariate regression.
#' @param type the format of the covariance matrix. \code{type = "cov"} runs the regression using the mean-centered covariance matrix. \code{type = "cor"} runs the regression using the mean-centered, standard-deviation-scaled correlation matrix.
#'
#' @export rrr

rrr <- function(x, y, gamma_matrix, rank, type = "cov"){
	if(type == "cov"){
		x_organize <- organize(x)
		y_organize <- organize(y)
	} else if(type == "cor"){
		x_organize <- organize(x, scale = TRUE)
		y_organize <- organize(y, scale = TRUE)
	} else {
		stop("type input not recognized")
	}
	sig_xx <- cov_matrix(x_organize, x_organize)
	sig_yx <- cov_matrix(y_organize, x_organize)
	sqrtm <- sqrt_matrix(gamma_matrix)
	weighted_matrix <- sqrtm %*%
		sig_yx %*%
		solve(sig_xx) %*%
		t(sig_yx) %*%
		sqrtm
	V_t <- eigen(weighted_matrix)$vectors[,1:rank] %>% 
		as.matrix(ncol = rank)
	A_t <- solve(sqrtm) %*% V_t
	B_t <- t(V_t) %*% 
		sqrtm %*% 
		sig_yx %*% 
		solve(sig_xx)
	C_t <- A_t %*% B_t
	mu_y <- mu_vars(y_organize)
	mu_x <- mu_vars(x_organize)
	mu_t <- mu_y - C_t %*% mu_x
	list(mean = mu_t[,1], A = A_t, B = B_t, C = C_t)
}

rrr_predict <- function(x, y, x_new, gamma_matrix, rank){
	regress <- rrr(x, y, gamma_matrix, rank)
	regress$mean + regress$C %*% organize(x_new)				
}

rrr_residuals <- function(x, y, gamma_matrix, rank){
	resid <- rrr_predict(x, y, gamma_matrix, rank) - organize(y)
	mse <- sum(resid^2) / (dim(resid)[1] * dim(resid)[2])
	list(MSE = mse, Residuals = resid)
}
