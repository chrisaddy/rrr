#' Reduced-Rank Regression
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x 
#' @param y
#' @param gamma_mat
#' @param rank rank of the coefficient matrix to estimate. \code{rank = full} is standard multivariate regression.
#'
#' @export rrr

rrr <- function(x, y, gamma_matrix, rank){
		x_organize <- organize(x)
		y_organize <- organize(y)
		sig_xx <- cov_matrix(x_organize, x_organize)
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
		A_t
}