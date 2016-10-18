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

rrr <- function(x, y, gamma_matrix, rank = "full"){
		if(rank == "full"){
			reduce_rank <- min(dim(x)[2], dim(x)[2])
		} else {
			reduce_rank <- rank
		}
		x_organize <- organize(x)
		y_organize <- organize(y)
		weighted_matrix <- sqrt_matrix(gamma_matrix) %>%
								cov_matrix(y_organize, x_organize) %*%
								solve(cov_matrix(x_organize, x_organize)) %*%
								cov_matrix(x_organize, y_organize) %*%
								sqrt_matrix(gamma_mat)
		V_t <- eigen(weighted_matrix)[,1:rank] %>% as.matrix(ncol = rank)
		A_t <- solve(sqrt_matrix(gamma_mat)) %*% 
								V_t(x_organize, y_organize, gamma_mat, reduce_rank)
		B_t <- t(V_t(x_organize, organize(y), gamma_mat, reduce_rank)) %*% 
								sqrt_matrix(gamma_mat) %*% 
								cov_matrix(organize(y), organize(x)) %*% 
								solve(cov_matrix(organize(x), organize(x)))
		A_t
}