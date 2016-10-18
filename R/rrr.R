V_t <- function(x, y, gamma_mat, rank){
				if(rank == "full"){
					reduce_rank <- min(dim(x)[2], dim(y)[2])
				} else {
					reduce_rank <- rank
				}
				w <- weighted_mat(x, y, gamma_mat)
				eigen(w)$vectors[,1:rank] %>%
				as.matrix(ncol = rank)
}

A_t <- function(x, y, gamma_mat, rank = "full"){
		if(rank == "full"){
			reduce_rank <- min(dim(x)[2], dim(y)[2])
		} else {
			reduce_rank <- rank
		}
		solve(sqrt_matrix(gamma_mat)) %*% 
			V_t(x, y, gamma_mat, reduce_rank)
}

B_t <- function(x, y, gamma_mat, rank = "full"){
		if(rank == "full"){
			reduce_rank <- min(dim(x)[2], dim(y)[2])
		} else {
			reduce_rank <- rank
		}
		t(V_t(organize(x), organize(y), gamma_mat, reduce_rank)) %*% sqrt_matrix(gamma_mat) %*% cov_mat(organize(y), organize(x)) %*% solve(cov_mat(organize(x), organize(x)))
}

C_t <- function(x, y, gamma_mat, rank) {
		if(rank == "full"){
			reduce_rank <- min(dim(x)[2], dim(y)[2])
		} else {
			reduce_rank <- rank
		}
		A_t(x, y, gamma_mat, reduce_rank) %*% 
			B_t(x, y, gamma_mat, reduce_rank)
}


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
		weighted_matrix <- sqrt_matrix(gamma_matrix)

}

weighted_mat <- function(x, y, gamma_mat){
			x_organize <- organize(x)
			y_organize <- organize(y)
			sqrt_matrix(gamma_mat) %*%
				cov_mat(y_organize, x_organize) %*%
				solve(cov_mat(x_organize, x_organize)) %*%
				cov_mat(x_organize, y_organize) %*%
				sqrt_matrix(gamma_mat)
}