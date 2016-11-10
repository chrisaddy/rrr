#' Fit Reduced-Rank Regression Model
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x data frame of input variables
#' @param y data frame of response variables
#' @param gamma_matrix weight matrix
#' @param rank of the coefficient matrix to estimate. Default \code{rank = full} is standard multivariate regression.
#' @param type the format of the covariance matrix. \code{type = "cov"} runs the regression using the mean-centered covariance matrix. \code{type = "cor"} runs the regression using the mean-centered, standard-deviation-scaled correlation matrix.
#' @param k small number to add to the ridge.
#'
#' @references Izenman, A.J. (2008) Modern Multivariate Statistical Techniques. Springer.
#' @export

rrr <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	if(rank == "full"){
		reduce_rank <- min(dim(x)[2], dim(y)[2])
	} else {
		reduce_rank <- rank
	}
	x_organize <- organize(x, type)
	y_organize <- organize(y, type)
	cov_x <- cov(x) #+ k * diag(1, dim(x_organize)[1])
	cov_yx <- cov(y, x)
    cov_y <- cov(y) #+ k * diag(1, dim(y_organize)[1])
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