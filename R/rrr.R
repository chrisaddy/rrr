#' Fit Reduced-Rank Regression Model
#'
#' \code{rrr} fits a reduced-rank regression model.
#'
#' @param x data frame or matrix of predictor variables
#' @param y data frame or matrix of response variables
#' @param gamma_matrix weight matrix
#' @param rank of the coefficient matrix to estimate. Default \code{rank = full} prodcues the standard multivariate regression technique.
#' @param type the format of the covariance matrix. \code{type = "cov"} runs the regression using the mean-centered covariance matrix. \code{type = "cor"} runs the regression using the mean-centered, standard-deviation-scaled correlation matrix.
#' @param k small constant to add to the ridge of covariance matrices to avoid singularity.
#'
#' @return list containing: matrix of means; matrix \eqn{\mathbf{A}} with \code{rank} number of columns; matrix \eqn{\mathbf{B}} with \code{rank} number of rows; matrix \eqn{\mathbf{C}} of coefficients; vector of eigenvalues of weight matrix.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' rrr(tobacco_x, tobacco_y, gamma, rank = 1)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
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
	eigens <- eigen(weighted_matrix)
	eigen_values <- eigens[["values"]]
	V_t <- eigens[["vectors"]][,1:reduce_rank] %>%
		as.matrix(ncol = reduce_rank)
	A_t <- solve(sqrtm) %*% V_t
	rownames(A_t) <- names(y)
	B_t <- t(V_t) %*%
		sqrtm %*%
		cov_yx %*%
		solve(cov_x)
	C_t <- A_t %*% B_t
	mu_y <- colMeans(y)
	mu_x <- colMeans(x)
	mu_t <- mu_y - C_t %*% mu_x
	list(mean = mu_t, A = A_t, B = B_t, C = C_t, eigen_values = eigen_values)
}