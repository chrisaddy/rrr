#' Reduced-Rank PCA
#' 
#' \code{pca} carries out reduced-rank principled component analysis on a
#' matrix of input data.
#' 
#' 
#' @export pca
pca <- function(x, rank, ridge = 0){
		x_ridge <- x
		gamma <- diag(1, dim(x)[1])
		mu <- mu_t(x_ridge, x_ridge, gamma, rank)
		coefficients <- C_t(x_ridge, x_ridge, gamma, rank)
		estimates <- mu + coefficients %*% x_ridge
		residuals <- estimates - x_ridge
		list(coefficients, estimates, residuals)
}
