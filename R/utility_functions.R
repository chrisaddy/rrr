cov_matrix <- function(var1, var2){
		var1 %*% t(var2) / (dim(var1)[2] - 1)
}

#' Square Root Matrix
#'
#' \code{sqrt_matrix} calculates the square root of a positive semi-definite matrix through spectral decomposition.
#'
#' @param matr a positive semi-definite matrix
#'
#' @export

sqrt_matrix <- function(matr){
			e <- eigen(matr)
			vecs <- e$vectors
			vals <- e$values
			vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs)
}

#' Organize Data for RRR
#'
#' \code{organize} mean-center and transpose a matrix
#'
#' @param vars data frame or matrix
#' @param scale logical. If \code{TRUE}, matrix is scaled by standard deviations.
#'
#' @export
organize <- function(vars, scale = FALSE){
	t(scale(as.matrix(vars), center = TRUE, scale))
}

#mu_vars <- function(var){
#				replicate(dim(var)[2],
#				rowMeans(var))
#}
