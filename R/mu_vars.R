#' Mean Matrix
#'
# \code{mu_vars} creates a matrix of sample means from a matrix \eqn{X} with the same dimensions as \eqn{X}. 

mu_vars <- function(var_matrix){
				replicate(dim(var_matrix)[2],
				rowMeans(var_matrix))
}
