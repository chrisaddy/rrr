#' Square Root Matrix
#'
#' \code{sqrt_matrix} calculates the square root of a positive semi-definite matrix through spectral decomposition.
#'
#' @param matr a positive semi-definite matrix
#'
#' @export

sqrt_matrix <- function(matr){
			eigens <- eigen(matr)
			vecs <- eigens$vectors
			vals <- eigens$values
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
	as.matrix(vars) %>%
		scale(center = TRUE, scale) %>%
		t()
}

binary_matrix <- function(class) {
    class <- as.matrix(dplyr::mutate_if(class, is.factor, as.character))
    mat <- stats::model.matrix(~ class -1)
##    mat <- dplyr::select(mat, -dim(mat)[2])
##    names(mat) <- substring(names(mat), 4, 100)
   # mat[,dim(mat)[2]] <- 0
    mat
}