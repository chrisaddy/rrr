#' @importFrom magrittr %>%
#' @export
magrittr::'%>%'

#' @importFrom dplyr as_data_frame bind_cols
#' @export

#' @importFrom stats cov
#' @export

#' @importFrom MASS ginv
#' @export
MASS::ginv

sqrt_matrix <- function(matr){
			eigens <- eigen(matr)
			vecs <- eigens$vectors
			vals <- eigens$values
			vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs)
}

organize <- function(vars, type = "cov"){
	matr <- as.matrix(vars)
	if(type == "cov"){
		matr %>%
		scale(scale = FALSE) %>%
		t()
		} else if(type == "cor"){
			matr %>%
			scale() %>%
			t()
		} else {
			stop("type not recognized")
		}
}

cov_matrix <- function(x, y, type = "cov"){
	n <- dim(x)[1]
	x_org <- organize(x, type)
	y_org <- organize(y, type)
	x %*% t(y) / n
}

binary_matrix <- function(class) {
    class <- as.matrix(dplyr::mutate_if(class, is.factor, as.character))
    mat <- stats::model.matrix(~ class -1)
##    mat <- dplyr::select(mat, -dim(mat)[2])
##    names(mat) <- substring(names(mat), 4, 100)
   # mat[,dim(mat)[2]] <- 0
    mat
}