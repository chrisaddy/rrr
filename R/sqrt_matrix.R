sqrt_matrix <- function(matr){
			e <- eigen(matr)
			vecs <- e$vectors
			vals <- e$values
			vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs)
}
