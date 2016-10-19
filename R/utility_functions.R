cov_matrix <- function(var1, var2){
		var1 %*% t(var2) / (dim(var1)[2] - 1)
}

sqrt_matrix <- function(matr){
			e <- eigen(matr)
			vecs <- e$vectors
			vals <- e$values
			vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs)
}

organize <- function(vars, scale = FALSE){
				vars %>%
				as.matrix() %>%
				scale(center = TRUE, scale) %>%
				t()
}

mu_vars <- function(var){
				replicate(dim(var)[2],
				rowMeans(var))
}