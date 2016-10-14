mu_vars <- function(var_matrix){
				replicate(dim(var_matrix)[2],
				rowMeans(var_matrix))
}
