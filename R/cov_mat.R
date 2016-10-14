cov_mat <- function(var1, var2){
		var1 %*% t(var2) / (dim(var1)[2] - 1)
}
