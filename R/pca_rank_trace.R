delta_C <- function(var_x, var_y, gamma_mat, rank){
		theta <- theta_full(var_x, var_y, gamma_mat)
		C <- C_t(var_x, var_y, gamma_mat, rank)
		norm(theta - C, type = "F") / norm(theta, type = "F")
}

delta_EE <- function(var_x, var_y, gamma_mat, rank){
		full <- sigma_ee_full(var_x, var_y, gamma_mat)	
		t <- sigma_ee_t(var_x, var_y, gamma_mat, rank)
		sig_YY <- cov_mat(var_y, var_y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}

#` Rank Trace for PCA
#`
#` \code{pca_rank_trace} is a wrapper function for \code{rank_trace} for Principle Components Analysis
#`
#` @param x
#`
#` export pca_rank_trace

pca_rank_trace <- function(x) {
					rank_trace(x, x, diag(1, dim(x)[1]))
}