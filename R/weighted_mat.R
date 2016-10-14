#` @export

weighted_mat <- function(var_x, var_y, gamma_mat){ 
			sig_XX <- cov_mat(var_x, var_x)
			sig_XY <- cov_mat(var_x, var_y)
			sig_YX <- cov_mat(var_y, var_x) 
			sqrt_matrix(gamma_mat) %*%
				sig_YX %*%
				solve(sig_XX) %*%
				sig_XY %*%
				sqrt_matrix(gamma_mat)
}
