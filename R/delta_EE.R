delta_EE <- function(var_x, var_y, gamma_mat, rank){
		full <- sigma_ee_full(var_x, var_y, gamma_mat)	
		t <- sigma_ee_t(var_x, var_y, gamma_mat, rank)
		sig_YY <- cov_mat(var_y, var_y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}