mu_t <- function(var_x, var_y, gamma_mat, rank){
		mu_vars(var_y) - 
			C_t(var_x, var_y, gamma_mat, rank) %*% 
				mu_vars(var_x)
}