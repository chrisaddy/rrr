#` @export

B_t <- function(var_x, var_y, gamma_mat, rank){
		t(V_t(var_x, var_y, gamma_mat, rank)) %*%
			sqrt_matrix(gamma_mat) %*%
			cov_mat(var_y, var_x) %*%
			solve(cov_mat(var_x, var_x))
}
