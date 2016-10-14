A_t <- function(var_x, var_y, gamma_mat, rank){
		solve(sqrt_matrix(gamma_mat)) %*% V_t(var_x,
											  var_y,
											  gamma_mat,
											  rank)
}