#` A Matrix
#`
#` @param var_x matrix of independent variables
#` @param var_y matrix of response variables
#` @param gamma_mat weight matrix
#` @param rank rank of the coefficient matrix
#`
#` @return None
#` @export

A_t <- function(var_x, var_y, gamma_mat, rank){
		solve(sqrt_matrix(gamma_mat)) %*% V_t(var_x,
											  var_y,
											  gamma_mat,
											  rank)
}