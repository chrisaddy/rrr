#` @export

delta_C <- function(var_x, var_y, gamma_mat, rank){
		theta <- theta_full(var_x, var_y, gamma_mat)
		C <- C_t(var_x, var_y, gamma_mat, rank)
		norm(theta - C, type = "F") / norm(theta, type = "F")
}