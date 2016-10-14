#` @export

theta_full <- function(var_x, var_y, gamma_mat){
		s <- dim(var_y)[1]
		C_t(var_x, var_y, gamma_mat, s)
}