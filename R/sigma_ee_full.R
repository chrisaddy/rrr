#` @export

sigma_ee_full <- function(var_x, var_y, gamma_mat){
       			s <- dim(var_y)[1]
			sigma_ee_t(var_x, var_y, gamma_mat, s)
}