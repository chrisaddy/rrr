sigma_ee_t <- function(var_x, var_y, gamma_mat, rank){
		e <- var_y - C_t(var_x, var_y, gamma_mat, rank) %*% var_x
		e %*% t(e) / (dim(var_x)[2] - dim(var_x)[1])
}