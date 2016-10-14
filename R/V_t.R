V_t <- function(var_x, var_y, gamma_mat, rank){
	       		eigen(weighted_mat(var_x,
						         						   var_y,
													          						   gamma_mat))$vectors[,1:rank]
}
