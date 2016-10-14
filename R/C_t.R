#` @export

C_t <- function(var_x, var_y, gamma_mat, rank) {
		A_t(var_x,
			var_y,
			gamma_mat,
			rank) %*% B_t(var_x,
					   	  var_y,
						  gamma_mat,
						  rank)
}