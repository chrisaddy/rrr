#` @export

rank_trace_x <- function(var_x, var_y, gamma_mat){
			rt <- c()
			fr <- min(dim(var_x)[1], dim(var_y)[1])
			for(i in 1:fr){
				rt[i] <- delta_C(var_x, var_y, gamma_mat, i)
			}
			c(1,rt)
}