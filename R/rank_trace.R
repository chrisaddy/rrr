sigma_ee_t <- function(var_x, var_y, gamma_mat, rank){
		e <- var_y - C_t(var_x, var_y, gamma_mat, rank) %*% var_x
		e %*% t(e) / (dim(var_x)[2] - dim(var_x)[1])
}

sigma_ee_full <- function(var_x, var_y, gamma_mat){
       			s <- dim(var_y)[1]
				sigma_ee_t(var_x, var_y, gamma_mat, s)
}

delta_C <- function(x, y, gamma_mat, rank){
		theta <- theta_full(x, y, gamma_mat)
		C <- C_t(x, y, gamma_mat, rank)
		norm(theta - C, type = "F") / norm(theta, type = "F")
}	

delta_EE <- function(x, y, gamma_mat, rank){
		full <- sigma_ee_full(x, y, gamma_mat)	
		t <- sigma_ee_t(x, y, gamma_mat, rank)
		sig_YY <- cov_mat(y, y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}

rank_trace_y <- function(x, y, gamma_mat){
			x_organize 
			rt <- c()
			fr <- min(dim(x)[1], dim(y)[1])
			for(i in 1:fr){
				rt[i] <- delta_EE(x, y, gamma_mat, i)
			}
			c(1,rt)
}

rank_trace_x <- function(x, y, gamma_mat){
			rt <- c()
			fr <- min(dim(x)[1], dim(y)[1])
			for(i in 1:fr){
				rt[i] <- delta_C(x, y, gamma_mat, i)
			}
			c(1,rt)
}

#' Rank Trace
#'
#' \code{rank trace}
#'
#' @param x data frame of input variables
#' @param y data frame of response variables
#' @param gamma_mat weight for \eqn{R} matrix
#' @export rank_trace_x
rank_trace <- function(x, y, gamma_mat) {
				rank <- 0:min(dim(x)[1], dim(y)[1])
				data_frame(rank = rank,
					       deltaC = rank_trace_x(x,
												 y,
												 gamma_mat),
						   deltaEE = rank_trace_y(x,
							  					  y,
							  					  gamma_mat))
}