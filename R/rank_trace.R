theta_full <- function(var_x, var_y, gamma_mat){
		s <- dim(var_y)[1]
		C_t(var_x, var_y, gamma_mat, s)
}

sigma_ee_t <- function(x, y, gamma_mat, rank){
		e <- y - rrr(x, y, gamma_mat, rank)$C %*% x
		e %*% t(e) / (dim(x)[2] - dim(x)[1])
}

sigma_ee_full <- function(x, y, gamma_mat){
       			s <- dim(y)[1]
				sigma_ee_t(x, y, gamma_mat, s)
}

delta_C <- function(x, y, gamma_mat, rank){
		theta <- theta_full(x, y, gamma_mat)
		C <- C_t(x, y, gamma_mat, rank)
		norm(theta - C, type = "F") / norm(theta, type = "F")
}	

delta_EE <- function(x, y, gamma_mat, rank){
		full <- sigma_ee_full(x, y, gamma_mat)	
		t <- sigma_ee_t(x, y, gamma_mat, rank)
		sig_YY <- cov_matrix(y, y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}


#' Rank Trace
#'
#' \code{rank trace}
#'
#' @param x data frame of input variables
#' @param y data frame of response variables
#' @param gamma_mat weight for \eqn{R} matrix
#' @export rank_trace

rank_trace <- function(x, y, gamma_mat) {
				x_org <- organize(x)
				y_org <- organize(y)
				full_rank <- min(dim(x_org)[1], dim(y_org)[1])
				rt_x <- c()
				rt_y <- c()
				for(i in 1:full_rank){
					rt_x[i] <- delta_C(x_org, y_org, gamma_mat, i)
					rt_y[i] <- delta_EE(x_org, y_org, gamma_mat, i)
				}
				rt_x <- c(1, rt_x)
				rt_y <- c(1, rt_y)
				rank <- 0:min(dim(x)[2], dim(y)[2])
				data_frame(rank = rank,
					       deltaC = rt_x,
						   deltaEE = rt_y)
}

#' Rank Trace Plot
#'
#' Plot of rank trace to determine suitable rank of coefficient matrix.
#'
#' @param x data frame of input variables
#' @param y data frame of response variables
#' @param gamma_mat weight matrix
#'
#' @export

rank_trace_plot <- function(x, y, gamma_mat){
			trace <- rank_trace(x, y, gamma_mat)
			ggplot() +
			aes(x = trace$deltaC,
				y = trace$deltaEE,
				label = trace$rank) +
			lims(x = c(0,1), y = c(0,1)) +
			geom_line(color = "red") +
			geom_text(check_overlap = TRUE, size = 5) +
			geom_label() + 
			labs(x = "dC", y = "dE") +
			ggtitle(expression(paste("Rank Trace Plot for ",
						 Theta^(0),
						 " to ",
						 Theta^(s) )))
}