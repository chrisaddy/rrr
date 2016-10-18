sigma_ee_t <- function(x, y, gamma_mat, rank){
		e <- y - C_t(x_organize, y_organize, gamma_mat, rank) %*% x
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
		sig_YY <- cov_mat(y, y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}

rank_trace_y <- function(x, y, gamma_mat){
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
				x_organize <- organize(x)
				y_organize <- organize(y)
				rank <- 0:min(dim(x)[2], dim(y)[2])
				data_frame(rank = rank,
					       deltaC = rank_trace_x(x_organize,
												 y_organize,
												 gamma_mat),
						   deltaEE = rank_trace_y(x_organize,
							  					  y_organize,
							  					  gamma_mat))
}

#' Rank Trace Plot
#
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