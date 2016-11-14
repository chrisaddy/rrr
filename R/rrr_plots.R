delta_EE <- function(x, y, gamma_mat, rank){
	full_rank <- min(dim(x)[2], dim(y)[2])
	y_organize <- organize(y)
	x_organize <- organize(x)
	sig_YX <- cov_matrix(y_organize, x_organize)
	sig_XX <- cov_matrix(x_organize, x_organize)
	sig_YY <- cov_matrix(y_organize, y_organize)
	sigma_ee_t <- y_organize - rrr(x, y, gamma_mat, rank)$C %*% x_organize
	e_full <- rrr_residuals(x, y, gamma_mat, full_rank)
	sigma_ee_full <- e_full %*% t(e_full)
	e_t <- rrr_residuals(x, y, gamma_mat, rank)
	sigma_ee_t <- e_t %*% t(e_t)
	sum((sigma_ee - sigma_ee_t)^2) / sum((sigma_ee - sig_YY)^2)
}

#' @export

delta_C <- function(x, y, gamma_mat, type = "cov", k = 0){
	theta_full <- rrr(x, y, gamma_mat, rank = "full", type, k)$C
	reduced_rank <- min(dim(x)[2], dim(y)[2])
	delta <- c()
	for(i in 1:reduced_rank){
		delta[i] <- sum((theta_full - rrr(x, y, gamma_mat, rank = i, type, k)$C)^2) / sum(theta_full^2)
	}
	c(1, delta)
}

#' @export

delta_EE <- function(x, y, gamma_mat, type = "cov", k = 0){
	cov_y <- cov(y) + k * diag(1, dim(y_organize)[1])
	cov_EE_full <- cov_y -  rrr(x, y, gamma_mat, rank = "full", type, k)
	reduced_rank <- min(dim(x)[2], dim(y)[2])
#	for(i in 1:reduced_rank){
#		sum(cov_EE_full - ) / sum((cov_EE_full - cov_y)^2)
#	}
}

#' Rank Trace
#'
#' Rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @inheritParams rrr
#'
#' @export

rank_trace <- function(x, y, gamma_mat) {
	s <- min(dim(x)[2], dim(y)[2])
	rt_x <- c()
	rt_y <- c()
	for(i in 1:s){
		rt_x[i] <- delta_C(x, y, gamma_mat, i)
		rt_y[i] <- delta_EE(x, y, gamma_mat, i)
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
#' @inheritParams rrr
#' @param gamma_mat weight matrix
#'
#' @export

rank_trace_plot <- function(x, y, gamma_mat){
	trace <- rank_trace(x, y, gamma_mat)
	ggplot() +
	aes(x = rank_trace_x(x, y, gamma_mat), #trace$deltaC,
		y = rank_trace_y(x, y, gamma_mat), #trace$deltaEE,
		label = 0:3) +
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