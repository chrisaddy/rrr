### Rank Trace Plot
###################

delta_coeff <- function(x, y, gamma_matrix, type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	coeff_full <- rrr(x, y, gamma_matrix, rank = full_rank, type, k)$C	
	delta_coeff_norm <- c()
	for(i in 1:full_rank){
		delta_coeff_norm[i] <-sqrt(sum((coeff_full - rrr(x, y, gamma_matrix, i, type, k)$C)^2))
	}
	delta_c <- delta_coeff_norm / sqrt(sum(coeff_full^2))
	data_frame(dC = c(1, delta_c))
}

delta_error <- function(x, y, gamma_matrix, type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	resid_full <- cov(rrr_residuals(x, y, gamma_matrix, full_rank, type, k))
	delta_resid_norm <- c()
	for(i in 1:full_rank){
		delta_resid_norm[i] <- sqrt(sum((resid_full - cov(rrr_residuals(x, y, gamma_matrix, i, type, k)))^2))
	}
	delta_EE <- delta_resid_norm / sqrt(sum((resid_full - cov(y))^2))
	data_frame(dEE = c(1, delta_EE))
}

#' Rank Trace
#'
#' Rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @inheritParams rrr
#'
#' @export

rank_trace <- function(x, y, gamma_matrix, type = "cov", k = 0){
	dC <- delta_coeff(x, y, gamma_matrix, type, k)
	dEE <- delta_error(x, y, gamma_matrix, k)
	ranks <- dplyr::data_frame(ranks = 0:(dim(dC)[1] - 1))
	dplyr::bind_cols(ranks, dC, dEE)

}	

#' Rank Trace Plot
#'
#' Plot of rank trace to determine suitable rank of coefficient matrix.
#'
#' @inheritParams rank_trace
#'
#' @export

rank_trace_plot <- function(x, y, gamma_matrix, type = "cov", k = 0){
	trace <- rank_trace(x, y, gamma_matrix, type = "cov", k)
	ggplot(trace, aes(dC, dEE, label = ranks)) +
	lims(x = c(0,1), y = c(0,1)) +
	geom_line(color = "red") +
	geom_text(check_overlap = TRUE, size = 5) +
	labs(x = "dC", y = "dE") +
	ggtitle(paste("Rank Trace Plot, k =  ", k, sep =""))
}
				 
#' RRR Residuals Plots
#' 
#' \code{rrr_residuals_plot}
#' 
#' @inheritParams rrr
#'
#' @export

rrr_residuals_plot <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	residuals <- rrr_residuals(x, y, gamma_matrix, rank, type, k)
	index <- as_data_frame(index = 1:dim(residuals)[1])
	df <- bind_cols(index, residuals)
	GGally::ggpairs(df)
}