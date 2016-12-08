### Rank Trace Plot
###################

delta_coeff <- function(x, y, gamma_matrix, type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	coeff_full <- rrr(x, y, gamma_matrix, rank = full_rank, type, k)[["C"]]	
	delta_coeff_norm <- c()
	for(i in 1:full_rank){
		delta_coeff_norm[i] <-sqrt(sum((coeff_full - rrr(x, y, gamma_matrix, i, type, k)[["C"]])^2))
	}
	delta_c <- delta_coeff_norm / sqrt(sum(coeff_full^2))
	data_frame(dC = c(1, delta_c))
}

delta_error <- function(x, y, gamma_matrix, type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	resid_full <- cov(rrr_residual(x, y, gamma_matrix, full_rank, type, k))
	delta_resid_norm <- c()
	for(i in 1:full_rank){
		delta_resid_norm[i] <- sqrt(sum((resid_full - cov(rrr_residual(x, y, gamma_matrix, i, type, k)))^2))
	}
	delta_EE <- delta_resid_norm / sqrt(sum((resid_full - cov(y))^2))
	data_frame(dEE = c(1, delta_EE))
}

#' Rank Trace Plot
#'
#' \code{rank_trace} is a plot used to determine the effective dimensionality, i.e., \eqn{t = \mathrm{rank}\left(\mathbf{C}\right)},
#' of the reduced-rank regression equation.
#'
#' @inheritParams rrr
#' @param plot if FALSE, returns data frame of rank trace coordinates. 
#' @param interactive if TRUE, creates an interactive plotly graphic.
#'
#' @return plot of rank trace coordinates if \code{plot = TRUE} or data frame of rank trace coordinates if \code{plot = FALSE}.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' rank_trace(tobacco_x, tobacco_y, gamma)
#' rank_trace(tobacco_x, tobacco_y, gamma, plot = FALSE)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rank_trace <- function(x, y, gamma_matrix, type = "cov", k = 0, plot = TRUE, interactive = FALSE){
	dC <- delta_coeff(x, y, gamma_matrix, type, k)
	dEE <- delta_error(x, y, gamma_matrix, k)
	ranks <- dplyr::data_frame(ranks = 0:(dim(dC)[1] - 1))
	trace <- dplyr::bind_cols(ranks, dC, dEE)
	if(plot == FALSE){
		trace
	} else {
	rt_plot <- ggplot(trace, aes(dC, dEE, label = ranks)) +
		lims(x = c(0,1), y = c(0,1)) +
		geom_line(color = "red") +
		geom_point(size = 5) + 
		geom_text(check_overlap = TRUE, size = 4, color = "white") +
		labs(x = "dC", y = "dE") +
		ggtitle(paste("Rank Trace Plot, k =  ", k, sep =""))
		if(interactive == TRUE){
			plotly::ggplotly(rt_plot)
			} else {
				rt_plot
			}
		}
	}
				 
#' Plot Residuals of Reduced-Rank Regression
#' 
#' \code{rrr_residual_plot} is a scatter plot matrix used for diagnostics of the reduced-rank regression model.
#' @inheritParams rrr
#'
#' @return ggplot object, scatterplot matrix.
#'
#' @examples
#' data(tobacco)
#' tobacco_x <- tobacco[,4:9]
#' tobacco_y <- tobacco[,1:3]
#' gamma <- diag(1, dim(tobacco_y)[2])
#' rrr_residual_plot(tobacco_x, tobacco_y, gamma)
#' gamma2 <- solve(cov(tobacco_y))
#' rrr_residual_plot(tobacco_x, tobacco_y, gamma2)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrr_residual_plot <- function(x, y, gamma_matrix, rank = "full", type = "cov", k = 0){
	residuals <- rrr_residual(x, y, gamma_matrix, rank, type, k)
	index <- dplyr::as_data_frame(1:dim(residuals)[1])
	names(index) <- "index"
	df <- bind_cols(index, residuals)
	plot <- GGally::ggpairs(df)
	plot[1,1] <- ggplot()
	plot
}