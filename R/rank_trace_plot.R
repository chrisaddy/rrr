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