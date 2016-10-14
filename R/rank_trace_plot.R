#` @export

rank_trace_plot <- function(var_x, var_y, gamma_mat){
			rx <- rank_trace_x(var_x, var_y, gamma_mat)
			ry <- rank_trace_y(var_x, var_y, gamma_mat)
			ggplot() +
			aes(x = rx, y = ry) +
			lims(x = c(0,1), y = c(0,1)) +
			geom_point() +
			geom_line(color = "red") +
			labs(x = "dC", y = "dE") +
			ggtitle(expression(paste("Rank Trace Plot for ",
						 Theta^(0),
						 " to ",
						 Theta^(s) )))
}