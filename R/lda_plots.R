#' LD Pairwise Plot
#'
#' \code{ld_pairwise_plot}
#'
#' @inheritParams ld_scores
#' @param ld_x integer. Linear discriminant function plotted along the x-axis.
#' @param ld_y integer. Linear discriminant function plotted along the y-axis.
#'
#' @export

ld_pairwise_plot <- function(x, class, ld_x = 1, ld_y = 2, rank = "full", type = "cov", k = 0){
	scores_object <- ld_scores(x, class, rank, type, k)[["scores"]]
	x_axis <- scores_object[,ld_x]
	y_axis <- scores_object[,ld_y]
	df <- dplyr::bind_cols(x_axis, y_axis, class)
	ggplot(df, aes_q(x = as.name(names(df)[1]),
					 y = as.name(names(df)[2]))) + 
		geom_point() + 
		labs(title = "LD Pairwise Plot")
}