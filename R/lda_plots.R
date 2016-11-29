#' LD Pairwise Plot
#'
#' \code{lda_pairwise_plot}
#'
#' @inheritParams lda_scores
#' @param lda_x integer. Linear discriminant function plotted along the x-axis.
#' @param lda_y integer. Linear discriminant function plotted along the y-axis.
#'
#' @export

lda_pairwise_plot <- function(x, class, lda_x = 1, lda_y = 2, rank = "full", type = "cov", k = 0){
	class <- as_data_frame(class)
	scores_object <- lda_scores(x, class, rank, type, k)
	scores <- scores_object[["scores"]]
	means <- scores_object[["class_means"]]
	scores_x_axis <- scores[,lda_x]
	scores_y_axis <- scores[,lda_y]
	means_x_axis <- means[,lda_x]
	means_y_axis <- means[,lda_y]
	df <- dplyr::bind_cols(scores_x_axis, scores_y_axis, class)
	#ggplot(df, aes_q(x = as.name(names(df)[1]),
	#				 y = as.name(names(df)[2]))) + 
	#	geom_point() + 
		#geom_point(means, aes_q(x = as.name(names(means)[1]),
		#					    y = as.name(names(means)[2]))) + 
	#	labs(title = "LD Pairwise Plot")
}