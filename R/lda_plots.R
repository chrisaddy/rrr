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
	class <- as_data_frame(factor(as_data_frame(class)[[1]]))
	scores_object <- lda_scores(x, class, rank, type, k)
	scores <- scores_object[["scores"]]
	means <- scores_object[["class_means"]]
	scores_x_axis <- scores[,lda_x]
	scores_y_axis <- scores[,lda_y]
	scores_df <- dplyr::bind_cols(scores_x_axis, scores_y_axis, class)
	names(scores_df) <- c("scores_x", "scores_y", "class")
	means_x_axis <- means[,lda_x]
	means_y_axis <- means[,lda_y]
	#means_df <- dplyr::bind_cols(means_x_axis, means_y_axis, )
	ggplot(scores_df, aes(x = scores_x,
				   		  y = scores_y)) + 
		geom_point(aes(color = class)) + 
		#geom_point(means_df, aes(x = ,
		#					    y = as.name(names(means)[2]))) + 
		labs(title = "LD Pairwise Plot",
			 x = paste("LD", lda_x, sep = ""),
			 y = paste("LD", lda_y, sep = ""))
}