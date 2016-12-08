#' LDA Pairwise Plot
#'
#' \code{lda_pairwise_plot} creates a plot with the linear discriminant scores of one linear discriminant
#' function along the \eqn{X}-axis against the linear discriminant scores of another linear discriminant function
#' alonge the \code{Y}-axis.
#'
#' @inheritParams lda_scores
#' @inheritParams pca_pairwise_plot
#' @param lda_x integer. Linear discriminant function plotted along the x-axis.
#' @param lda_y integer. Linear discriminant function plotted along the y-axis.
#'
#' @return ggplot object if \code{interactive = FALSE}, the default, or an interactive plotly plot if \code{interactive = TRUE}.
#'
#' @examples
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' lda_pairwise_plot(iris_x, iris_y)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

lda_pairwise_plot <- function(x, class, lda_x = 1, lda_y = 2, rank = "full", type = "cov", k = 0, quadratic = FALSE, interactive = FALSE){
	#class <- as_data_frame(factor(as_data_frame(class)[[1]]))
	scores_x <- scores_y <- NULL
	scores_object <- lda_scores(x, class, rank, type, k, quadratic)
	scores <- scores_object[["scores"]]
	scores_x_axis <- scores[,lda_x]
	scores_y_axis <- scores[,lda_y]
	scores_class <- scores[, dim(scores)[2]]
	scores_df <- dplyr::bind_cols(scores_x_axis, scores_y_axis, scores_class)
	names(scores_df) <- c("scores_x", "scores_y", "class")
	means <- scores_object[["class_means"]]
	means_x_axis <- means[,lda_x]
	means_y_axis <- means[,lda_y]
	means_class <- means[, dim(means)[2]]
	means_df <- dplyr::bind_cols(means_x_axis, means_y_axis, means_class)
	names(means_df) <- c("means_x", "means_y", "class")
	static_plot <- return(ggplot() +
		geom_point(data = scores_df, aes(x = scores_x, y = scores_y, color = class)) + 
		#geom_text(data = means_df, aes(x = means_x,
		#					     	    y = means_y,
		#					     	    label = abbreviate(class))) + 
		labs(title = "LD Pairwise Plot",
			 x = paste("LD", lda_x, sep = ""),
			 y = paste("LD", lda_y, sep = "")))
	if(interactive == TRUE){
		plotly::ggplotly(static_plot)
	} else {
		static_plot
	}
}

#' LDA 3D Plot
#'
#' \code{lda_3D_plot} creates an interactive, html plotly plot that can be manipulated by the viewer. 
#'
#' @inheritParams lda_pairwise_plot
#' @inheritParams pca_3D_plot
#' @param lda_z integer. Linear discriminant function plotted along the z-axis.
#'
#' @return plotly object.
#' 
#' @export

lda_3D_plot <- function(x, class, lda_x = 1, lda_y = 2, lda_z = 3, rank = "full", type = "cov", k = 0, quadratic = FALSE, point_size = 3){
	class <- as_data_frame(factor(as_data_frame(class)[[1]]))
	scores <- lda_scores(x, class, rank, type, k, quadratic)[["scores"]]
	if(dim(scores)[2] <= 4){
		stop("too few linear discriminant functions for lda_3D_plot")
	}
	scores_x <- scores[,lda_x]
	scores_y <- scores[,lda_y]
	scores_z <- scores[,lda_z]
	scores_tbl <- dplyr::bind_cols(scores_x, scores_y, scores_z, class)
	names(scores_tbl) <- c("x_coord", "y_coord", "z_coord", "class")
	ptsize <- point_size
	plotly::plot_ly(scores_tbl,
					x = ~x_coord,
					y = ~y_coord,
					z = ~z_coord,
					type = "scatter3d",
					mode = "markers",
					color = ~class,
					marker = list(size = ptsize),
					name = "LDA 3D Scatter Plot")
}