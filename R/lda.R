# \code{lda} fits a reduced-rank linear discriminant analysis model by fitting a reduced-rank canonical variate
# model to the predictor variables and a response-indicator matrix.
# 
# @inheritParams rrr
# @param class vector or one-column data frame of type factor or character that are the class labels of the observations.
# 
# @return list containing: a named vector of prior probabilities; a data frame of canonical-variate coefficients of the predictor variables;
# a data frame of canonical-variate coefficients for the reponse variables.


lda <- function(x, class, rank = "full", k = 0){
	ordered <- lda_organize(x, class)
	x_organize <- ordered[["features_ordered"]]
	y_organize <- ordered[["classes_ordered"]]
    full_rank <- min(dim(x_organize)[2], dim(y_organize)[2])
    if(rank == "full"){
    reduce_rank <- full_rank
    } else if(rank <= full_rank){
    reduce_rank <- rank
    } else {
    stop("rank out of bounds")
    }
	cva_object <- cva(x_organize, y_organize, rank, k)
    G <- as_data_frame(t(cva_object$G))
    names(G) <- paste("LD", 1:reduce_rank, sep = "")
    H <- as_data_frame(cva_object$H)
    names(H) <- names(G)
    class_names <- class %>% 
        as_data_frame() %>%
        select(class = 1) %>%
        mutate(class = as.character(class)) %>%
        distinct() %>%
        arrange(class) %>%
        as.matrix() %>% 
        as.vector()
    n <- colSums(y_organize)
    num_classes <- dim(y_organize)[1]
    mean_y <- n / num_classes
    prior <- c(mean_y, 1 - sum(mean_y))
    names(prior) <- class_names
    list(G = G, H = H)
}

# \code{lda_scores} calculates canonical-variate scores. In the case of linear discriminant analysis, the canonical-variate
# scores of the predictor variables are the linear discriminant functions evaluated at the values of predictor variables, while the scores for the
# response variables correspond to class means.
#
# @inheritParams lda
#
# @return list containing: data frame of linear discriminant function scores; data frame of class means.
#
# @examples
# data(iris)
# lda_scores(iris[,1:4], iris[,5])

lda_scores <- function(x, class, rank = "full", k = 0){
    class_label <- as_data_frame(factor(as_data_frame(class)[[1]]))
    names(class_label) <- "class"
    lda_object <- lda(x, class_label, rank, k)
    ordered <- lda_organize(x, class_label)
    x_organize <- ordered[["features_ordered"]]
    y_organize <- ordered[["classes_ordered"]]
    xi <- as_data_frame(t(t(as.matrix(lda_object[["G"]])) %*% t(x_organize)))
    names(xi) <- paste("LD", 1:dim(xi)[2], sep = "")
    omega <- as_data_frame(t(as.matrix(lda_object[["H"]]) %*% t(y_organize)))
    names(omega) <- names(xi)
    list(scores = dplyr::bind_cols(xi, class_label), class_means = dplyr::distinct(dplyr::bind_cols(omega, class_label))) 
}

# \code{lda_pairwise_plot} creates a plot with the linear discriminant scores of one linear discriminant
# function along the \eqn{X}-axis against the linear discriminant scores of another linear discriminant function
# alonge the \code{Y}-axis.
#
# @return ggplot object if \code{interactive = FALSE}, the default, or an interactive plotly plot if \code{interactive = TRUE}.
#
# @examples
# data(iris)
# iris_x <- iris[,1:4]
# iris_y <- iris[5]
# lda_pairwise_plot(iris_x, iris_y)

lda_pairwise_plot <- function(x, class, lda_x = 1, lda_y = 2, rank = "full", k = 0, interactive = FALSE){
    #class <- as_data_frame(factor(as_data_frame(class)[[1]]))
    scores_x <- scores_y <- NULL
    scores_object <- lda_scores(x, class, rank, k)
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
        #                               y = means_y,
        #                               label = abbreviate(class))) + 
        labs(title = "LD Pairwise Plot",
             x = paste("LD", lda_x, sep = ""),
             y = paste("LD", lda_y, sep = "")))
    if(interactive == TRUE){
        plotly::ggplotly(static_plot)
    } else {
        static_plot
    }
}

# LDA 3D Plot
#
# \code{lda_3D_plot} creates an interactive, html plotly plot that can be manipulated by the viewer. 
#
# @inheritParams lda_pairwise_plot
# @inheritParams pca_3D_plot
# @param lda_z integer. Linear discriminant function plotted along the z-axis.
#
# @examples
# data(pendigits)
# digits_features <- pendigits[,1:34]
# digits_class <- pendigits[,35]
# lda_3D_plot(digits_features, digits_class, k = 0.0001)
#
# @return plotly object.
# 
# @export

lda_3D_plot <- function(x, class, lda_x = 1, lda_y = 2, lda_z = 3, rank = "full", k = 0, point_size = 3){
    #x_coord <- y_coord <- z_coord <- NULL
    class <- as_data_frame(factor(as_data_frame(class)[[1]]))
    scores <- lda_scores(x, class, rank, k)[["scores"]]
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
                    marker = list(size = ptsize))
}