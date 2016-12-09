pca <- function(x, rank = "full", k = 0){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    s_xx <- cov(x) + k * diag(1, dim(x)[2])
    eigens <- eigen(s_xx)
    A <- eigens[["vectors"]][,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    vals <- eigens[["values"]]
    total_var <- sum(vals)
    gof <- c()
    for(i in 1:length(vals)){
        gof[i] <- sum(vals[(i + 1):length(vals)]) / total_var
    }
    gof[length(gof)] <- 0
    gof <- gof[1:reduce_rank]
    names(gof) <- colnames(A)
    list(means = means, C = as_data_frame(A %*% t(A)), PC = as_data_frame(A),goodness_of_fit = gof)
}

pca_scores <- function(x, rank = "full", k = 0){
    pca <- pca(x, rank, k)
    scores <- t(t(as.matrix(pca[["PC"]])) %*% organize(x)) %>%
        as_data_frame()
    names(scores) <- paste("PC", 1:dim(scores)[2], sep = "")
    scores
}

pca_rank_trace <- function(x, k = 0, plot = TRUE, interactive = FALSE){
    eigenvecs <- pca(x, rank = "full", k)$A
    eigens <- eigen(cov(x) + k * diag(1, dim(x)[2]))
    full_rank <- dim(eigens$vectors)[2]
    delta_C <- function(t){
        sqrt(1 - t / full_rank)
    }
    delta_residuals <- function(t){
        sqrt((sum(eigens$values[(t+1):full_rank])^2) / sum((eigens$values))^2)
    }
    pca_rank_trace_x <- c()
    pca_rank_trace_y <- c()
    for(i in 1:full_rank){
        pca_rank_trace_x[i] <- delta_C(i)
        pca_rank_trace_y[i] <- delta_residuals(i)
    }
    pca_rank_trace_y[length(pca_rank_trace_y)] <- 0
    trace <- dplyr::data_frame(rank = 0:length(pca_rank_trace_x),
        delta_C = c(1, pca_rank_trace_x),
        delta_residuals = c(1, pca_rank_trace_y))
    if(plot == FALSE){
        trace
    } else {
        tuner <- k
        static_plot <- ggplot2::ggplot(trace, ggplot2::aes(x = delta_C,
                y = delta_residuals,
                label = rank)) +
            lims(x = c(0,1), y = c(0,1)) +
            geom_line(color = "red") +
            geom_point(size = 5) + 
            geom_text(check_overlap = TRUE, size = 4, color = "white") + 
            labs(x = "dC", y = "dE") +
            ggtitle(paste("PCA Rank Trace Plot, k = ", tuner, sep = ""))
        if(interactive == TRUE){
                plotly::ggplotly(static_plot )
            } else {
                static_plot
            }
    }
}

pca_pairwise <- function(x, pca_x, pca_y, rank = "full", k = 0){
    scores <- pca_scores(x, rank, k)
    score_1 <- paste("PC", pca_x, sep = "")
    score_2 <- paste("PC", pca_y, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2))
}

# \code{pca_pairwise_plot} plots the scores of one principal component along the \eqn{X}-axis 
# against the scores of another along the \eqn{Y}-axis in a two-dimensional scatterplot.
# @return ggplot object if \code{interactive = FALSE}, the default, or an interactive plotly plot if \code{interactive = TRUE}.
# @examples
# data(pendigits)
# digits_features <- pendigits[,1:34]
# digits_class <- pendigits[,35]
# pca_pairwise_plot(digits_features, pc_x = 1, pc_y = 3)
#
 
pca_pairwise_plot <- function(x, pc_x = 1, pc_y = 2, class_labels = NULL, rank = "full", k = 0, interactive = FALSE, point_size = 2.5){
    pairs <- pca_pairwise(x, pc_x, pc_y, rank, k)
    if(is.null(class_labels)){
        pairs_tbl <- pairs
        names(pairs_tbl) <- c("pc_x", "pc_y")
        static_plot <- ggplot2::ggplot(pairs_tbl,
                    aes(pc_x,
                        pc_y)) +
            geom_point() +
            labs(x = paste("PC", pc_x, sep = ""),
                          y = paste("PC", pc_y, sep = "")) +
            labs(title = "PC Pairwise") +
            theme(legend.title = element_blank())
    } else {
        class <- dplyr::as_data_frame(class_labels)
        pairs_tbl <- bind_cols(pairs, class)
        names(pairs_tbl) <- c("pc_x", "pc_y", "class")
        ptsize <- point_size
        static_plot <- ggplot2::ggplot(pairs_tbl,
                                       aes(pc_x,
                                           pc_y)) +
            geom_point(aes(color = factor(class)), size = (ptsize - 1.5)) +
            labs(x = paste("PC", pc_x, sep = ""),
                          y = paste("PC", pc_y, sep = "")) +
            labs(title = "PC Pairwise") +
            theme(legend.title = element_blank())
    }
    if(interactive == TRUE){
        interactive_plot <- plotly::ggplotly(static_plot)
        for(i in 1:dim(distinct(class_labels))[1]){
            interactive_plot[["x"]][["data"]][[i]][["marker"]][["size"]] <- ptsize
        }
        interactive_plot
    } else {
        print(static_plot)
    }
}

pca_threewise <- function(x, pca_x, pca_y, pca_z, rank = "full"){
    scores <- pca_scores(x, rank)
    score_1 <- paste("PC", pca_x, sep = "")
    score_2 <- paste("PC", pca_y, sep = "")
    score_3 <- paste("PC", pca_z, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2), ends_with(score_3))
}


# \code{pca_3D_plot} creates an interactive, html plotly plot that can be manipulated by the viewer.
#
# @inheritParams pca
# @param pca_x integer number of the principal component used for the x-axis
# @param pca_y integer number of the principal component used for the y-axis
# @param pca_z integer number of the principal component used for the z-axis
# @param class_labels data frame or vector of class labels
# @param rank rank of coefficient matrix
# @param point_size size of points in scatter
#
#@return plotly object.

pca_3D_plot <- function(x, class_labels = "none", pca_x = 1, pca_y = 2, pca_z = 3, rank = "full", point_size = 3){
    if(class_labels == "none"){
        class <- dplyr::as_data_frame(rep("class = none", dim(x)[1]))
    } else {
        class <- dplyr::as_data_frame(class_labels)
    }
    threes <- pca_threewise(x, pca_x, pca_y, pca_z, rank)
    threes_tbl <- bind_cols(threes, class)
    names(threes_tbl) <- c("x_coord", "y_coord", "z_coord", "class")
    ptsize <- point_size
    plotly::plot_ly(threes_tbl,
        x = ~x_coord,
        y = ~y_coord,
        z = ~z_coord,
        type = "scatter3d",
        mode = "markers",
        color = ~factor(class),
        marker = list(size = ptsize))
}

# plots all pairs of principal components
#
# @inheritParams pca_scores
# @inheritParams pca_pairwise_plot

pca_allpairs_plot <- function(x, rank, k = 0, class_labels = NULL){
    pc <- pca_scores(x, rank, k)
    if(is.null(class_labels)){
    GGally::ggpairs(pc, upper = list(continuouse = "blank"))
    } else {
    class <- as_data_frame(as.factor(class_labels[[1]]))
    df <- dplyr::bind_cols(pc, class)
    names(df)[dim(df)[2]] <- "class"
    GGally::ggpairs(df, upper = list(continuous = "blank"), aes(color = class))
    }
}