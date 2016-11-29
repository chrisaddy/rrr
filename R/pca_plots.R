#' PCA Rank Trace
#'
#' Rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @inheritParams pca
#' 
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[,1:34]
#' pca_rank_trace(digits_features, k = 0.001)
#'
#' @export

pca_rank_trace <- function(x, type = "cov", k = 0){
    eigenvecs <- pca(x, rank = "full", type, k)$A
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
    dplyr::data_frame(rank = 0:length(pca_rank_trace_x),
                      delta_C = c(1, pca_rank_trace_x),
                      delta_residuals = c(1, pca_rank_trace_y))
}

#' PCA Rank Trace Plot
#'
#' Plot of rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @inheritParams pca
#' @param interactive logical. If \code{TRUE} prints an interactive Plotly graphic.
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[,1:34]
#' pca_rank_trace_plot(digits_features)
#'
#' @seealso \code{pca_rank_trace} \code{rank_trace_plot}
#' 
#' @export

pca_rank_trace_plot <- function(x, type = "cov", k = 0, interactive = FALSE){
    tuner <- k
    rt <- pca_rank_trace(x)
    static_plot <- ggplot2::ggplot(rt, ggplot2::aes(x = delta_C,
                                                    y = delta_residuals,
                                                    label = rank)) +
        lims(x = c(0,1), y = c(0,1)) +
        geom_line(color = "red") +
        geom_text(check_overlap = TRUE, size = 5) +
        labs(x = "dC", y = "dE") +
        ggtitle(paste("PCA Rank Trace Plot, k = ", tuner, sep = ""))
	if(interactive == TRUE){
        	plotly::ggplotly(static_plot)
    	} else {
        	static_plot
    	}
}

pca_pairwise <- function(x, pca_x, pca_y, rank = "full", type = "cov"){
    scores <- pca_scores(x, rank, type)
    score_1 <- paste("PC", pca_x, sep = "")
    score_2 <- paste("PC", pca_y, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2))
}

#' PC Pairwise Plot
#'
#' \code{pca_pairwise_plot}
#'
#' @inheritParams pca
#' @param pc_x principal component for the x-axis
#' @param pc_y principal component for the y-axis
#' @param class_labels data frame or vector of class labels
#' @param interactive logical. If \code{TRUE} prints an interactive Plotly graphic.
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[,1:34]
#' digits_class <- pendigits[,35]
#' pca_pairwise_plot(digits_features, pc_x = 1, pc_y = 3)
#'
#' @export
 
pca_pairwise_plot <- function(x, pc_x = 1, pc_y = 2, class_labels = NULL, rank = "full", type = "cov", interactive = FALSE){
    pairs <- pca_pairwise(x, pc_x, pc_y, rank, type)
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
        class <- as_data_frame(class_labels)
        pairs_tbl <- bind_cols(pairs, class)
        names(pairs_tbl) <- c("pc_x", "pc_y", "class")
        static_plot <- ggplot2::ggplot(pairs_tbl,
                                       aes(pc_x,
                                           pc_y)) +
            geom_point(aes(color = factor(class))) +
            labs(x = paste("PC", pc_x, sep = ""),
                          y = paste("PC", pc_y, sep = "")) +
            labs(title = "PC Pairwise") +
            theme(legend.title = element_blank())
    }
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }
}

pca_threewise <- function(x, pca_x, pca_y, pca_z, rank = "full", type = "cov"){
    scores <- pca_scores(x, rank, type)
    score_1 <- paste("PC", pca_x, sep = "")
    score_2 <- paste("PC", pca_y, sep = "")
    score_3 <- paste("PC", pca_z, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2), ends_with(score_3))
}


#' 3D Principal Component Plot
#'
#' code{pca_plot_3D}
#'
#' @param x data frame or matrix of predictor variables
#' @param pca_x integer number of the principal component used for the x-axis
#' @param pca_y integer number of the principal component used for the y-axis
#' @param pca_z integer number of the principal component used for the z-axis
#' @param class_labels data frame or vector of class labels
#' @param rank rank of coefficient matrix
#' @param type type of covariance matrix
#'
#' @export

pca_plot_3D <- function(x, pca_x = 1, pca_y = 2, pca_z = 3, class_labels = NULL, rank = "full", type = "cov"){
    class <- as_data_frame(class_labels)
    threes <- pca_threewise(x, pca_x, pca_y, pca_z, rank, type)
    threes_tbl <- bind_cols(threes, class)
    names(threes_tbl) <- c("x_coord", "y_coord", "z_coord", "class")
    plot_ly(threes_tbl, x = ~x_coord, y = ~y_coord, z = ~z_coord, color = ~factor(class)) #%>%
        #layout(title = "PC Scatter 3D",
         #scene = list(
         #  xaxis = list(title = names(threes)[1]), 
         #  yaxis = list(title = names(threes)[2]), 
        #   zaxis = list(title = names(threes)[3])))
}


#' Scatterplot Matrix of Principle Components
#'
#' plots all pairs of principal components
#'
#' @inheritParams pca_scores
#' @inheritParams pca_pairwise_plot
#'
#' @export

#pca_allpairs_plot <- function(x, rank, type = "cov", k = 0, class_labels = NULL){
#	class <- as.factor(class_labels[[1]])
#	pc <- pca_scores(x, rank, type, k)
#	df <- dplyr::bind_cols(pc, as_data_frame(class))
#	names(df)[dim(df)[2]] <- "class"
#	GGally::ggpairs(df, aes(color = class))
#}