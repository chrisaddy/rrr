#' Fit Reduced-Rank PCA
#' 
#' \code{pca} carries out reduced-rank principled component analysis on a
#' matrix of input data.
#' 
#' @param x data frame or matrix of input variables
#' @param rank rank of coefficient matrix. Default \code{="full"}.
#' @param type 
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#' @export pca

pca <- function(x, rank = "full", type = "cov"){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    A <- eigen(cov(x, x))$vectors[,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    B <- t(A)
    list(means = means, A = A, B = B, C = A %*% B)
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pc_scores} returns data frame of principle component scores from reduced-rank PCA
#'
#' @param x a data frame or matrix of predictor variables.
#'
#' @export

pc_scores <- function(x, rank = "full", type = "cov"){
    pca_object <- pca(x, rank, type)
    t(pca_object$B %*% organize(x))
}


#' Reduced-Rank PCA Prediction
#'
#' \code{pca_predict} predicts
#'
#' @param pca_object a reduced-rank PCA object from \code{pca()}.
#' @param new_x a data frame or matrix inputs to be predicted.
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#' @export pca_predict

pca_predict <- function(pca_object, new_x){
    pca_object$C %*% new_x
}

#' PCA Goodness of Fit
#'
#' @param x data frame or matrix of observations used in a principal components analysis
#' 
#' @export

pca_gof <- function(x) {
			x_organize <- organize(x)
			eigens <- eigen(cov_matrix(x_organize, x_organize))$values
			total_var <- sum(eigens)
			gof <- c()
			for(i in 1:length(eigens)){
				gof[i] <- sum(eigens[(i + 1):length(eigens)]) / total_var
			}
			gof
}

pca_rank_trace <- function(x){
    eigenvecs <- pca(x)$A
    eigens <- eigen(cov(x, x))
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
#' @param x data frame or matrix of input variables
#'
#' @export

pca_rank_trace_plot <- function(x, interactive = FALSE){
    rt <- pca_rank_trace(x)
    static_plot <- ggplot2::ggplot(rt) +
        ggplot2::aes(x = delta_C,
                     y = delta_residuals,
                     label = rank) +
        ggplot2::lims(x = c(0,1), y = c(0,1)) +
        ggplot2::geom_line(color = "red") +
        ggplot2::geom_text(check_overlap = TRUE, size = 5) +
        ggplot2::labs(x = "dC", y = "dE") +
        ggplot2::ggtitle(expression(paste("Rank Trace Plot for ",
                                          Theta^(0),
                                          " to ",
                                          Theta^(s) )))
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }
}

pc_scores <- function(x, rank = "full", type = "cov"){
    dplyr::as_data_frame(t(pca(x, rank, type)$B %*% organize(x)))
}


pc_pairwise <- function(x, pc_1, pc_2, rank = "full", type = "cov"){
    scores <- pc_scores(x, rank, type)
    score_1 <- paste("PC", pc_1, sep = "")
    score_2 <- paste("PC", pc_2, sep = "")
    dplyr::select(scores, ends_with(score_1), ends_with(score_2))
}

pc_threewise <- function(x, pc_x, pc_y, pc_z, rank = "full", type = "cov"){
    scores <- pc_scores(x, rank, type)
    score_1 <- paste("PC", pc_x, sep = "")
    score_2 <- paste("PC", pc_y, sep = "")
    score_3 <- paste("PC", pc_z, sep = "")
    dplyr::select(scores, ends_with(score_1), ends_with(score_2), ends_with(score_3))
}


#' PC Pairwise Plot
#'
#' \code{pc_pairwise_plot}
#'
#' @param x data frame or matrix of predictor variables
#' @param pc_1 principal component for the x-axis
#' @param pc_2 principal component for the y-axis
#' @param class_labels data frame or vector of class labels
#' @param rank rank of coefficient matrix
#' @param type type of covariance matrix
#'
#' @export
 
pc_pairwise_plot <- function(x, pc_1 = 1, pc_2 = 2, class_labels = NULL, rank = "full", type = "cov", interactive = FALSE){
    class <- dplyr::as_data_frame(class_labels)
    names(class) <- c("class")
    pairs <- pc_pairwise(x, pc_1, pc_2, rank, type)
    pairs_tbl <- dplyr::bind_cols(pairs, class)
    static_plot <- ggplot2::ggplot(pairs_tbl,
                    aes_string(colnames(pairs_tbl)[1],
                               colnames(pairs_tbl)[2])) +
        geom_point(aes(color = factor(class))) +
        ggplot2::labs(x = paste("PC", pc_1, sep = ""),
                      y = paste("PC", pc_2, sep = "")) +
        labs(title = "PC Pairwise") +
        theme(legend.title = element_blank())
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }   
}


#' 3D Principal Component Plot
#'
#' \code{pc_plot_3D}
#'
#' @param x data frame or matrix of predictor variables
#' @param pc_x integer number of the principal component used for the x-axis
#' @param pc_y integer number of the principal component used for the y-axis
#' @param pc_z integer number of the principal component used for the z-axis
#' @param class_labels data frame or vector of class labels
#' @param rank rank of coefficient matrix
#' @param type type of covariance matrix
#'
#' @export
pc_plot_3D <- function(x, pc_x = 1, pc_y = 2, pc_z = 3, class_labels = NULL, rank = "full", type = "cov"){
    class <- dplyr::as_data_frame(class_labels)
    threes <- pc_threewise(x, pc_x, pc_y, pc_z, rank, type)
    threes_tbl <- dplyr::bind_cols(threes, class)
    names(threes_tbl) <- c("x_coord", "y_coord", "z_coord", "class")
    plot_ly(threes_tbl, x = ~x_coord, y = ~y_coord, z = ~z_coord, color = ~factor(class)) %>%
        layout(title = "PC Scatter 3D",
         scene = list(
           xaxis = list(title = names(threes)[1]), 
           yaxis = list(title = names(threes)[2]), 
           zaxis = list(title = names(threes)[3])))
}
