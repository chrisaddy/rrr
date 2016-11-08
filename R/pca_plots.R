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
    data_frame(rank = 0:length(pca_rank_trace_x),
                      delta_C = c(1, pca_rank_trace_x),
                      delta_residuals = c(1, pca_rank_trace_y))
}

#' PCA Rank Trace Plot
#'
#' Plot of rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @param x data frame or matrix of input variables
#' @param interactive logical. If \code{TRUE} prints an interactive Plotly graphic.
#'
#' @export

pca_rank_trace_plot <- function(x, interactive = FALSE){
    rt <- pca_rank_trace(x)
    static_plot <- ggplot(rt) +
        aes(x = delta_C,
                     y = delta_residuals,
                     label = rank) +
        lims(x = c(0,1), y = c(0,1)) +
        geom_line(color = "red") +
        geom_text(check_overlap = TRUE, size = 5) +
        labs(x = "dC", y = "dE") +
        ggtitle(expression(paste("Rank Trace Plot for ",
                                          Theta^(0),
                                          " to ",
                                          Theta^(s) )))
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }
}
