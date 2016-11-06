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

pc_pairwise <- function(x, pc_1, pc_2, rank = "full", type = "cov"){
    scores <- dplyr::as_data_frame(pca(x, rank, type)$A)
    score_1 <- paste("PC", pc_1, sep = "")
    score_2 <- paste("PC", pc_2, sep = "")
    dplyr::select(scores, ends_with(score_1), ends_with(score_2))
}

pca_pairwise_plot <- function(x, pc_1, pc_2, class_labels = NULL, rank = "full", type = "cov"){
    pairs <- pc_pairwise(x, pc_1, pc_2, rank, type)
    ggplot2::ggplot(pairs,
                    aes_string(colnames(pairs)[1],
                               colnames(pairs)[2],
                               label = class_labels)) +
        geom_point() +
        ggplot2::labs(x = paste("PC", pc_1, sep = ""),
                      y = paste("PC", pc_2, sep = ""))
}
