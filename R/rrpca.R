#' Fit Reduced-Rank PCA
#' 
#' \code{pca} carries out reduced-rank principled component analysis on a
#' matrix of input data.
#' 
#' @param x data frame or matrix of input variables
#' @param rank rank of coefficient matrix. Default \code{="full"}.
#' @param type type of covariance matrix.
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
    A <- eigen(cov(x))$vectors[,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    B <- t(A)
    list(means = means, A = A, B = B, C = A %*% B)
}