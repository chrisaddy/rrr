#' Fit Reduced-Rank PCA Model
#' 
#' \code{pca} fits reduced-rank principle component analysis model.
#' 
#' @inheritParams rrr
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#'
#' @export

rrpca <- function(x, rank = "full", type = "cov", k = 0){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    s_xx <- cov(x)
    A <- eigen(s_xx)$vectors[,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    B <- t(A)
    list(means = means, A = A, B = B, C = A %*% B)
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pc_scores} returns data frame of principle component scores from reduced-rank PCA
#'
#' @inheritParams rrpca
#'
#' @export

pc_scores <- function(x, rank = "full", type = "cov", k = 0){
    pca <- rrpca(x, rank, type, k)
    scores <- t(pca$A %*% organize(x)) %>%
        as_data_frame()
    names(scores) <- paste("PC", 1:dim(scores)[2], sep = "")
    scores
}