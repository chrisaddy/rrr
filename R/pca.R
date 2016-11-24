#' Fit Reduced-Rank PCA Model
#' 
#' \code{pca} fits reduced-rank principle component analysis model.
#' 
#' @inheritParams rrr
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' rrpca(digits_features, rank = 3)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrpca <- function(x, rank = "full", type = "cov", k = 0){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    s_xx <- cov(x) + k * diag(1, dim(x)[2])
    A <- eigen(s_xx)$vectors[,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    list(means = as_data_frame(means), C = as_data_frame(A %*% t(A)), PC = as_data_frame(A))
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pc_scores} returns data frame of principle component scores from reduced-rank PCA
#'
#' @inheritParams rrpca
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' pc_scores(digits_features, rank = 3)
#'
#' @return data frame with \code{rank} number of columns, each of which represent the principal component scores of the observations.
#'
#' @export

pc_scores <- function(x, rank = "full", type = "cov", k = 0){
    pca <- rrpca(x, rank, type, k)
    scores <- t(t(as.matrix(pca$PC)) %*% organize(x)) %>%
        as_data_frame()
    names(scores) <- paste("PC", 1:dim(scores)[2], sep = "")
    scores
}

#' Predict via Reduced-Rank Principle Components Analysis
#' 
#' \code{rrpca_predict} 
#' 
#' @param rrpca_object `list` object obtained from `rrpca()`
#' @param x_new data frame or matrix of new observations to predict. 
#'
#'

pca_predict <- function(rrpca_object, x_new){
	rrpca_object$C %*% organize(x_new)
}

#' Reduced-Rank PCA Error
#'
#' \code{rrpca_error}
#'
#' @inheritParams pca_predict
#'
#' 

pca_error <- function(rrpca_object, x_new){
	x_new - pca_predict(rrpca_object, x_new)
}	

#' Reduced-Rank PCA Residuals
#'
#' \code{rrpca_residuals}
#'
#' @inheritParams rrpca
#'
#'

pca_residuals <- function(x, rank = "full", type = "cov", k = 0){
	object <- rrpca(x, rank, type, k)
	pca_error(object, x)
}
