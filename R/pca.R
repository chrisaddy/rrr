#' Fit Reduced-Rank Principal Component Model
#' 
#' \code{pca} fits reduced-rank principal component analysis model.
#' 
#' @inheritParams rrr
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' pca(digits_features, rank = 3)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

pca <- function(x, rank = "full", type = "cov", k = 0){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    s_xx <- cov(x) + k * diag(1, dim(x)[2])
    A <- eigen(s_xx)[["vectors"]][,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    list(means = as_data_frame(means), C = as_data_frame(A %*% t(A)), PC = as_data_frame(A))
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pca_scores} returns data frame of principal component scores from reduced-rank PCA
#'
#' @inheritParams pca
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' pca_scores(digits_features, rank = 3)
#'
#' @return data frame with \code{rank} number of columns, each of which represent the principal component scores of the observations.
#'
#' @export

pca_scores <- function(x, rank = "full", type = "cov", k = 0){
    pca <- pca(x, rank, type, k)
    scores <- t(t(as.matrix(pca[["PC"]])) %*% organize(x)) %>%
        as_data_frame()
    names(scores) <- paste("PC", 1:dim(scores)[2], sep = "")
    scores
}

#' Predict via Reduced-Rank Principal Component Analysis
#' 
#' \code{pca_predict} 
#' 
#' @param pca_object `list` object obtained from `pca()`
#' @param x_new data frame or matrix of new observations to predict. 
#'
#'

pca_predict <- function(pca_object, x_new){
	pca_object[["C"]] %*% organize(x_new)
}

#' Reduced-Rank PCA Error
#'
#' \code{pca_error}
#'
#' @inheritParams pca_predict
#'
#' 

pca_error <- function(pca_object, x_new){
	x_new - pca_predict(pca_object, x_new)
}	

#' Reduced-Rank PCA Residuals
#'
#' \code{pca_residuals}
#'
#' @inheritParams pca
#'
#'

pca_residuals <- function(x, rank = "full", type = "cov", k = 0){
	object <- pca(x, rank, type, k)
	pca_error(object, x)
}
