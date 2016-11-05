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
    A <- t(eigen(cov(x, x))$vectors[,1:reduce_rank])
    B <- t(A)
    list(means = means, A = A, B = B, C = A %*% B)
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

pca_predict <- function(pca_object, new_x) {
    pca_object$C %*% 
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
