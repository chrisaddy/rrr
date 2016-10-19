#' Reduced-Rank PCA
#' 
#' \code{pca} carries out reduced-rank principled component analysis on a
#' matrix of input data.
#' 
#' @param x data frame of input variables
#' @param rank of coefficient matrix. Default \code{="full"}.
#' 
#' @export pca

pca <- function(x, rank = "full", type = "cov"){
	rrr(x, x, diag(1, dim(x)[2]), rank, type)
}

#' Reduced-Rank PCA Prediction
#'
#' \code{pca_predict} predicts 
#'
#' @export pca_predict

pca_predict <- function(x, rank = "full", type = "cov") {
	rrr_predict(x, x, x, diag(1, dim(x)[2]), type = "cov")
}

#' PCA Goodness of Fit
#'
#' @export

pca_gof <- function(x) {
			x_organize <- organize(x)
			eigens <- eigen(cov_mat(x_organize, x_organize))$values
			total_var <- sum(eigens)
			gof <- c()
			for(i in 1:length(eigens)){
				gof[i] <- sum(eigens[(i + 1):length(eigens)]) / total_var
			}
			gof
}
