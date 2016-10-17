pc_gamma <- function(x){
		diag(1, dim(x)[2])
}

pc <- function(x, rank){
		gamma <- pc_gamma(x)
		A_t(x, x, gamma, rank)
}

#' Reduced-Rank PCA
#' 
#' \code{pca} carries out reduced-rank principled component analysis on a
#' matrix of input data.
#' 
#' 
#' @export pca

pca <- function(x, rank){
		pc(x, rank) %>%
			as_data_frame()
		pc
}

#` Reduced Rank PCA Prediction
#`
#` \code{pca_predict} predicts 
#`
#` @export pca_predict

pca_predict <- function(x, rank) {
		gamma <- diag(1, dim(x)[2])
		coefficients <- pc(x, rank) %>% as.matrix()
		mean <- mu_t(x, x, gamma, rank)
		mean + coefficients %*% B_t(x, x, gamma, rank) %*% organize(x)
}