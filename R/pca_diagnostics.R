#' PCA Goodness of Fit
#'
#' @inheritParams rrpca
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