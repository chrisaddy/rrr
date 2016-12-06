#' PCA Goodness of Fit
#'
#' \code{pca_gof} is a measure of the goodness of fit of the \eqn{i}th principal component.
#'
#' @inheritParams pca
#' 
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

pca_gof <- function(x) {
			x_organize <- organize(x)
			eigens <- eigen(cov_matrix(x_organize, x_organize))[["values"]]
			total_var <- sum(eigens)
			gof <- c()
			for(i in 1:length(eigens)){
				gof[i] <- sum(eigens[(i + 1):length(eigens)]) / total_var
			}
			gof[length(gof)] <- 0
			gof
}