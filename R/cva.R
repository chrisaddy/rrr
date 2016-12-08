#' Fit Reduced-Rank CVA Model
#'
#' \code{cva} fits a reduced-rank canonical variate/correlation model. This is a special case of reduced-rank regression
#' with the weight matrix set to the covariance matrix of \eqn{Y}, i.e., \eqn{\mathbf{\Gamma} = \mathbf{\Sigma}_{YY}}.
#' Canonical variate analysis creates a set of new predictor variables that are linear combinations of the orignal predictors,
#' and a set of new resonse variables that are linear combinations of the original responses, such that each pair of new predictor and new
#' response maximizes correlation between the two and all pairs of canonical variates are independent of each other.
#'
#' @inheritParams rrr
#'
#' @return list containing: matrix of means; matrix of canonical variate coefficients of predictor variables; matrix of canonical variate coefficients of response variables; vector of canonical correlations. 
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva(galaxy_x, galaxy_y, rank = 2)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @seealso \code{\link{rrr}}
#'
#' @export

cva <- function(x, y, rank = "full", type = "cov", k = 0) {
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	rrr_object <- rrr(x, y, gamma, rank, type, k)
	H <- ginv(rrr_object[["A"]])
	colnames(H) <- names(y)
    list(mean = rrr_object[["mean"]], G = rrr_object[["B"]], H = H, canonical_corr = rrr_object[["eigen_values"]])
}

#' Canonical Variate Scores
#'
#' \code{cva_scores} creates linear combinations of predictor and response variables from the coefficients of reduced-rank canonical variate analysis.
#' 
#' @inheritParams cva
#'
#' @return list containing: data frame of canonical variate scores of predictor variables; data frame of canonical variate scores of response variables; vector of canonical correlations.
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_scores(galaxy_x, galaxy_y, rank = 2)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @seealso \code{\link{cva}}
#'
#' @export

cva_scores <- function(x, y, rank = "full", type = "cov", k = 0){
	cva_object <- cva(x, y, rank, type, k)
	correlation <- cva_object[["canonical_corr"]]
	xi <- as_data_frame(t(cva_object[["G"]] %*% organize(x)))
	names(xi) <- paste("xi", 1:dim(xi)[2], sep = "")
	omega <- as_data_frame(t(cva_object[["H"]] %*% organize(y)))
	names(omega) <- paste("omega", 1:dim(omega)[2], sep = "")
	list(xi = xi, omega = omega, canonical_corr = correlation)
}