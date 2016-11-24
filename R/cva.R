#' Fit Reduced-Rank CVA Model
#'
#' \code{cva} fits a reduced-rank canonical variate/correlation model.
#'
#' @inheritParams rrr
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' rrcva(galaxy_x, galaxy_y, rank = 2)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrcva <- function(x, y, rank = "full", type = "cov", k = 0) {
	gamma <- solve(cov(y) + k * diag(1, dim(y)[2]))
	rrr_object <- rrr(x, y, gamma, rank, type, k)
        rrr_object$H <- ginv(rrr_object$A)
    list(mean = rrr_object$mean, G = rrr_object$B, H = ginv(rrr_object$A), canonical_corr = rrr_object$eigen_values)
}

#' Canonical Variate Scores
#' 
#' @inheritParams rrcva
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cv_scores(galaxy_x, galaxy_y, rank = 2)
#'
#' @export

cv_scores <- function(x, y, rank = "full", type = "cov", k = 0){
	cva_object <- rrcva(x, y, rank, typ, k)
	correlation <- cva_object[["canonical_corr"]]
	xi <- as_data_frame(t(cva_object$G %*% organize(x)))
	names(xi) <- paste("xi", 1:dim(xi)[2], sep = "")
	omega <- as_data_frame(t(cva_object$H %*% organize(y)))
	names(omega) <- paste("omega", 1:dim(omega)[2], sep = "")
	list(xi = xi, omega = omega, canonical_corr = correlation)
}