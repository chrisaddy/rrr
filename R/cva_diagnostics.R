#' Residuals of Reduced-Rank Canonical Variate
#'
#' \code{cva_residuals} returns the multivariate residuals of the reduced-rank CVA regression
#'
#' @inheritParams rrcva
#'
#' @export

cva_residuals <- function(x, y, rank = "full", type = "cov"){
	cva_object <- rrcva(x, y, rank, type)
	as_data_frame(t(cva_object$H %*% organize(y) - cva_object$G %*% organize(x)))
}