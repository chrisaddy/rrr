#' Fit Reduced-Rank CVA Model
#'
#' \code{cva} fits a reduced-rank canonical variate/correlation model.
#'
#' @inheritParams rrr
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#'
#' @export

rrcva <- function(x, y, rank = "full", type = "cov", k = 0) {
	y_c <- organize(y, type)
	gamma <- cov(y)
	rrr_object <- rrr(x, y, gamma, rank, type)
        rrr_object$H <- ginv(rrr_object$A)
    list(mean = rrr_object$mean, C = rrr_object$C, G = rrr_object$B, H = ginv(rrr_object$A))
}