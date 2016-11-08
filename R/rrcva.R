#' Reduced-Rank Canonical Variate Analysis
#'
#' \code{cva} fits a reduced-rank canonical variate/correlation model.
#'
#' @inheritParams rrr
#'
#' @export

cva <- function(x, y, rank = "full", type = "cov", k = 0) {
	if(type == "cov"){
		y_c <- scale(y, center = TRUE, scale = FALSE)
	} else if(type == "cor"){
		y_c <- scale(y, center = TRUE, scale = TRUE)
	} else {
		stop("argument type not recognized")
	}
	gamma <- cov(y, y)
	rrr_object <- rrr(x, y, gamma, rank, type)
        rrr_object$H <- MASS::ginv(rrr_object$A)
        rrr_object
}

