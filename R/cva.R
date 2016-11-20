#' Fit Reduced-Rank CVA Model
#'
#' \code{cva} fits a reduced-rank canonical variate/correlation model.
#'
#' @inheritParams rrr
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrcva <- function(x, y, rank = "full", type = "cov", k = 0) {
	y_c <- organize(y, type)
	gamma <- cov(y)
	rrr_object <- rrr(x, y, gamma, rank, type)
        rrr_object$H <- ginv(rrr_object$A)
    list(mean = rrr_object$mean, G = rrr_object$B, H = ginv(rrr_object$A), eigen_values = rrr_object$eigen_values)
}


identity_rank <- function(rank){
	diag(1, rank)
}

lambda_rank <- function(x, y, rank, k = 0){
	cov_x <- cov(x) + k * diag(dim(x)[2])
	cov_y <- cov(y) + k * diag(dim(y)[2])
	cov_xy <- cov(x, y)
	cov_yx <- t(cov_xy)
	r_star <- solve(sqrt_matrix(cov_x)) %*%
		cov_xy %*%
		solve(cov_y) %*%
		cov_yx %*%
		solve(sqrt_matrix(cov_x)) 
	eigens <- Re(eigen(r_star)$values)[1:rank]
	diag(eigens)
}

#' Canonical Covariance Matrix
#'
#' \code{canonical_cov} 
#' @inheritParams rrcva
#' 
#' @examples
#' data(COMBO17)
#' 
#'
#' @export

canonical_cov <- function(x, y, rank = "full", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
        	if(rank == "full"){
                	reduce_rank <- full_rank
        	} else if(rank <= full_rank){
                	reduce_rank <- rank
        	} else {
                	stop("rank out of bounds")
        	}
	lambda <- lambda_rank(x, y, reduce_rank, k)
	identity <- identity_rank(reduce_rank)
	cov_mat <- rbind(cbind(lambda, lambda), cbind(lambda, identity))
	mat_names <-c(paste("xi", 1:reduce_rank, sep = ""),
		      paste("omega", 1:reduce_rank, sep = ""))
	rownames(cov_mat) <- mat_names
	colnames(cov_mat) <- mat_names
	cov_mat 		
}



#' Canonical Correlation Matrix
#'
#' \code{canonical_corr}
#'
#' @inheritParams rrcva
#'
#' @export

canonical_corr <- function(x, y, rank = "full", type = "cov", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
        	if(rank == "full"){
                	reduce_rank <- full_rank
        	} else if(rank <= full_rank){
                	reduce_rank <- rank
        	} else {
                	stop("rank out of bounds")
        	}
	identity <- identity_rank(reduce_rank)
	lambda <- sqrt(lambda_rank(x, y, reduce_rank, k))
	cov_mat <- rbind(cbind(identity, lambda), cbind(lambda, identity))
	mat_names <-c(paste("xi", 1:reduce_rank, sep = ""),
		      paste("omega", 1:reduce_rank, sep = ""))
	rownames(cov_mat) <- mat_names
	colnames(cov_mat) <- mat_names
	cov_mat
}
