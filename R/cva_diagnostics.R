#' Error of Reduced-Rank Canonical Variate Analysis
#'
#' @inheritParams rrcva
#' @param x_new new x
#' @param y_new new y
#'
#' @export

cva_error <- function(x, y, x_new, y_new, rank = "full", type = "cov", k = 0){
	cva_object <- rrcva(x, y, rank, type, k)
	index <- data_frame(index = 1:dim(y_new)[1])
	error <- as_data_frame(t(cva_object$H %*% organize(y_new) - cva_object$G %*% organize(x_new)))
	names(error) <- paste("CV", 1:dim(error)[2], sep = "") 
	dplyr::bind_cols(index, error)
}

#' Residuals of Reduced-Rank Canonical Variate Analysis
#'
#' \code{cva_residuals} returns the multivariate residuals of the reduced-rank CVA regression
#'
#' @inheritParams rrcva
#'
#' @export

cva_residuals <- function(x, y, rank = "full", type = "cov", k = 0){
	cva_object <- rrcva(x, y, rank, type, k)
	cva_error(cva_object, x, y)	
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