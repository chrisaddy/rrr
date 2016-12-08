#' Error of Reduced-Rank Canonical Variate Analysis
#'
#' \code{cva_error} calculates the errors of a canonical-variate regression model built on a training set and applied to a test set.
#'
#' @inheritParams cva
#' @inheritParams rrr_error
#'
#' @return data frame of error.
#'
#' @examples
#' set.seed(12345)
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_train <- 
#' sample_size <- floor(0.66 * nrow(galaxy))
#' train_ind <- sample(seq_len(nrow(galaxy)), size = sample_size)
#' train <- galaxy[train_ind, ]
#' test <- galaxy[-train_ind, ]
#' train_x <- select(train, -Rmag:-chi2red)
#' train_y <- select(train, Rmag:chi2red)
#' test_x <- select(test, -Rmag:-chi2red)
#' test_y <- select(test, Rmag:chi2red)
#' cva_error(train_x, train_y, test_x, test_y, rank = 2, k = 0.001)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#' @export

cva_error <- function(x, y, x_new, y_new, rank = "full", type = "cov", k = 0){
	cva_object <- cva(x, y, rank, type, k)
	index <- data_frame(index = 1:dim(y_new)[1])
	error <- as_data_frame(t(cva_object[["H"]] %*% organize(y_new) - cva_object[["G"]] %*% organize(x_new)))
	names(error) <- paste("CV", 1:dim(error)[2], sep = "") 
	dplyr::bind_cols(index, error)
}

#' Residuals of Reduced-Rank Canonical Variate Analysis
#'
#' \code{cva_residual} returns the multivariate residuals of the reduced-rank CVA regression.
#'
#' @inheritParams cva
#'
#' @return data frame of residuals.
#'
#' @examples
#' library(dplyr)
#' data(COMBO17)
#' galaxy <- as_data_frame(COMBO17)
#' galaxy <- select(galaxy, -starts_with("e."), -Nr, -UFS:-IFD)
#' galaxy <- na.omit(galaxy)
#' galaxy_x <- select(galaxy, -Rmag:-chi2red)
#' galaxy_y <- select(galaxy, Rmag:chi2red)
#' cva_residual(galaxy_x, galaxy_y, rank = 2, k = 0.001)
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

cva_residual <- function(x, y, rank = "full", type = "cov", k = 0){
	cva_error(x, y, x, y, rank, type, k)	
}

# identity_rank <- function(rank){
# 	diag(1, rank)
# }

# lambda_rank <- function(x, y, rank, k = 0){
# 	cov_x <- cov(x) + k * diag(dim(x)[2])
# 	cov_y <- cov(y) + k * diag(dim(y)[2])
# 	cov_xy <- cov(x, y)
# 	cov_yx <- t(cov_xy)
# 	r_star <- solve(sqrt_matrix(cov_x)) %*%
# 		cov_xy %*%
# 		solve(cov_y) %*%
# 		cov_yx %*%
# 		solve(sqrt_matrix(cov_x)) 
# 	eigens <- Re(eigen(r_star)[["values"]])[1:rank]
# 	diag(eigens)
# }

#' Canonical Covariance Matrix
#'
#' \code{canonical_cov} 
#' @inheritParams cva
#' 
#' @examples
#' data(COMBO17)
#' 
#'
#'

#canonical_cov <- function(x, y, rank = "full", k = 0){
#	full_rank <- min(dim(x)[2], dim(y)[2])
#        	if(rank == "full"){
#                	reduce_rank <- full_rank
#        	} else if(rank <= full_rank){
#                	reduce_rank <- rank
#        	} else {
#                	stop("rank out of bounds")
#        	}
#	lambda <- lambda_rank(x, y, reduce_rank, k)
#	identity <- identity_rank(reduce_rank)
#	cov_mat <- rbind(cbind(lambda, lambda), cbind(lambda, identity))
#	mat_names <-c(paste("xi", 1:reduce_rank, sep = ""),
#		      paste("omega", 1:reduce_rank, sep = ""))
#	rownames(cov_mat) <- mat_names
#	colnames(cov_mat) <- mat_names
#	cov_mat 		
#}

#' Canonical Correlation Matrix
#'
#' \code{canonical_corr}
#'
#' @inheritParams cva
#'
#'

#canonical_corr <- function(x, y, rank = "full", type = "cov", k = 0){
#	full_rank <- min(dim(x)[2], dim(y)[2])
#        	if(rank == "full"){
#                	reduce_rank <- full_rank
#        	} else if(rank <= full_rank){
#                	reduce_rank <- rank
#        	} else {
#                	stop("rank out of bounds")
#        	}
#	identity <- identity_rank(reduce_rank)
#	lambda <- sqrt(lambda_rank(x, y, reduce_rank, k))
#	cov_mat <- rbind(cbind(identity, lambda), cbind(lambda, identity))
#	mat_names <-c(paste("xi", 1:reduce_rank, sep = ""),
#		      paste("omega", 1:reduce_rank, sep = ""))
#	rownames(cov_mat) <- mat_names
#	colnames(cov_mat) <- mat_names
#	cov_mat
#}