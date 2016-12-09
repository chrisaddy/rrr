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


#cva_allpairs_plot <- function(x, y, rank, type = "cov", k = 0){
#	scores_object <- cva_scores(x, y, rank, type, k)
#	all_pairs <- dplyr::bind_cols(scores_object[["xi"]], scores_object[["omega"]])
#	GGally::ggpairs(all_pairs)
#}


#' Predict via Reduced-Rank Principal Component Analysis
#' 
#' \code{pca_predict} 
#' 
#' @param pca_object list object obtained from \code{pca()}
#' @param x_new data frame or matrix of new observations to predict. 
#'
#'

#pca_predict <- function(pca_object, x_new){
#	pca_object[["C"]] %*% organize(x_new)
#}

#' Reduced-Rank PCA Error
#'
#' \code{pca_error}
#'
#' @inheritParams pca_predict
#'
#' 

#pca_error <- function(pca_object, x_new){
#	x_new - pca_predict(pca_object, x_new)
#}	

#' Reduced-Rank PCA Residuals
#'
#' \code{pca_residuals}
#'
#' @inheritParams pca
#'
#'

#pca_residuals <- function(x, rank = "full", type = "cov", k = 0){
#	object <- pca(x, rank, type, k)
#	pca_error(object, x)
#}

lda <- function(x, class, rank = "full", type = "cov", k = 0, quadratic = FALSE){
    class_names <- class %>% 
        as_data_frame() %>%
        select(class = 1) %>%
        mutate(class = as.character(class)) %>%
        distinct() %>%
        arrange(class) %>%
        as.matrix() %>% 
        as.vector()
    if(quadratic == TRUE){
    x <- expand_feature_space(x)
    }
    ordered <- lda_organize(x, class)
    x_ordered <- ordered$features_ordered
    y_ordered <- ordered$classes_ordered
    full_rank <- min(dim(x_ordered)[2], dim(y_ordered)[2])
    if(rank == "full"){
    reduce_rank <- full_rank
    } else if(rank <= full_rank){
    reduce_rank <- rank
    } else {
    stop("rank out of bounds")
    }
    n <- colSums(y_ordered)
    num_classes <- dim(y_ordered)[1]
    mean_y <- n / num_classes
    prior <- c(mean_y, 1 - sum(mean_y))
    names(prior) <- class_names
    n_last <- n[length(n)]
    #x_organize <- organize(x_ordered)
    #y_organize <- organize(y_ordered) 
    cov_y_inv <- solve(cov(y_ordered) + k * diag(1, dim(y_ordered)[2]))#diag(1 / n) + 1 / n_last * matrix(1, length(n), length(n))
    cov_x_inv_sqrt <- sqrt_matrix(solve(cov(x_ordered) + k * diag(1, dim(x_ordered)[2])))
    cov_xy <- cov(x_ordered, y_ordered) 
    cov_yx <- t(cov_xy)
    r_star <- cov_x_inv_sqrt %*% cov_xy %*% cov_y_inv %*% cov_yx %*% cov_x_inv_sqrt
    eigens <- eigen(r_star)
    vecs <- Re(eigens$vectors)[,1:reduce_rank]
    gam <- Re(cov_x_inv_sqrt %*% vecs)
    h <- t(gam) %*% cov_xy %*% cov_y_inv
    G <- as_data_frame(gam)
    names(G) <- paste("LD", 1:reduce_rank, sep = "") 
    H <- as_data_frame(t(h))
    names(H) <- names(G) 
    eigen_portion <- Re(eigens$values[1:reduce_rank] / sum(eigens$values[1:reduce_rank]))
    list(prior = prior, G = G, H = H, eigen_portion = eigen_portion) 
}
