lda_organize <- function(features, classes){
    names(classes) <- "class"
    if(length(classes$class) <= 2){
	stop("too few classes for multiclass lda")
    }
    combine_df <- dplyr::bind_cols(features, classes)
    arrange_df <- dplyr::arrange(combine_df, class)
    features_ordered <- dplyr::select(arrange_df, -class)
    classes_ordered <- dplyr::select(arrange_df, class)
    list(features_ordered = features_ordered, 
	 classes_ordered = binary_matrix(classes_ordered))
}

#' Fit Reduced-Rank LDA Model
#'
#' \code{lda} fits a reduced-rank linear discriminant analysis model.
#' 
#' @inheritParams rrr
#' @param quadratic `logical`. If `TRUE` fits a reduced-rank linear discriminant model on an expanded feature space that includes all the input variables, their squares, and all of pairwise cross-products.
#' 
#' @return `list` containing: a data frame of class
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrlda <- function(x, y, rank = "full", type = "cov", k = 0, quadratic = FALSE){
    if(quadratic == TRUE){
	x <- expand_feature_space(x)
    }
    ordered <- lda_organize(x, y)
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
    n_last <- n[length(n)]
    x_organize <- organize(x_ordered)
    y_organize <- organize(y_ordered) 
    cov_y_inv <- solve(y_organize %*% t(y_organize))#diag(1 / n) + 1 / n_last * matrix(1, length(n), length(n))
    cov_x_inv_sqrt <- x_organize %*% t(x_organize) + k * diag(1, dim(x_ordered)[2])
    cov_xy <- x_organize %*% t(y_organize) 
    cov_yx <- t(cov_xy)
    r_star <- cov_x_inv_sqrt %*% cov_xy %*% cov_y_inv %*% cov_yx %*% cov_x_inv_sqrt
    eigens <- eigen(r_star)
    vecs <- eigens$vectors[,1:reduce_rank]
    gam <- cov_x_inv_sqrt %*% vecs
    h <- t(gam) %*% cov_xy %*% cov_y_inv
    G <- as_data_frame(gam)
    names(G) <- paste("LD", 1:reduce_rank, sep = "") 
    H <- as_data_frame(t(h))
    names(H) <- names(G) 
    eigen_portion <- eigens$values[1:reduce_rank] / sum(eigens$values[1:reduce_rank])
    list(G = G, H = H, eigen_portion = eigen_portion) 
}

rrlda2 <- function(x, y, rank = "full", type = "cov", k = 0, quadratic = FALSE){
	ordered <- lda_organize(x, y)
	x_ordered <- ordered$features_ordered
	y_ordered <- ordered$features_ordered
	rrcva(x_ordered, y_ordered, rank, type, k)
}	
#' Linear Discriminant Scores
#'
#' \code{ld_scores} 
#'
#' @inheritParams rrlda
#'
#' @export

ld_scores <- function(x, y, rank = "full", type = "cov", k = 0){
    class <- y
    names(class) <- "class"
    lda_object <- rrlda(x, y, rank, type, k)
    ordered <- lda_organize(x, y)
    x_organize <- organize(ordered$features_ordered)
    y_organize <- organize(ordered$classes_ordered)
    xi <- as_data_frame(t(t(lda_object$G) %*% x_organize))
    omega <- as_data_frame(t(t(lda_object$H) %*% y_organize))
    #list(ldf = dplyr::bind_cols(xi, class), 
#	 class_mean = unique(dplyr::bind_cols(omega, class)))
    list(scores = dplyr::bind_cols(xi, class), class_means = unique(dplyr::bind_cols(omega, class))) 
}

