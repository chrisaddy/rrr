lda_organize <- function(features, classes){
    classes <- as_data_frame(classes)
    names(classes) <- "class"
    if(dim(unique(classes))[1] <= 2){
    stop("too few classes for multiclass lda")
    }
    combine_df <- dplyr::bind_cols(features, classes)
    #arrange_df <- combine_df #dplyr::arrange(combine_df, class)
    features_ordered <- dplyr::select(combine_df, -class)
    classes_ordered <- dplyr::select(combine_df, class)
    list(features_ordered = features_ordered, 
	 classes_ordered = binary_matrix(classes_ordered))
}

#' Fit Reduced-Rank LDA Model
#'
#' \code{lda} fits a reduced-rank linear discriminant analysis model.
#' 
#' @inheritParams rrr
#' @param class vector or one-column data frame of type `factor` or `character` that are the class labels of the observations.
#' @param quadratic `logical`. If `TRUE` fits a reduced-rank linear discriminant model on an expanded feature space that includes all the input variables, their squares, and all of pairwise cross-products.
#' 
#' @return `list` containing: a data frame of class
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[,-35:-36]
#' digits_class <- pendigits[,35]
#' lda(digits_features, digits_class, rank = 3, k = 0.001)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

lda <- function(x, class, rank = "full", type = "cov", k = 0, quadratic = FALSE){
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
    gam <- cov_x_inv_sqrt %*% vecs
    h <- t(gam) %*% cov_xy %*% cov_y_inv
    G <- as_data_frame(gam)
    names(G) <- paste("LD", 1:reduce_rank, sep = "") 
    H <- as_data_frame(t(h))
    names(H) <- names(G) 
    eigen_portion <- eigens$values[1:reduce_rank] / sum(eigens$values[1:reduce_rank])
    list(G = G, H = H, eigen_portion = eigen_portion) 
}


#lda <- function(x, class, rank = "full", type = "cov", k = 0, quadratic = FALSE){
#	ordered <- lda_organize(x, class)
#	x_organize <- ordered[["features_ordered"]]
#	y_organize <- ordered[["classes_ordered"]]
#	cva(x_organize, y_organize, rank, type, k)
#}

#' Linear Discriminant Scores
#'
#' \code{lda_scores} 
#'
#' @inheritParams lda
#'
#' @examples
#' data(iris)
#' lda_scores(iris[,1:4], iris[,5])
#'
#' @export

lda_scores <- function(x, class, rank = "full", type = "cov", k = 0){
    class_label <- as_data_frame(factor(as_data_frame(class)[[1]]))
    names(class_label) <- "class"
    lda_object <- lda(x, class_label, rank, type, k)
    ordered <- lda_organize(x, class_label)
    x_organize <- organize(ordered[["features_ordered"]])
    y_organize <- organize(ordered[["classes_ordered"]])
    xi <- as_data_frame(t(t(as.matrix(lda_object[["G"]])) %*% x_organize))
    names(xi) <- paste("LD", 1:dim(xi)[2], sep = "")
    omega <- as_data_frame(t(as.matrix(lda_object[["H"]]) %*% y_organize))
    names(omega) <- names(xi)
    list(scores = dplyr::bind_cols(xi, class_label), class_means = dplyr::distinct(dplyr::bind_cols(omega, class_label))) 
}