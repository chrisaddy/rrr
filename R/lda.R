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
#' data(iris)
#' iris_x <- iris[,1:4]
#' iris_y <- iris[5]
#' lda(iris_x, iris_y)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

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

lda <- function(x, class, rank = "full", type = "cov", k = 0, quadratic = FALSE){
	ordered <- lda_organize(x, class)
	x_organize <- ordered[["features_ordered"]]
	y_organize <- ordered[["classes_ordered"]]
    full_rank <- min(dim(x_organize)[2], dim(y_organize)[2])
    if(rank == "full"){
    reduce_rank <- full_rank
    } else if(rank <= full_rank){
    reduce_rank <- rank
    } else {
    stop("rank out of bounds")
    }
	cva_object <- cva(x_organize, y_organize, rank, type, k)
    G <- as_data_frame(t(cva_object$G))
    names(G) <- paste("LD", 1:reduce_rank, sep = "")
    H <- as_data_frame(cva_object$H)
    names(H) <- names(G)
    class_names <- class %>% 
        as_data_frame() %>%
        select(class = 1) %>%
        mutate(class = as.character(class)) %>%
        distinct() %>%
        arrange(class) %>%
        as.matrix() %>% 
        as.vector()
    n <- colSums(y_organize)
    num_classes <- dim(y_organize)[1]
    mean_y <- n / num_classes
    prior <- c(mean_y, 1 - sum(mean_y))
    names(prior) <- class_names
    list(prior = prior, G = G, H = H)
}

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

lda_scores <- function(x, class, rank = "full", type = "cov", k = 0, quadratic = FALSE){
    if(quadratic == TRUE){
        x <- expand_feature_space(x)
    }
    class_label <- as_data_frame(factor(as_data_frame(class)[[1]]))
    names(class_label) <- "class"
    lda_object <- lda(x, class_label, rank, type, k)
    ordered <- lda_organize(x, class_label)
    x_organize <- ordered[["features_ordered"]]
    y_organize <- ordered[["classes_ordered"]]
    xi <- as_data_frame(t(t(as.matrix(lda_object[["G"]])) %*% t(x_organize)))
    names(xi) <- paste("LD", 1:dim(xi)[2], sep = "")
    omega <- as_data_frame(t(as.matrix(lda_object[["H"]]) %*% t(y_organize)))
    names(omega) <- names(xi)
    list(scores = dplyr::bind_cols(xi, class_label), class_means = dplyr::distinct(dplyr::bind_cols(omega, class_label))) 
}