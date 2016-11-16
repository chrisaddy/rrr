lda_organize <- function(features, classes){
    names(classes) <- "class"
    combine_df <- dplyr::bind_cols(features, classes)
    arrange_df <- dplyr::arrange(combine_df, class)
    features_ordered <- dplyr::select(arrange_df, -class)
    classes_ordered <- dplyr::select(arrange_df, class)
    list(features_ordered = features_ordered, classes_ordered = binary_matrix(classes_ordered))
}

#' Fit Reduced-Rank LDA Model
#'
#' \code{lda} fits a reduced-rank linear discriminant analysis model.
#' 
#' @inheritParams rrr
#' 
#' @return `list` containing: a data frame of class
#'
#' @examples
#' rrlda()
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

rrlda <- function(x, y, rank = "full", type = "cov", k = 0){
    ordered <- lda_organize(x, y)
    x_ordered <- ordered$features_ordered
    y_ordered <- ordered$classes_ordered
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
    lda_object <- rrlda(x, y, rank, type, k)
    ordered <- lda_organize(x, y)
    x_organize <- organize(ordered$features_ordered)
    #y_organize <- organize(ordered$classes_ordered)

    #xi <- as_data_frame(t(lda_object$G %*% x_organize)
   # omega <- as_data_frame(t(lda_object$H %*% y_organize)
    #list(ldf = dplyr::bind_cols(xi, class), class_mean = unique(dplyr::bind_cols(omega, class)))
}

### Deprecated version of rrlda
#rrlda <- function(x, y, rank = 2, type = "cov", k = 0) {
### build Y matrix
#    y_binary <- binary_matrix(y)
### centering matrix
#    n <- dim(y_binary)[1]
#    centering_matrix <- diag(1, n) - matrix(1, n, n) / n
#    y_center <- t(y_binary) %*% centering_matrix
#    n_class <- as.vector(colSums(y_binary))
#    x_center <- organize(x)
#    s_yy_inv <- y_center %*% t(y_center)###diag(1 / n_class) + (1 / n_class[length(n_class)]) * matrix(1, length(n_class), length(n_class))
#    s_xy <- x_center %*% t(y_center)
#    s_yx <- t(s_xy)
#    s_xx_inv_sqrt <- solve(sqrt_matrix(x_center %*% t(x_center) + k * diag(1, dim(x_center)[1])))
#    R_star <- s_xx_inv_sqrt %*%
#        s_xy %*% s_yy_inv %*%
#        s_yx %*%
#        s_xx_inv_sqrt
#    R_star
#    V_star <- eigen(R_star)$vectors[,1:rank]
#    gamma <- s_xx_inv_sqrt %*% Re(V_star)
#    xi <- t(gamma) %*% x_center
#    omega <- as_data_frame(t(gamma) %*% s_xy %*% s_yy_inv %*% y_center)
#    xi_df <- as_data_frame(t(xi))
#    omega_df <- as_data_frame(t(omega))
#    list(gamma = gamma,
#         xi = bind_cols(xi_df, y),
#         class_means = bind_cols(omega_df, y))
#}

#' @export
