#' Fit Reduced-Rank LDA Model
#'
#' \code{lda} fits a reduced-rank linear discriminant analysis model.
#' 
#' @inheritParams rrr
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#'
#' @export

rrlda <- function(x, y, rank, type = "cov", k = 0){
    class <- y
    y_binary <- binary_matrix(y)
    cva_object <- rrcva(x, y_binary, rank, type, k)
    list(class = class, mean = cva_object$mean, G = cva_object$G, H = cva_object$H)
}

#' Linear Discriminant Scores
#'
#' \code{ld_scores} 
#'
#' @inheritParams rrlda
#'
#' @export

ld_scores <- function(x, y, rank, type = "cov", k = 0){
    class <- y
    y_binary <- binary_matrix(y)
    lda_object <- rrlda(x, y, rank, type, k)
    xi <- as_data_frame(t(lda_object$G %*% organize(x)))
    omega <- as_data_frame(t(lda_object$H %*% organize(y_binary)))
    list(ldf = dplyr::bind_cols(xi, class), class_mean = unique(dplyr::bind_cols(omega, class)))
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
