#' Reduced-Rank Linear Discriminant Analysis
#'
#' \code{lda} produces a linear discriminant analysis as a reduced-rank regression.
#' 
#' @inheritParams rrr
#'
#' @export

rrlda <- function(x, y, k = 0, ld_num = 2) {
### build Y matrix
    y_binary <- binary_matrix(y)
### centering matrix
    n <- dim(y_binary)[1]
    centering_matrix <- diag(1, n) - matrix(1, n, n) / n
    y_center <- t(y_binary) %*% centering_matrix
    n_class <- as.vector(colSums(y_binary))
    x_center <- organize(x)
    s_yy_inv <- y_center %*% t(y_center)###diag(1 / n_class) + (1 / n_class[length(n_class)]) * matrix(1, length(n_class), length(n_class))
    s_xy <- x_center %*% t(y_center)
    s_yx <- t(s_xy)
    s_xx_inv_sqrt <- solve(sqrt_matrix(x_center %*% t(x_center) + k * diag(1, dim(x_center)[1])))
    R_star <- s_xx_inv_sqrt %*%
        s_xy %*% s_yy_inv %*%
        s_yx %*%
        s_xx_inv_sqrt
    R_star
    V_star <- eigen(R_star)$vectors[,1:ld_num]
    gamma <- s_xx_inv_sqrt %*% Re(V_star)
    xi <- t(gamma) %*% x_center
    omega <- as_data_frame(t(gamma) %*% s_xy %*% s_yy_inv %*% y_center)
    xi_df <- as_data_frame(t(xi))
    omega_df <- as_data_frame(t(omega))
    list(gamma = gamma,
         xi = bind_cols(xi_df, y),
         class_means = bind_cols(omega_df, y))
}