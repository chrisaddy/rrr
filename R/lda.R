binary_matrix <- function(class) {
    class <- dplyr::mutate_if(class, is.factor, as.character) %>% as.matrix()
    mat <- model.matrix(~ class -1)
##    mat <- dplyr::select(mat, -dim(mat)[2])
##    names(mat) <- substring(names(mat), 4, 100)
    mat
}

#' Reduced-Rank Linear Discriminant Analysis
#'
#' \code{lda} produces a linear discriminant analysis as a reduced-rank regression.
#'
#' @param x data frame of input variables
#' @param y data frame of categorical response variables
#' @param k parameter to ensure non-singular covariance matrices
#' @param ld_num number of linear discriminants to split the data
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
    s_yy_inv <- diag(1 / n_class) + (1 / n_class[length(n_class)]) * matrix(1, length(n_class), length(n_class))
    s_xy <- x_center %*% t(y_center)
    s_yx <- t(s_xy)
    s_xx_inv_sqrt <- solve(sqrt_matrix(x_center %*% t(x_center) + k * diag(1, dim(x_center)[1])))
    R_star <- s_xx_inv_sqrt %*%
        s_xy %*% s_yy_inv %*%
        s_yx %*%
        s_xx_inv_sqrt
    R_star
    V_star <- eigen(R_star)$vectors[,1:ld_num]
    gamma <- s_xx_inv_sqrt %*% V_star
    Xi <- t(gamma) %*% x_center
    Omega <- as_data_frame(t(gamma) %*% s_xy %*% s_yy_inv %*% y_center)
    Xi_df <- as_data_frame(t(Xi))
    Omega_df <- as_data_frame(t(Omega))
    list(Xi = cbind(Xi_df, y),
         Omega = Omega_df)
}

#' Plot of original classes along LD axes
#'
#'
#' @export

lda_original_plot <- function(rrlda_object){
    Xi <- rrlda_object$Xi
    Omega <- rrlda_object$Omega
    ggplot(Xi,
           aes(V1, V2, color = class)) +
        geom_point() +
        geom_point(aes(Omega$V1,
                       Omega$V2),
                   color = "black",
                   shape = "M",
                   size = 4)
}

#' Plot of LDA Classifications along LD Axes
#'
#' \code{lda_plot}
#'
#' 


#lda_plot <- function(rrlda_object){
#}



# Reduced-Rank Quadratic Linear Discriminant Analysis
#
# \code{qlda} produces a linear discriminant analysis with quadratic bounds by introducing squares and cross-products of all variables into the feature space.
#
#

#qlda <- function(x, y){
  #  x_expand <- 
#}
