binary_matrix <- function(vec) {
    #vec <- dplyr::mutate_if(vec, is.factor, is.character)
    mat <- dplyr::as_data_frame(model.matrix(~ vec -1))
    mat <- dplyr::select(mat, -dim(mat)[2])
    names(mat) <- substring(names(mat), 4, 100)
    mat
}

#' Reduced-Rank Linear Discriminant Analysis
#'
#' \code{lda} produces a linear discriminant analysis as a reduced-rank regression.
#'
#' @param formula
#' @param data
#' @param rank
#'
#' @export

lda <- function(x, y, rank = "full") {
### build Y matrix
    y_binary <- binary_matrix(y)
### centering matrix
    n <- dim(y_binary)[2]
    centering_matrix <- diag(1, n) - matrix(1, n, n) / n
    x_center <- x %*% centering_matrix
    y_center <- y %*% centering_matrix
    gamma <- cov(y_binary, y_binary)
    beta_tau <- y_center %*% t(x_center)
}
