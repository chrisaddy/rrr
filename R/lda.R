binary_matrix <- function(vec) {
    mat <- dplyr::as_data_frame(model.matrix(~ vec -1))
    mat <- dplyr::select(mat, -dim(mat)[2])
    names(mat) <- substring(names(mat), 4, 100)
    mat
}

#' Reduced-Rank Linear Discriminant Analysis
#'
#' \code{lda} produces a linear discriminant analysis as a reduced-rank regression.
#'
#'
#'
#' @export

lda <- function(x, y, rank = "full") {
### build Y matrix
    y_binary <- binary_matrix(y)
### centering matrix
    gamma <- cov(y_binary, y_binary)
    rrr(x, y_binary, gamma, rank, type = "cov")
}
