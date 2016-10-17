#' Mean Matrix
#' 
#' \code{mu_vars} produces a mean matrix.
#' 
#' @param var a data frame of 
#' 
#' @export mu_vars
mu_vars <- function(var){
				replicate(dim(var)[2],
				rowMeans(var))
}