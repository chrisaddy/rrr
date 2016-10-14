#' Mean Matrix
#' 
#' Mean Matrix
#' 
#' 
#' @export mu_vars
mu_vars <- function(var_matrix){
				replicate(dim(var_matrix)[2],
				rowMeans(var_matrix))
}
