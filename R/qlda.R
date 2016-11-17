#' Reduced-rank Linear Discriminant Analysis with Quadratic Terms
#'
#' \code{rrqlda} fits a linear discriminant analysis after expanding the feature space to include the squares and cross products of all of the input variables.
#'
#' @inheritParams rrlda
#'
#' @export

rrqlda <- function(x, y, rank = "full", type = "cov", k = 0){
	x_expanded <- expand_feature_space(x)
	rrlda(x_expanded, y, rank, type, k)
}

#' Linear Discriminant Scores of Feature Space with Quadratic Terms
#' 
#' \code{qld_scores}
#' 
#' @inheritParams ld_scores
#' 
#' @export

qld_scores <- function(x, y, rank = "full", type = "cov", k = 0){
	x_expanded <- expand_feature_space(x)
	ld_scores(x_expanded, y, rank, type, k)
}
