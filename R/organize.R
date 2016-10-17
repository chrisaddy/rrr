#' Organize data for RRR
#'
#' Organize data in format suitable for reduced-rank regression
#'
#' @params vars matrix of variables to be organized
#'
#' @return centered matrix whose rows are variables and columns observations
#' @export organize

organize <- function(vars, scale = FALSE){
				vars %>%
				as.matrix() %>%
				scale(center = TRUE, scale) %>%
				t()
}
