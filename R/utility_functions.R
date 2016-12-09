#' @importFrom magrittr %>%
#' @export
magrittr::'%>%'

#' @import dplyr
#' @export

#' @import ggplot2
#' @export

#' @importFrom stats cov
#' @export

#' @importFrom plotly ggplotly plot_ly
#' @export

#' @importFrom MASS ginv
#' @export
MASS::ginv

sqrt_matrix <- function(matr){
			eigens <- eigen(matr)
			vecs <- eigens[["vectors"]]
			vals <- eigens[["values"]]
			vecs %*% diag(sqrt(vals)) %*% t(vecs)
}

organize <- function(vars, type = "cov"){
	matr <- as.matrix(vars)
	if(type == "cov"){
		matr %>%
		scale(scale = FALSE) %>%
		t()
		} else if(type == "cor"){
			matr %>%
			scale() %>%
			t()
		} else {
			stop("type not recognized")
		}
}

cov_matrix <- function(x, y, type = "cov"){
	n <- dim(x)[1]
	x_org <- organize(x, type)
	y_org <- organize(y, type)
	x %*% t(y) / n
}

# Binary Indicator Matrix
#
# @param class vector or data frame of type character
#
# @return matrix of binary indicator matrix where each column represents  

binary_matrix <- function(class) {
	class <- dplyr::as_data_frame(class)
	num_class <- dim(unique(class))[1]
    class <- as.matrix(dplyr::mutate_if(class, is.factor, as.character))
    mat <- stats::model.matrix(~ class -1)[,-num_class]
    colnames(mat) <- gsub("class", "", colnames(mat))
    mat
}

square <- function(x){
	x^2
}

expand_feature_space <- function(feature_space){
	feature_space <- dplyr::as_data_frame(feature_space)
	feature_squares <- feature_space %>%
		mutate_each(funs(square))
	names(feature_squares) <- paste(names(feature_space), "_squared", sep = "")
	dplyr::bind_cols(feature_space, feature_squares)
}	

lda_organize <- function(features, classes){
    classes <- as_data_frame(classes)
    names(classes) <- "class"
    if(dim(unique(classes))[1] <= 2){
    stop("too few classes for multiclass lda")
    }
    combine_df <- dplyr::bind_cols(features, classes)
    arrange_df <- dplyr::arrange(combine_df, class)
    features_ordered <- dplyr::select(combine_df, -class)
    classes_ordered <- dplyr::select(combine_df, class)
    list(features_ordered = scale(features_ordered, scale = FALSE), 
	 classes_ordered = scale(binary_matrix(classes_ordered)), scale = FALSE)
}