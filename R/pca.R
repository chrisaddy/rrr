#' Fit Reduced-Rank Principal Component Model
#' 
#' \code{pca} fits a reduced-rank principal component model. This is a special case of reduced-rank
#' regression with the weight matrix set the identity, i.e., \eqn{\mathbf{\Gamma} = \mathbf{I}}, and
#' response variables set to the predictor variables, i.e, \eqn{\mathbf{X} \equiv \mathbf{Y}}. Principal
#' componennt analysis creates a new variable, the first principal component, that is a linear combination
#' of the origianl predictor variables that maximizes variance. The second principal component is constructed
#' as a linear combination of the original predictor variables while being uncorrelated with the first principal
#' component. There are, at most, as many principal components as original predictor variables, each of which are
#' uncorrelated with each \code{i}-th principal component capture less variance than the \code{i - 1} principal
#' components.
#' 
#' @inheritParams rrr
#'
#' @return list containing: named vector of means; data frame of reduced-rank regression coefficients; data frame of principal components; named matrix of goodness-of-fit statistics for the principal components.
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' pca(digits_features, rank = 3)
#'
#' @references Izenman, A. J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

pca <- function(x, rank = "full", type = "cov", k = 0){
    if(rank == "full"){
        reduce_rank <- dim(x)[2]
    } else {
        reduce_rank <- rank
    }
    means <- colMeans(x)
    s_xx <- cov(x) + k * diag(1, dim(x)[2])
    eigens <- eigen(s_xx)
    A <- eigens[["vectors"]][,1:reduce_rank]
    colnames(A) <- paste("PC", 1:reduce_rank, sep = "")
    vals <- eigens[["values"]]
    total_var <- sum(vals)
    gof <- c()
    for(i in 1:length(vals)){
        gof[i] <- sum(vals[(i + 1):length(vals)]) / total_var
    }
    gof[length(gof)] <- 0
    gof <- gof[1:reduce_rank]
    names(gof) <- colnames(A)
    list(means = means, C = as_data_frame(A %*% t(A)), PC = as_data_frame(A),goodness_of_fit = gof)
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pca_scores} returns data frame of principal component scores from reduced-rank PCA.
#'
#' @inheritParams pca
#'
#' @return data frame.
#'
#' @examples
#' data(pendigits)
#' digits_features <- pendigits[, -35:-36]
#' pca_scores(digits_features, rank = 3)
#'
#' @return data frame with \code{rank} number of columns, each of which represent the principal component scores of the observations.
#'
#' @references Izenman, A.J. (2008) \emph{Modern Multivariate Statistical Techniques}. Springer.
#'
#' @export

pca_scores <- function(x, rank = "full", type = "cov", k = 0){
    pca <- pca(x, rank, type, k)
    scores <- t(t(as.matrix(pca[["PC"]])) %*% organize(x)) %>%
        as_data_frame()
    names(scores) <- paste("PC", 1:dim(scores)[2], sep = "")
    scores
}

#' Predict via Reduced-Rank Principal Component Analysis
#' 
#' \code{pca_predict} 
#' 
#' @param pca_object list object obtained from \code{pca()}
#' @param x_new data frame or matrix of new observations to predict. 
#'
#'

#pca_predict <- function(pca_object, x_new){
#	pca_object[["C"]] %*% organize(x_new)
#}

#' Reduced-Rank PCA Error
#'
#' \code{pca_error}
#'
#' @inheritParams pca_predict
#'
#' 

#pca_error <- function(pca_object, x_new){
#	x_new - pca_predict(pca_object, x_new)
#}	

#' Reduced-Rank PCA Residuals
#'
#' \code{pca_residuals}
#'
#' @inheritParams pca
#'
#'

#pca_residuals <- function(x, rank = "full", type = "cov", k = 0){
#	object <- pca(x, rank, type, k)
#	pca_error(object, x)
#}
