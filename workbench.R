#' @importFrom magrittr %>%
#' @export
magrittr::"%>%"

#' @importFrom dplyr select as_data_frame
#' @export

#' @import ggplot2
#' @export

#' @importFrom plotly plot_ly ggplotly
#' @export

#' @importFrom stats cov model.matrix
#' @export

cv_errors <- function(x, y, gamma_mat, k = 5) {
				folds <- cvFolds(dim(y)[2], k)
				l <- c()
				for(j in 1:dim(y)[1]){
					l[j] <- rank_error_rate(x, y, gamma_mat, j, folds)
				}
				l
}


error <- function(x, y, gamma_mat, rank, folds, i){
		predict <- train_test(x, y, gamma_mat, rank, folds, i)
		e <- test_set(x, y, folds, i)$y - predict
		sum(e^2) / (dim(e)[1] * dim(e)[2])
}

fold <- function(data, k){
		full_set <- data[sample(nrow(data)),]
		folds <- cut(seq(1, nrow(full_set)),
			     breaks = k,
			     labels = FALSE)
		for(i in 1:k){
			test_i <- which(folds == i, arr.ind = TRUE)
			test <- data[test_i, ]
			train <- data[-test_i, ]
			train	
		}
}


#' Ridge Regression
#'
#' \code{ridge} fits a generalized ridge regression model by minimizing sum of squared error subject to an elliptical restriction on model parameters.
#'
#' @param x data frame of input variables
#' @param y vector of response variables
#' @param k small constant to augment the diagonal of the covariance matrix \eqn{\Sigma_{XX}}
#' @radius radius of 

ridge <- function(x, y, k = 0){
			x_mat <- as.matrix(x)
			y_mat <- as.matrix(y)
			solve((t(x_mat) %*% x + k * diag(1, dim(x_mat)[2]))) %*% t(x_mat) %*% y_mat
}

test_set <- function(x,y,folds,i){
		test_x <- x[,which(folds$which == i)]
		test_y <- y[,which(folds$which == i)]
		list(x = test_x, y = test_y)
}

train_set <- function(x, y, folds, i){
		train_x <- x[,which(folds$which != i)]
		train_y <- y[,which(folds$which != i)]
		list(x = train_x, y = train_y)
}

train_test <- function(x, y, gamma_mat, rank, folds, i){
		train <- train_set(x, y, folds, i)
		train_x <- train$x
		train_y <- train$y
		test <- test_set(x, y, folds, i)
		test_x <- test$x
		mu <- mu_t(train_x, train_y, gamma_mat, rank)
		C <- C_t(train_x, train_y, gamma_mat, rank)
		mu[,1:dim(test_x)[2]] + C %*% test_x
} 

#` Penalized Regression
#`
#` \code{penalized} is used to fit least-squares models subject to constraints on the parameters.
#`
#` @param 
#`
#`



rrr_predict <- function(x, y, x_new, gamma_matrix, rank = "full", type = "cov"){
	regress <- rrr(x, y, gamma_matrix, rank, type)
	regress$mean + regress$C %*% organize(x_new)
}

rrr_residuals <- function(x, y, gamma_matrix, rank){
	resid <- rrr_predict(x, y, gamma_matrix, rank) - organize(y)
	mse <- sum(resid^2) / (dim(resid)[1] * dim(resid)[2])
	list(MSE = mse, Residuals = resid)
}

#' Reduced-rank Principal Component Scores
#'
#' \code{pc_scores} returns data frame of principle component scores from reduced-rank PCA
#'
#' @param x a data frame or matrix of predictor variables.
#' @param rank rank of coefficient matrix.
#' @param type type of covariance matrix.
#'
#' @export

pc_scores <- function(x, rank = "full", type = "cov"){
    pca_object <- pca(x, rank, type)
    t(pca_object$B %*% organize(x))
}


#' Reduced-Rank PCA Prediction
#'
#' \code{pca_predict} predicts
#'
#' @param pca_object a reduced-rank PCA object from \code{pca()}.
#' @param new_x a data frame or matrix inputs to be predicted.
#'
#' @references Izenman, A. J. (2008) Modern Multivariate Statistical Techniques. Springer.
#' @export pca_predict

pca_predict <- function(pca_object, new_x){
    pca_object$C %*% new_x
}

#' PCA Goodness of Fit
#'
#' @param x data frame or matrix of observations used in a principal components analysis
#' 
#' @export

pca_gof <- function(x) {
			x_organize <- organize(x)
			eigens <- eigen(cov_matrix(x_organize, x_organize))$values
			total_var <- sum(eigens)
			gof <- c()
			for(i in 1:length(eigens)){
				gof[i] <- sum(eigens[(i + 1):length(eigens)]) / total_var
			}
			gof
}

pca_rank_trace <- function(x){
    eigenvecs <- pca(x)$A
    eigens <- eigen(cov(x, x))
    full_rank <- dim(eigens$vectors)[2]
    delta_C <- function(t){
        sqrt(1 - t / full_rank)
    }
    delta_residuals <- function(t){
        sqrt((sum(eigens$values[(t+1):full_rank])^2) / sum((eigens$values))^2)
    }
    pca_rank_trace_x <- c()
    pca_rank_trace_y <- c()
    for(i in 1:full_rank){
        pca_rank_trace_x[i] <- delta_C(i)
        pca_rank_trace_y[i] <- delta_residuals(i)
    }
    pca_rank_trace_y[length(pca_rank_trace_y)] <- 0
    data_frame(rank = 0:length(pca_rank_trace_x),
                      delta_C = c(1, pca_rank_trace_x),
                      delta_residuals = c(1, pca_rank_trace_y))
}

#' PCA Rank Trace Plot
#'
#' Plot of rank trace to determine number of principle components to use in reduced-rank PCA.
#'
#' @param x data frame or matrix of input variables
#' @param interactive logical. If \code{TRUE} prints an interactive Plotly graphic.
#'
#' @export

pca_rank_trace_plot <- function(x, interactive = FALSE){
    rt <- pca_rank_trace(x)
    static_plot <- ggplot(rt) +
        aes(x = delta_C,
                     y = delta_residuals,
                     label = rank) +
        lims(x = c(0,1), y = c(0,1)) +
        geom_line(color = "red") +
        geom_text(check_overlap = TRUE, size = 5) +
        labs(x = "dC", y = "dE") +
        ggtitle(expression(paste("Rank Trace Plot for ",
                                          Theta^(0),
                                          " to ",
                                          Theta^(s) )))
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }
}

pc_scores <- function(x, rank = "full", type = "cov"){
    as_data_frame(t(pca(x, rank, type)$B %*% organize(x)))
}


pc_pairwise <- function(x, pc_1, pc_2, rank = "full", type = "cov"){
    scores <- pc_scores(x, rank, type)
    score_1 <- paste("PC", pc_1, sep = "")
    score_2 <- paste("PC", pc_2, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2))
}

pc_threewise <- function(x, pc_x, pc_y, pc_z, rank = "full", type = "cov"){
    scores <- pc_scores(x, rank, type)
    score_1 <- paste("PC", pc_x, sep = "")
    score_2 <- paste("PC", pc_y, sep = "")
    score_3 <- paste("PC", pc_z, sep = "")
    select(scores, ends_with(score_1), ends_with(score_2), ends_with(score_3))
}


#' PC Pairwise Plot
#'
#' \code{pc_pairwise_plot}
#'
#' @param x data frame or matrix of predictor variables
#' @param pc_1 principal component for the x-axis
#' @param pc_2 principal component for the y-axis
#' @param class_labels data frame or vector of class labels
#' @param rank rank of coefficient matrix
#' @param type type of covariance matrix
#' @param interactive logical. If \code{TRUE} prints an interactive Plotly graphic.
#'
#' @export
 
pc_pairwise_plot <- function(x, pc_1 = 1, pc_2 = 2, class_labels = NULL, rank = "full", type = "cov", interactive = FALSE){
    if(is.null(class_labels)){

    }
    class <- as_data_frame(class_labels)
    names(class) <- c("class")
    pairs <- pc_pairwise(x, pc_1, pc_2, rank, type)
    pairs_tbl <- bind_cols(pairs, class)
    static_plot <- ggplot(pairs_tbl,
                    aes_q(colnames(pairs_tbl)[1],
                               colnames(pairs_tbl)[2])) +
        geom_point(aes(color = factor(class))) +
        labs(x = paste("PC", pc_1, sep = ""),
                      y = paste("PC", pc_2, sep = "")) +
        labs(title = "PC Pairwise") +
        theme(legend.title = element_blank())
    if(interactive == TRUE){
        ggplotly(static_plot)
    } else {
        static_plot
    }   
}


#' 3D Principal Component Plot
#'
#' \code{pc_plot_3D}
#'
#' @param x data frame or matrix of predictor variables
#' @param pc_x integer number of the principal component used for the x-axis
#' @param pc_y integer number of the principal component used for the y-axis
#' @param pc_z integer number of the principal component used for the z-axis
#' @param class_labels data frame or vector of class labels
#' @param rank rank of coefficient matrix
#' @param type type of covariance matrix
#'
#' @export
pc_plot_3D <- function(x, pc_x = 1, pc_y = 2, pc_z = 3, class_labels = NULL, rank = "full", type = "cov"){
    class <- as_data_frame(class_labels)
    threes <- pc_threewise(x, pc_x, pc_y, pc_z, rank, type)
    threes_tbl <- bind_cols(threes, class)
    names(threes_tbl) <- c("x_coord", "y_coord", "z_coord", "class")
    plot_ly(threes_tbl, x = ~x_coord, y = ~y_coord, z = ~z_coord, color = ~factor(class)) %>%
        layout(title = "PC Scatter 3D",
         scene = list(
           xaxis = list(title = names(threes)[1]), 
           yaxis = list(title = names(threes)[2]), 
           zaxis = list(title = names(threes)[3])))
}


##predict_rrlda <- function(rrlda_object, new_x, new_y){
##    coefs <- rrlda(rrlda_object)$
##}
##' Plot of original classes along LD axes
##'
##'
#' 

#lda_original_plot <- function(rrlda_object){
#    xi <- rrlda_object$xi
#    omega <- rrlda_object$class_means
#    ggplot(xi,
#           aes(V1, V2, color = class)) +
#        geom_point() #+
#        geom_point(aes(omega$V1,
#                       omega$V2),
#                   color = "black",
#                   shape = "M",
#                   size = 4)
#}

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

#' Canonical Variate Scores
#'
#' \code{cv_scores} calculates the canonical variate scores for \eqn{\mathbf{X}} and \eqn{\mathbf{Y}}
#'
#' @param cva_object A reduced-rank CVA object obtained by \code{cva}
#' @param x data frame or matrix of predictor variables
#' @param y data frame or matrix of response variables
#'
#' @export

cv_scores <- function(cva_object, x, y) {
    xi <- cva_object$B %*% organize(x)
    omega <- cva_object$H %*% organize(y)
}
