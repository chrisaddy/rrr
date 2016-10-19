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