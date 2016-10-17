#` Ridge Regression
#`
#` \code{ridge} fits a generalized ridge regression model by minimizing sum of squared error subject to an elliptical restriction on model parameters.
#`
#` @param x data frame of input variables
#` @param y vector of response variables
#` @param tuner small constant to augment the diagonal of 
#` @radius radius of 

ridge <- function(x, y, tuner, radius = 1){
			solve((t(x) %*% x + k * diag(1, dim(x)[2]))) %*% t(x) %*% y
}