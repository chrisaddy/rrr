cv_errors <- function(x, y, gamma_mat, k = 5) {
				folds <- cvFolds(dim(y)[2], k)
				l <- c()
				for(j in 1:dim(y)[1]){
					l[j] <- rank_error_rate(x, y, gamma_mat, j, folds)
				}
				l
}