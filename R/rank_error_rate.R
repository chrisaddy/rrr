#` @export

rank_error_rate <- function(x, y, gamma_mat, rank, folds){
			k <- folds$K
			errors <- c()
			for(i in 1:k){
				errors[i] <- error(x, y, gamma_mat, rank, folds, i)
			}
			mean(errors) / (dim(y)[2])
}
