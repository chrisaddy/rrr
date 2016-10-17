#` @export

error <- function(x, y, gamma_mat, rank, folds, i){
		predict <- train_test(x, y, gamma_mat, rank, folds, i)
		e <- test_set(x, y, folds, i)$y - predict
		sum(e^2) / (dim(e)[1] * dim(e)[2])
}