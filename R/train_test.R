#` @export

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