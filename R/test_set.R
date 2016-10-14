test_set <- function(x,y,folds,i){
		test_x <- x[,which(folds$which == i)]
		test_y <- y[,which(folds$which == i)]
		list(x = test_x, y = test_y)
}