#` @export

train_set <- function(x, y, folds, i){
		train_x <- x[,which(folds$which != i)]
		train_y <- y[,which(folds$which != i)]
		list(x = train_x, y = train_y)
}