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