#### load packages

require(plyr)
require(MMST)
require(ggplot2)
require(cvTools)
require(matrixcalc)

### format data into x,y centered matrices

center_data <- function(data){
					scale(data, scale = FALSE)
}

organize <- function(vars){
				vars %>%
				center_data() %>%
				as.matrix() %>%	
				t()
}

## estimates for variables means

mu_vars <- function(var_matrix){
				replicate(dim(var_matrix)[2],
				rowMeans(var_matrix))
}


## estimates for the covariance matrices

cov_mat <- function(var1, var2){
		var1 %*% t(var2) / (dim(var1)[2] - 1)
}

sqrt_matrix <- function(matr){
			e <- eigen(matr)
			vecs <- e$vectors
			vals <- e$values
			vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs)
}

weighted_mat <- function(var_x, var_y, gamma_mat){ 
			sig_XX <- cov_mat(var_x, var_x)
			sig_XY <- cov_mat(var_x, var_y)
			sig_YX <- cov_mat(var_y, var_x) 
			sqrt_matrix(gamma_mat) %*%
				sig_YX %*%
				solve(sig_XX) %*%
				sig_XY %*%
				sqrt_matrix(gamma_mat)
}

V_t <- function(var_x, var_y, gamma_mat, rank){
       		eigen(weighted_mat(var_x,
       						   var_y,
       						   gamma_mat))$vectors[,1:rank]
}

A_t <- function(var_x, var_y, gamma_mat, rank){
		solve(sqrt_matrix(gamma_mat)) %*% V_t(var_x,
											  var_y,
											  gamma_mat,
											  rank)
}

B_t <- function(var_x, var_y, gamma_mat, rank){
		t(V_t(var_x, var_y, gamma_mat, rank)) %*%
			sqrt_matrix(gamma_mat) %*%
			cov_mat(var_y, var_x) %*%
			solve(cov_mat(var_x, var_x))
}

C_t <- function(var_x, var_y, gamma_mat, rank) {
		A_t(var_x,
			var_y,
			gamma_mat,
			rank) %*% B_t(var_x,
					   	  var_y,
						  gamma_mat,
						  rank)
}

theta_full <- function(var_x, var_y, gamma_mat){
		s <- dim(var_y)[1]
		C_t(var_x, var_y, gamma_mat, s)
}

mu_t <- function(var_x, var_y, gamma_mat, rank){
		mu_vars(var_y) - 
			C_t(var_x, var_y, gamma_mat, rank) %*% 
				mu_vars(var_x)
}	       

sigma_ee_t <- function(var_x, var_y, gamma_mat, rank){
		e <- var_y - C_t(var_x, var_y, gamma_mat, rank) %*% var_x
		e %*% t(e) / (dim(var_x)[2] - dim(var_x)[1])
}

sigma_ee_full <- function(var_x, var_y, gamma_mat){
       			s <- dim(var_y)[1]
			sigma_ee_t(var_x, var_y, gamma_mat, s)
}

delta_C <- function(var_x, var_y, gamma_mat, rank){
		theta <- theta_full(var_x, var_y, gamma_mat)
		C <- C_t(var_x, var_y, gamma_mat, rank)
		norm(theta - C, type = "F") / norm(theta, type = "F")
}	

delta_EE <- function(var_x, var_y, gamma_mat, rank){
		full <- sigma_ee_full(var_x, var_y, gamma_mat)	
		t <- sigma_ee_t(var_x, var_y, gamma_mat, rank)
		sig_YY <- cov_mat(var_y, var_y)
		norm(full - t, type = "F") / norm(full - sig_YY, type = "F")
}

rank_trace_x <- function(var_x, var_y, gamma_mat){
			rt <- c()
			fr <- min(dim(var_x)[1], dim(var_y)[1])
			for(i in 1:fr){
				rt[i] <- delta_C(var_x, var_y, gamma_mat, i)
			}
			c(1,rt)
}

rank_trace_y <- function(var_x, var_y, gamma_mat){
			rt <- c()
			fr <- min(dim(var_x)[1], dim(var_y)[1])
			for(i in 1:fr){
				rt[i] <- delta_EE(var_x, var_y, gamma_mat, i)
			}
			c(1,rt)
}

rank_trace_plot <- function(var_x, var_y, gamma_mat){
			rx <- rank_trace_x(var_x, var_y, gamma_mat)
			ry <- rank_trace_y(var_x, var_y, gamma_mat)
			ggplot() +
			aes(x = rx, y = ry) +
			lims(x = c(0,1), y = c(0,1)) +
			geom_point() +
			geom_line(color = "red") +
			labs(x = "dC", y = "dE") +
			ggtitle(expression(paste("Rank Trace Plot for ",
						 Theta^(0),
						 " to ",
						 Theta^(s) )))
}
				

### Cross-validation

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

### Cross Validation
####################

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

error <- function(x, y, gamma_mat, rank, folds, i){
		predict <- train_test(x, y, gamma_mat, rank, folds, i)
		e <- test_set(x, y, folds, i)$y - predict
		frobenius.norm(e) / rank
}

rank_error_rate <- function(x, y, gamma_mat, rank, folds){
			k <- folds$K
			errors <- c()
			for(i in 1:k){
				errors[i] <- error(x, y, gamma_mat, rank, folds, i)
			}
			mean(errors) / (dim(y)[2])
}

cv_errors <- function(x, y, gamma_mat, k = 5) {
				folds <- cvFolds(dim(y)[2], k)
				l <- c()
				for(j in 1:dim(y)[1]){
					l[j] <- rank_error_rate(x, y, gamma_mat, j, folds)
				}
				l
}
