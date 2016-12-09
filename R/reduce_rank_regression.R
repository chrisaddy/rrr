reduce_rank_regression <- function(x, y, gamma_matrix, rank = "full", k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	if(rank == "full"){
		reduce_rank <- full_rank
	} else if(rank <= full_rank){
		reduce_rank <- rank
	} else {
		stop("rank out of bounds")
	}
	cov_x <- cov(x) + k * diag(1, dim(x)[2])
	cov_yx <- cov(y, x)
	cov_y <- cov(y) + k * diag(1, dim(y)[2])
    	cov_xy <- t(cov_yx)
	sqrtm <- sqrt_matrix(gamma_matrix)
	weighted_matrix <- sqrtm %*%
		cov_yx %*%
		solve(cov_x) %*%
		cov_xy %*%
		sqrtm
	eigens <- eigen(weighted_matrix)
	eigen_values <- eigens[["values"]]
	V_t <- eigens[["vectors"]][,1:reduce_rank] %>%
		as.matrix(ncol = reduce_rank)
	A_t <- solve(sqrtm) %*% V_t
	rownames(A_t) <- names(y)
	B_t <- t(V_t) %*%
		sqrtm %*%
		cov_yx %*%
		solve(cov_x)
	C_t <- A_t %*% B_t
	mu_y <- colMeans(y)
	mu_x <- colMeans(x)
	mu_t <- mu_y - C_t %*% mu_x
	list(mean = mu_t, A = A_t, B = B_t, C = C_t, eigen_values = eigen_values)
}

# Predict Multivariate Responses via Reduced-Rank Regression
#
# \code{rrr_predict} predicts a matrix of responses from the coefficients of a \code{rrr} object.
#
# @inheritParams rrr
# @param rrr_object an object of type list that contains the means and coefficients of the reduced-rank regression.
# @param x_new data frame or matrix of input variables used to predict the response matrix.
#
# @return data frame of predicted values of resonse variables.
# @examples
# data(tobacco)
# set.seed(123)
# train_index <- sample(dim(tobacco)[1], 18)
# tobacco_x <- tobacco[,4:9]
# tobacco_y <- tobacco[,1:3]
# tobacco_train_x <- tobacco_x[train_index, ]
# tobacco_train_y <- tobacco_y[train_index, ]
# tobacco_test_x <- tobacco_x[-train_index, ]
# gamma <- diag(1, dim(tobacco_y)[2])
# tobacco_rrr <- rrr(tobacco_train_x, tobacco_train_y, gamma, rank = 1)
# rrr_predict(tobacco_rrr, tobacco_test_x)


rrr_predict <- function(rrr_object, x_new){
	num_obs <- dim(x_new)[1]
	coeffs <- rrr_object[["C"]]
	means <- matrix(rep(rrr_object[["mean"]], num_obs), ncol = num_obs)
	df <- t(means + coeffs %*% organize(x_new)) %>% 
		as_data_frame()
	df
}

# Reduced-Rank Regression Error
#
# \code{rrr_error} calculates the error from predicting response variables from a test set using coefficients calculated from training data.
#
# @inheritParams rrr_predict
# @param y_new data frame or matrix of observed response variables.
#
# @return data frame of prediction errors.
#
# @examples
# data(tobacco)
# set.seed(123)
# train_index <- sample(dim(tobacco)[1], 18)
# tobacco_x <- tobacco[,4:9]
# tobacco_y <- tobacco[,1:3]
# tobacco_train_x <- tobacco_x[train_index, ]
# tobacco_train_y <- tobacco_y[train_index, ]
# tobacco_test_x <- tobacco_x[-train_index, ]
# tobacco_test_y <- tobacco_y[-train_index, ]
# gamma <- diag(1, dim(tobacco_y)[2])
# tobacco_rrr <- rrr(tobacco_train_x, tobacco_train_y, gamma, rank = 1)
# rrr_error(tobacco_rrr, tobacco_test_x, tobacco_test_y)

rrr_error <- function(rrr_object, x_new, y_new){
	as_data_frame(y_new - rrr_predict(rrr_object, x_new))
}

### RRR Residuals

rrr_residual <- function(x, y, gamma_matrix, rank = "full", k = 0){
	object <- reduce_rank_regression(x, y, gamma_matrix, rank, k)
	rrr_error(object, x, y)
}

# Plot Residuals of Reduced-Rank Regression
# 
# \code{rrr_residual_plot} is a scatter plot matrix used for diagnostics of the reduced-rank regression model.
# @inheritParams rrr
#
# @return ggplot object, scatterplot matrix.

rrr_residual_plot <- function(x, y, gamma_matrix, rank = "full", k = 0, plot = TRUE){
	residuals <- rrr_residual(x, y, gamma_matrix, rank, k)
	index <- dplyr::as_data_frame(1:dim(residuals)[1])
	names(index) <- "index"
	df <- bind_cols(index, residuals)
	static_plot <- GGally::ggpairs(df)
	static_plot[1,1] <- ggplot()
	if(plot == TRUE){
		static_plot
	} else {
		residuals
	}
}

### Rank Trace Plot
###################

delta_coeff <- function(x, y, gamma_matrix, k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	coeff_full <- reduce_rank_regression(x, y, gamma_matrix, rank = full_rank, k)[["C"]]	
	delta_coeff_norm <- c()
	for(i in 1:full_rank){
		delta_coeff_norm[i] <-sqrt(sum((coeff_full - reduce_rank_regression(x, y, gamma_matrix, i, k)[["C"]])^2))
	}
	delta_c <- delta_coeff_norm / sqrt(sum(coeff_full^2))
	data_frame(dC = c(1, delta_c))
}

delta_error <- function(x, y, gamma_matrix, k = 0){
	full_rank <- min(dim(x)[2], dim(y)[2])
	resid_full <- cov(rrr_residual(x, y, gamma_matrix, full_rank, k))
	delta_resid_norm <- c()
	for(i in 1:full_rank){
		delta_resid_norm[i] <- sqrt(sum((resid_full - cov(rrr_residual(x, y, gamma_matrix, i, k)))^2))
	}
	delta_EE <- delta_resid_norm / sqrt(sum((resid_full - cov(y))^2))
	data_frame(dEE = c(1, delta_EE))
}


rrr_rank_trace <- function(x, y, gamma_matrix, k = 0, plot = TRUE, interactive = FALSE){
	dC <- delta_coeff(x, y, gamma_matrix, k)
	dEE <- delta_error(x, y, gamma_matrix, k)
	ranks <- dplyr::data_frame(ranks = 0:(dim(dC)[1] - 1))
	trace <- dplyr::bind_cols(ranks, dC, dEE)
	if(plot == FALSE){
		trace
	} else {
	rt_plot <- ggplot(trace, aes(dC, dEE, label = ranks)) +
		lims(x = c(0,1), y = c(0,1)) +
		geom_line(color = "red") +
		geom_point(size = 5) + 
		geom_text(check_overlap = TRUE, size = 4, color = "white") +
		labs(x = "dC", y = "dE") +
		ggtitle(paste("Rank Trace Plot, k =  ", k, sep =""))
		if(interactive == TRUE){
			plotly::ggplotly(rt_plot)
			} else {
				rt_plot
			}
		}
	}