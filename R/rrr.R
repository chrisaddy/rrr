weighted_mat <- function(x, y, gamma_mat){
			x_organize <- organize(x)
			y_organize <- organize(y)
			sqrt_matrix(gamma_mat) %*%
				cov_mat(y_organize, x_organize) %*%
				solve(cov_mat(x_organize, x_organize)) %*%
				cov_mat(x_organize, y_organize) %*%
				sqrt_matrix(gamma_mat)
}

V_t <- function(x, y, gamma_mat, rank){
				w <- weighted_mat(x, y, gamma_mat)
				eigen(w)$vectors[,1:rank] %>%
				as.matrix(ncol = rank)
}

A_t <- function(x, y, gamma_mat, rank){
		solve(sqrt_matrix(gamma_mat)) %*% 
			V_t(x, y, gamma_mat, rank)
}

B_t <- function(x, y, gamma_mat, rank){
		t(V_t(x, y, gamma_mat, rank)) %*%
			sqrt_matrix(gamma_mat) %*%
			cov_mat(organize(y), organize(x)) %*%
			solve(cov_mat(organize(x), organize(x)))
}

C_t <- function(x, y, gamma_mat, rank) {
		A_t(x, y, gamma_mat, rank) %*% 
			B_t(x, y, gamma_mat, rank)
}


#` Reduced-Rank Regression
#`
#` \code{rrr} fits a reduced-rank regression model.
#`
#` @param x 
#` @param y
#` @param gamma_mat
#` @param rank
#`
#` @export rrr

rrr <- function(x, y, gamma_mat, rank){
		x_c <- organize(x)
		y_c <- organize(y)
		C_t(x, y, gamma_mat, rank)
}