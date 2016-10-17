#` @export

mu_t <- function(x, y, gamma_mat, rank){
		x_organize <- organize(x)
		y_organize <- organize(y)
		mu_vars(y_organize) - C_t(x_organize, y_organize, gamma_mat, rank) %*% 
				mu_vars(x_organize)
}