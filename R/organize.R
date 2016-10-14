organize <- function(vars){
				vars %>%
				center_data() %>%
				as.matrix() %>%	
				t()
}
