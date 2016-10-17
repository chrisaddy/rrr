#` Rank Trace Plot for PCA
#`
#` \code{pca_rank_trace} is a wrapper for \code{rank_trace}
#`
#` export pca_rank_trace

pca_rank_trace_plot <- function(x) {
		rank_trace_plot(x, x, diag(1, dim(x)[1]))
}