simulation_fun <- function(n, condition =NULL, alpha = 0.05, 
	Iter = 1000, hist_data = data){

	## n -- sample size
	## alpha -- significance level
	## hist_data is training data, which is used to generate the empirical distributions

	## Distribution of explantory variables are fixed
	dist_prev <- table(data$previous.year)/nrow(data)
	
	#Simulate Data
}