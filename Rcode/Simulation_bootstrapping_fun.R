
simulation_bootstrapping_fun <- function(n, mv_cond =NULL, tissue_cond = NULL, alpha = 0.05, 
	Iter = 1000, hist_data = data, beta = beta, var = var_select){

	## n -- sample size
	## alpha -- significance level
	## hist_data is training data, which is used to generate the empirical distributions
	## beta -- the estimated coefficients
	## var -- pass the selected variables using Elastic Net

	## Distribution of explantory variables are fixed
	if(!is.null(mv_cond))
		hist_data <- subset(hist_data, eval(parse(text=mv_cond)))
	if(!is.null(tissue_cond))
		hist_data<- subset(hist_data, eval(parse(text=tissue_cond)))

	if(nrow(hist_data) < 10) return("data size is too small")
	p_value <- rep(NA, Iter)
	## Iteration
	for(i in 1:Iter){

		#Simulate Patient profiles using bootstrapping, i.e. resampling
 		index <- sample(1:nrow(hist_data), n, replace = TRUE)
		simu_data <- hist_data[index,]

		## Generate response variable -- Nosebleeds
		X <- model.matrix(~ treatment *(mucus.viscosity + tissue.use + previous.year+ country),data=simu_data)
		X <- X[,var]
		lambda <- exp(X %*% beta)
		nosebleeds <- rpois(n, lambda)

		## Test treatment effect
		#p_value[i] <- summary(glm(simu_data$nosebleeds ~ simu_data$treatment, family = poisson))$coefficients[2,4]
		p_value[i] <- summary(glm( nosebleeds ~ simu_data$treatment, family = poisson))$coefficients[2,4]
	}

	return(mean(p_value < alpha))
}