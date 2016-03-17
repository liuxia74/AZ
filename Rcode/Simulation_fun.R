
mround <- function(x,base){ 
        base*round(x/base) 
}

simulation_fun <- function(n, mv_cond =NULL, tissue_cond = NULL, alpha = 0.05, 
	Iter = 1000, hist_data = data, beta = beta, var = var_select){

	## n -- sample size
	## alpha -- significance level
	## hist_data is training data, which is used to generate the empirical distributions
	## beta -- the estimated coefficients
	## var -- pass the selected variables using Elastic Net

	## Distribution of explantory variables are fixed
	dist_prev <- prop.table(table(hist_data$previous.year))
	if(is.null(mv_cond))
		dist_mv <- prop.table(table( mround(hist_data$mucus.viscosity,.5)))
	if(!is.null(mv_cond))
		dist_mv <- prop.table(table( mround(subset(hist_data, eval(parse(text=mv_cond)))$mucus.viscosity,.5)))
	if(is.null(tissue_cond))
		dist_tissue <- prop.table(table(hist_data$tissue))
	if(!is.null(tissue_cond))
		dist_tissue <- prop.table(table(subset(hist_data, eval(parse(text=tissue_cond)))$tissue))
	## 
	dist_country <- prop.table(table(hist_data$country))

	p_value <- rep(NA, Iter)
	## Iteration
	for(i in 1:Iter){

		#Simulate Data from their empirical distribution with assumption that the variables are independent
		simu_data <- NULL
		simu_data$previous.year <- sample(as.numeric(names(dist_prev)), n, prob = as.vector(dist_prev), replace = TRUE)
		simu_data$mucus.viscosity <- sample(as.numeric(names(dist_mv)), n, prob = as.vector(dist_mv), replace = TRUE)
		simu_data$tissue.use <- sample( names(dist_tissue) , n, prob = as.vector(dist_tissue), replace = TRUE)
		simu_data$tissue.use <-  factor(simu_data$tissue.use, levels = names(dist_tissue))
		simu_data$country <- sample( names(dist_country) , n, prob = as.vector(dist_country), replace = TRUE)
		simu_data$treatment <- sample( c(0,1) , n, prob = c(.5,.5), replace = TRUE)
		simu_data <- as.data.frame(simu_data)

		## Generate response variable -- Nosebleeds
		X <- model.matrix(~ treatment *(mucus.viscosity + tissue.use + previous.year+ country),data=simu_data)
		X <- X[,var]
		lambda <- exp(X %*% beta)
		nosebleeds <- rpois(n, lambda)

		## Test treatment effect
		p_value[i] <- summary(glm(nosebleeds ~ simu_data$treatment, family = poisson))$coefficients[2,4]
	}

	return(mean(p_value < alpha))
}
