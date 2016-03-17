
#######################
## 0. Load data
#######################
data_folder <- "../Scenario/data/"

subject <- read.csv(paste(data_folder,"subject.csv", sep=''))
efficacy <- read.csv(paste(data_folder,"efficacy.csv", sep=''))
randomization <- read.csv(paste(data_folder,"randomization.csv", sep=''))


## merge the table

data <- merge(merge(subject, efficacy, by = ("subject")),
	 randomization, by = "subject")

## 
data$treatment <- ifelse(data$arm == 'ACTIVE', 1, 0)
 

## Only keep those patients with duration > 300 days
## reduce the sample size from 444 to 398
data <- subset(data, duration > 300 & !is.na(mucus.viscosity))

## impute eye colour
data$eye.colour <- as.character(data$eye.colour)
data$eye.colour[is.na(data$eye.colour)] <-"Unavailable"
data$eye.colour <- as.factor(data$eye.colour)

############################
## 1. Exploratory Analysis
############################

############################
## 1.0. Overall effectiveness
############################

overall_rslt <- glm(nosebleeds ~ treatment , data=data, family = poisson)

#############################
## 1.1. Effect of Mucus Viscosity
###############################
mucus_rslt <- glm(nosebleeds ~ treatment * mucus.viscosity, data=data, family = poisson)

## treatment effect of different level of mucus viscosity
coeff <- mucus_rslt$coefficient
cov_matrix <- summary(mucus_rslt)$cov.unscaled

mucus.viscosity = seq(0,6.5, by = .05)
mucus_effect <- coeff[2] + coeff[4] *  mucus.viscosity
mucus_stder <- sqrt(cov_matrix[2,2] + 2 * cov_matrix[2,4]*  mucus.viscosity + cov_matrix[4,4] *  mucus.viscosity ^2)

mucus_zscore <- mucus_effect / mucus_stder
## One side p_value P( z > |mucus_zscore| )
mucus_p_value <- 2 *(1- pnorm(abs(mucus_zscore), 0, 1))

## cutoff of mucus viscosity 

mucus_output <- data.frame(mucus.viscosity =  mucus.viscosity, mucus_p_value = mucus_p_value)

plot(mucus_p_value ~ mucus.viscosity, data=mucus_output, xlab = 'mucus viscosity', ylab = 'p_vlaue')
abline(0.05, 0, col = 2)

#plot(data$mucus.viscosity , mucus_effect )
#points(data$mucus.viscosity , mucus_effect + mucus_stder * 1.96 )
#points(data$mucus.viscosity , mucus_effect - mucus_stder * 1.96 )

#############################
## 1.2. Effect of Tissue usage
###############################

tissue_rslt <- glm(nosebleeds ~ treatment * tissue.use, data=data, family = poisson)

 

##############################
## 2. Statistical Model
##############################

###################################################################
## 2.1. use Elastic Net (with cross-validation) to do variable selection over the full model 
####################################################################
library(glmnet)

full_rslt <- glm(nosebleeds ~ treatment *(mucus.viscosity + tissue.use + previous.year+ country + eye.colour),data=data, family = poisson)

X <- model.matrix(~ treatment *(mucus.viscosity + tissue.use + previous.year+ country + eye.colour),data=data)

glmnet_rslt <- cv.glmnet(x=X[,-1], y = data$nosebleeds,  family = 'poisson')
glmnet_rslt5 <- cv.glmnet(x=X[,-1], y = data$nosebleeds, alpha=.5, family = 'poisson')

glmnet_rslt$cvm[which(glmnet_rslt$lambda == glmnet_rslt$lambda.1se)]
glmnet_rslt5$cvm[which(glmnet_rslt5$lambda == glmnet_rslt5$lambda.1se)]
## variables selected using lambda.1se and alpha = .5 
coef(glmnet_rslt, s = 'lambda.1se') 
coeff <- coef(glmnet_rslt5, s = 'lambda.1se') 

beta <- coeff[which(as.matrix(coeff)!=0)]
var_select <- rownames(coeff)[which(as.matrix(coeff)!=0)] 
var_select <- var_select[-1] ## remove intercept


## randomForest
rf_rslt <- randomForest(nosebleeds ~ treatment  + mucus.viscosity + tissue.use + previous.year+ country + eye.colour, data=data, importance=TRUE,
                        proximity=TRUE)

#########################################
## 3. Simulation
#########################################

source("./Simulation_fun.R")

## sample size = 800, general population
simulation_fun(n= 800, mv_cond =NULL, tissue_cond = NULL, alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

## sample size = 800, mucus viscosity > 1.2
simulation_fun(n= 800, mv_cond ="mucus.viscosity > 2", tissue_cond = NULL, alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

## sample size = 800, tissue.use == 'MEDIUM'
simulation_fun(n= 800, mv_cond = NULL, tissue_cond = "tissue.use== 'MEDIUM'", alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)


##
simulation_fun(n= 3000, mv_cond =NULL, tissue_cond = NULL, alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

## sample size = 3000, mucus viscosity > 1.2
simulation_fun(n= 3000, mv_cond ="mucus.viscosity > 1.2", tissue_cond = NULL, alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

## sample size = 3000, tissue.use == 'MEDIUM'
simulation_fun(n= 3000, mv_cond = NULL, tissue_cond = "tissue.use== 'MEDIUM'", alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

## sample size = 2000, mucus viscosity > 1.5
simulation_fun(n= 2000, mv_cond ="mucus.viscosity > 1.5", tissue_cond = NULL, alpha = 0.05, 
	Iter = 5000, hist_data = data, beta = beta, var = var_select)

