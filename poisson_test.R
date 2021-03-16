library(keras)
library(rstan)
library(ggplot2)
library(reshape2) # for melt()
library(BioStatR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory
setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/Influenza-github/")

# Input parameters
n_iter <- 10000
n_chains <- 7

# Generate data with generic lambda value
real_lambda = 100
n_data = 1
X <- rpois(n_data, real_lambda)
r <-  rnorm(n_data,0,1)
Y <- X + r # make sure none of Y values are negative!
targets <- read.csv("data/EndOfGenAPMC_Target_vector.csv")
targets <- t(targets$x) 
targets <- cbind(targets, 50,50,50,50) #cbind(targets, 182,182,182,182) # 1 is placeholder for peak date
colnames(targets) <- c("T1.12","T2.12","T3.12","T4.12","T1.13","T2.13","T3.13","T4.13","T1.14","T2.14","T3.14","T4.14","T1.15","T2.15","T3.15","T4.15","peak1","peak2","peak3","peak4" )

n_data = 16
Y = targets
# Stan data
stan.dat = list(
  n = n_data,
  y = as.vector(Y)
)

# Run stan
s0 <- stan(file = "code/poisson_test.stan", data = stan.dat, iter = n_iter, chains = n_chains, 
           pars = c("lambda") )#, control = list(adapt_delta = 0.9,max_treedepth = 15))
s_df <- extract(s0)
lambda <- s_df$lambda
lambda_df <- as.data.frame(lambda)

# Plot estimated labmda vs real lambda
ggplot(lambda_df, aes(lambda))+
  geom_density()+
  geom_vline(xintercept = targets, color='red')
