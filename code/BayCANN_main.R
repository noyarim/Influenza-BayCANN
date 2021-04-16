# Code last updated on: 2021. 03.19
# Original code: main_multiceptron_Influenza4_2BayCANN_unscaled_randSE(usethis).R
# Change made: 
# (1) 

#install.packages("keras")
#install_keras()
#install_keras(tensorflow = "gpu")
#install.packages("curl")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("BioStatR")
#install.packages("rstan")
library(keras)
library(rstan)
library(ggplot2)
library(reshape2) # for melt()
library(BioStatR)
library(MASS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory
setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/Influenza-github/")

# Functions
source("code/BayCANN_functions.R")

# Set transformation of outcomes for ANN
log = 1 # log transformation 

# =================================================
# PSA Samples
# training
param <- read.csv("output/PSA_output/x_train.csv")
outcome <- read.csv("output/PSA_output/y_train.csv")
param <- param[,-c(1)]
outcome <- outcome[,-c(1)]
# testing
param_t <- read.csv("output/PSA_output/x_test.csv")
outcome_t <- read.csv("output/PSA_output/y_test.csv")
param_t <- param_t[,-c(1)]
outcome_t <- outcome_t[,-c(1)]
# scale the PSA inputs and outputs for ANN 
prepared_data <- prepare_data(as.matrix(param), as.matrix(outcome), as.matrix(param_t), as.matrix(outcome_t), log=1)
# save to the environment
list2env(prepared_data, envir = .GlobalEnv)

# Target data
targets <- read.csv("data/EndOfGenAPMC_Target_vector.csv")
targets <- t(targets$x) 
targets <- cbind(targets, 50,50,50,50) #cbind(targets, 182,182,182,182) # 1 is placeholder for peak date
colnames(targets) <- c("T1.12","T2.12","T3.12","T4.12","T1.13","T2.13","T3.13","T4.13","T1.14","T2.14","T3.14","T4.14","T1.15","T2.15","T3.15","T4.15","peak1","peak2","peak3","peak4" )
true_targets_mean <- targets
#true_targets_mean <- log(targets) #c(log(targets[1,1:ncol(targets)-1]),targets[1,ncol(targets)])
#true_targets_se <- true_targets_mean*0.1 #TODO: standard error is arbitrary.

# =====================================
# Scale the targets and their SE
y_targets <- true_targets_mean
#y_targets <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1
#y_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
y_targets <- as.matrix(y_targets)
#y_targets_se <- 0.1 * y_targets
#y_targets <- exp(y_targets)
#y_targets_se <- matrix(0.07,ncol=ncol(y_targets),nrow=1) #t(as.matrix(y_targets_se)) #matrix(0.1,ncol=16,nrow=1) 
#y_targets_se[17:20] = 0.07 # smaller se for peak dates
#y_targets_se <- 0.1 * y_targets

# # ======================================
# # Verify the truth is within the bounds of the x's!
# par(mfrow = c(3,3), mar=c(1,1,1,1))
# for (i in 1:9){
#   hist(Xunscaled[,i], main = x_names[i])
#   abline(v=x_true_unscaled[i],col="red")
# }

# ============== TensorFlow Keras ANN Section ========================
# Input parameters
n_iter <- 10000
n_hidden_nodes <- 200
n_hidden_layers <- 5 - 1 # technically it is easier if we just subtract 1 frmo the numbrer of layers
n_epochs <- 10000
verbose <- 0
n_chains <- 4

n_outputs=1


model <- keras_model_sequential() 
model %>% 
  layer_dense(units = n_hidden_nodes, activation = 'tanh', input_shape = n_inputs)%>%#, kernel_regularizer = regularizer_l2(0.0000001)) %>% 
  layer_dense(units = n_hidden_nodes, activation = 'tanh')%>%#, kernel_regularizer = regularizer_l2(0.0000001)) %>%
  layer_dense(units = n_hidden_nodes, activation = 'tanh')%>%#, kernel_regularizer = regularizer_l2(0.0000001)) %>%
  layer_dense(units = n_hidden_nodes, activation = 'tanh')%>%#, kernel_regularizer = regularizer_l2(0.0000001)) %>%
  layer_dense(units = n_hidden_nodes, activation = 'tanh')%>%#, kernel_regularizer = regularizer_l2(0.0000001)) %>%
  layer_dense(units = n_outputs)
summary(model)
#  opt = optimizer_adam(lr=0.0001) # default is 0.001
model %>% compile(
  loss = 'mean_squared_error',
  optimizer = 'adam'
)
keras.time <- proc.time()
history <- model %>% fit(
  xtrain_scaled, ytrain_scaled[,2], 
  epochs = n_epochs, batch_size = n_train, 
  validation_split = 0,
  validation_data= list(xtest_scaled,ytest_scaled[,2]),
  verbose = verbose
  #callbacks=list(es)
)
weights <- get_weights(model) #get ANN weights
pred <- model %>% predict(xtest_scaled)
proc.time() - keras.time #keras ann fitting time


# plot ANN convergence
#  png(filename=paste0('output/ann_convergence_influenza_4lyrs_',paste0(n_hidden_nodes,'.png')))
#  plot(history)
#  dev.off()

loss <- history$metrics$loss
write.csv(loss, paste0('output/ann_convergence(loss)_4lyrs(2LHS,Raw)_',paste0(n_hidden_nodes,'.csv')))
val_loss <- history$metrics$val_loss
write.csv(val_loss, paste0('output/ann_convergence(valloss)_4lyrs(2LHS,Raw)_',paste0(n_hidden_nodes,'.csv')))

# plot ANN validation
model %>% evaluate(xtest_scaled, ytest_scaled[,1])
png(filename=paste0('output/ann_validation_vs_observed_influenza_4lyrs(2LHS,Raw)_',paste0(n_hidden_nodes,'.png')))
par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
for (o in 1:n_outputs){
  plot(ytest_scaled[,o], pred[,o], xlim=c(-1,1),ylim=c(-1,1)) #main = y_names[o])
}
dev.off()
# plot ANN fit (on scale)
png(filename=paste0('output/ann_fit_vs_observed_influenza_4lyrs(2LHS,Raw)_',paste0(n_hidden_nodes,'.png')))
pred_train <- model %>% predict(xtrain_scaled)
par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
for (o in 1:n_outputs){
  plot(ytrain_scaled[,o], pred_train[,o], xlim=c(-1,1),ylim=c(-1,1)) #main = y_names[o])
}
dev.off()
# plot ANN fit (unscaled)
y_train_unscaled <- outcome #outcome[train_ind,]
pred_train <- model %>% predict(xtrain_scaled)
pred_train_unscaled <- unscale_data(pred_train,ymins,ymaxs) #exp(unscale_data(pred_train,ymins,ymaxs))
png(filename=paste0('output/ann_fit_vs_observed_influenza_4lyrs(unscaled,2LHS,Raw)_',paste0(n_hidden_nodes,'.png')))
par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
for (o in 1:n_outputs){
  plot(y_train_unscaled[,o], pred_train_unscaled[,o], xlim=c(0,max(outcome[,o],pred_train_unscaled[,o])),ylim=c(0,max(outcome[,o],pred_train_unscaled[,o]))) #main = y_names[o])
}
dev.off()


# ======== STAN SECTION ====================
# pass the weights and biases to Stan for Bayesian calibration
n_layers <- length(weights)
weight_first <- weights[[1]]
beta_first <- 1 %*% weights[[2]] 
weight_last <- weights[[n_layers-1]]
beta_last <- 1 %*% weights[[n_layers]]
weight_middle <- array(0, c(n_hidden_layers, n_hidden_nodes, n_hidden_nodes))
beta_middle <- array(0, c(n_hidden_layers, 1, n_hidden_nodes))
for (l in 1:n_hidden_layers){
  weight_middle[l,,] <- weights[[l*2+1]]
  beta_middle[l,,] <- weights[[l*2+2]]
}
stan.dat=list(
  num_hidden_nodes = n_hidden_nodes, 
  num_hidden_layers= n_hidden_layers, 
  num_inputs=n_inputs,
  num_outputs=n_outputs,
  num_targets=1,
  y_maxs = t(as.matrix(ymaxs)),
  y_mins = t(as.matrix(ymins)),
  y_targets = y_targets,
#  y_targets_se = y_targets_se,
  beta_first = beta_first,
  beta_middle = beta_middle,
  beta_last = beta_last, 
  weight_first = weight_first, 
  weight_middle = weight_middle, 
  weight_last = weight_last)
#m <- stan_model("code/post_multi_perceptron_influenza.stan")
#stan.time <- proc.time()
#s <- sampling(m, data = stan.dat, iter = n_iter, chains = n_chains,pars = c("Xq"))
#proc.time() - stan.time # stan sampling time
#fitmat = as.matrix(s)
#Xq <- fitmat[,grep("Xq", colnames(fitmat))]

s0 <- stan(file = "code/post_multi_perceptron_influenza_peak_unscale(PoissonLP).stan", data = stan.dat, iter = n_iter, chains = n_chains, cores = 7, 
           pars = c("Xq") , control = list(adapt_delta = 0.9, max_treedepth = 15))
proc.time()
s_df <- extract(s0)
Xq  <- s_df$Xq
dim(Xq) <- c(n_iter*n_chains/2, n_inputs)
se_post <- s_df$y_targets_se # posterior of targets_se
se_prior <- matrix(rep(abs(rnorm(nrow(se_post),0,0.005)),n_outputs),nrow=n_iter*n_chains/2,n_outputs)
se_post <- cbind(replicate(16,se_post[,1]),replicate(4,se_post[,2])) # if only 2 SEs are sampled

# get ypred using stan output and the Keras ANN
Xqmeans <- colMeans(Xq)
pred_keras <- model %>% predict(Xq)
pred_keras_means <- colMeans(pred_keras)
plot(log(y_targets), pred_keras_means)

# Scale the posteriors
Xq_unscaled <- unscale_data(Xq, vec.mins = xmins, vec.maxs = xmaxs)
pred_y_unscaled<- unscale_data(pred_keras, vec.mins=ymins, vec.max=ymaxs)
Xq_prior <- unscale_data(xscaled,xmins, xmaxs)
Xq_post <- Xq_unscaled
# SAve the unscaled posterior samples
write.csv(Xq_unscaled, file = "output/calibrated_posteriors_influenza_randomSE_100k.csv")
write.csv(exp(pred_y_unscaled), file = "output/posterior_predicted_influenza_randomSE_100k.csv")

# Plot historgram for individual posteriors and compare to prior
plot_prpst <- function(priordf, postdf, type){
  priordf <- data.frame(priordf)
  postdf <- data.frame(postdf)
  if (type == 'Xq'){
    col_names <- x_names
  } else {
    col_names <- y_names}
  colnames(postdf) <- col_names
  colnames(priordf) <- col_names
  priordf$type <- 'prior'
  postdf$type <- paste('chain ', sort(rep(1:n_chains, n_iter/2)))
  priorpost <- rbind(priordf, postdf)
  melt_df <- melt(priorpost)
  #line_df <- data.frame(variable = x_names, intercept = cbind(x_true_unscaled))
  #colnames(line_df) <- c("variable","intercept")
  ggplot(melt_df, aes(value, fill=type, colour = type)) + 
    geom_density(alpha = 0.1) + 
    facet_wrap(~variable, scales="free")#+
  #geom_vline(data=b_param_t, aes(xintercept=V1), color='red')
  #ggsave("output/convergence_influenza_MSE_newBayCANN.png")  
}

b_param_t <- as.data.frame(t(b_param_unscaled))
b_param_t$variable <- x_names
plot_prpst(Xq_prior, Xq_post, 'Xq')

print(s0)
plot(s0)
traceplot(s0,inc_warmup=TRUE, pars=c('Xq'))
ainfo <- get_adaptation_info(s0)
cat(ainfo[[1]])
seed <- get_seed(s0)
sp <- get_sampler_params(s0)
sp2 <- get
# divergence
pairs(s0, pars=c("Xq"))


# Plot histogram for combined posterior and compare to the truth and the prior
bestparam_t <-as.data.frame(t(bestparam))
bestparam_t$variable <- colnames(param)
map <- data.frame(V1 = t(b_param_unscaled), V2 = colMedians(Xq_unscaled), variable =colnames(param))
melt_df_comb = melt_df
melt_df_comb[melt_df_comb == "chain  1" | melt_df_comb == "chain  2" |
               melt_df_comb == "chain  3" | melt_df_comb == "chain  4"| melt_df_comb == "chain  5"| melt_df_comb == "chain  6"| melt_df_comb == "chain  7" ] <- "post"
ggplot(melt_df_comb, aes(value, fill=type, colour = type)) + 
  geom_density(alpha = 0.1) + 
  facet_wrap(~variable, scales="free") + 
  geom_vline(data=bestparam_t, aes(xintercept=V1), color='black')+
  geom_vline(data=map, aes(xintercept=V1),color='red')+
  geom_vline(data=map, aes(xintercept=V2), color='red', linetype='dotted')
ggsave("output/prior_post_truth_influenza.png", height = 10, width = 12)
png(filename='output/joint_post_influenza.png')
pairs(Xq_unscaled, labels=x_names, diag.panel = panel.hist)
dev.off()

xscaled_t <- melt(xscaled)
ggplot(data=xscaled_t)+
  geom_density(aes(x=value))+
  facet_wrap(~variable)+
  scale_x_log10()
