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
setwd("/home/kel171/Influenza/")

# ==================
# Input parameters
n_obs <- 22000
n_iter <- 10000
n_hidden_nodes <- 100
n_hidden_layers <- 5 - 1 # technically it is easier if we just subtract 1 frmo the numbrer of layers
n_epochs <- 10000
verbose <- 0
n_batch_size <- 2000
validation_split <- 0.2
n_chains <- 7

# =================================================
# functions
scale_data <- function(unscaled_data){
  vec.maxs <- apply(unscaled_data, 2, max) 
  vec.mins <- apply(unscaled_data, 2, min)
  vec.ones <- matrix(1, nrow = nrow(unscaled_data), 1)
  mat.maxs <- vec.ones %*% vec.maxs
  mat.mins <- vec.ones %*% vec.mins
  scaled_data <- 2 * (unscaled_data - mat.mins) / (mat.maxs - mat.mins) - 1
  results <- list(scaled_data = scaled_data, vec.mins = vec.mins, vec.maxs = vec.maxs)
  return(results)
}
unscale_data <- function(scaled_data, vec.mins, vec.maxs){
  vec.ones <- matrix(1, nrow = nrow(scaled_data), 1)
  mat.mins <- vec.ones %*% vec.mins
  mat.maxs <- vec.ones %*% vec.maxs
  unscaled_data <- (scaled_data + 1) * (mat.maxs - mat.mins) / 2 + mat.mins
}


# PSA Samples
runID_list <- c(611,612,613,614,615,616,617,618,619,620)
param <- read.csv("output/PSA_output/EndOfGenAPMC__AllParamSets_Run#621.csv")
outcome1 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnOutcome_Run#621.csv")
outcome2 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakAll_Run#621.csv") # 1 binary variable
outcome3 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakDate_Run#621.csv") # peak days in 4 seasons
for (runID in runID_list){
  print(runID)
  this_param <- read.csv(paste0("output/PSA_output/EndOfGenAPMC__AllParamSets_Run#",paste0(runID,".csv")))
  param <-rbind(param,this_param)
  this_outcome1 <- read.csv(paste0("output/PSA_output/EndOfGenAPMC__AllSimnOutcome_Run#",paste0(runID,".csv")))
  outcome1 <-rbind(outcome1,this_outcome1)
  this_outcome2 <- read.csv(paste0("output/PSA_output/EndOfGenAPMC__AllSimnPeakAll_Run#",paste0(runID,".csv")))
  outcome2 <-rbind(outcome2,this_outcome2)
  this_outcome3 <- read.csv(paste0("output/PSA_output/EndOfGenAPMC__AllSimnPeakDate_Run#",paste0(runID,".csv")))
  outcome3 <-rbind(outcome3,this_outcome3)
}
outcome <- cbind(outcome1, outcome3) #cbind(outcome1, outcome2)
colnames(param) <- c("trans.T1","trans.T2","trans.T3","trans.T4","immun1","immun2","immun3","astn12","astn13","astn14","astn15")
colnames(outcome) <- c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15","peak1","peak2","peak3","peak4") 
# remove the samples peak == 366 (max) and == 1(min)
keep_idx <- which(outcome[,17]<366&outcome[,18]<366&outcome[,19]<366&outcome[,20]<366&outcome[,17]!=1&outcome[,18]!=1&outcome[,19]!=1&outcome[,20]!=1)
outcome <- outcome[keep_idx,]
param <- param[keep_idx,]
n_obs <- nrow(outcome) # update number of psa samples
# outcome transformation
# # if using log
# Yunscaled <- log(outcome)# log transformation of outcome data. 
# #if using box-cox transformation
# outcome_trans = outcome
# bc_lambda <- rep(NA, ncol(outcome))
# for (i in 1:ncol(outcome)){
#   temp <- as.data.frame(cbind(outcome[,i],param))
#   lm_model <- lm(outcome[,i]~param[,3], data=temp)
#   bc <- boxcox(lm_model)
#   bc_lambda[i] <- bc$x[which.max(bc$y)]
#   outcome_trans[,i] = (outcome[,i]^bc_lambda[i]-1)/bc_lambda[i]
# }
# # comparison of log transformation and box-cox transformation
# par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
# for (o in 1:n_outputs){
#   hist(log(outcome[,o]))
# }
# par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
# for (o in 1:n_outputs){
#   hist(outcome_trans[,o])
# }
# log transformation of outcome data. 
Yunscaled <- log(outcome)
#Yunscaled <- cbind(log(outcome[,1:ncol(outcome1)]),outcome[,(ncol(outcome1)+1):ncol(outcome)]) # logtransformation outcomes except peak dates
Xunscaled <- param

# scale the PSA inputs and outputs for ANN (this is for ANN)
xresults <- scale_data(Xunscaled) 
yresults <- scale_data(Yunscaled)
xscaled <- xresults$scaled_data 
yscaled <- yresults$scaled_data 
xmins <- xresults$vec.mins
xmaxs <- xresults$vec.maxs
ymins <- yresults$vec.mins
ymaxs <- yresults$vec.maxs
y_names <-colnames(outcome)#colnames(out_i_det_unif)
x_names <-colnames(param)#colnames(samp_i_unif)

# get the target data
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
N <- floor(n_obs * .8)
Nt <- n_obs-N
train_ind <- sample(n_obs,N)
test_ind <- setdiff(1:n_obs, train_ind) 
X <- xscaled[train_ind,]
Y <- yscaled[train_ind,]
Xt <- xscaled[test_ind,]
yt <- yscaled[test_ind,]
num_outputs <- ncol(Y)
indcol <- 1:num_outputs
x_train <- data.matrix(X)
y_train <- data.matrix(Y[,indcol])
x_test <- data.matrix(Xt)
y_test <- data.matrix(yt[,indcol])
n_outputs <- dim(y_test)[2]
n_inputs <- dim(x_test)[2]

compile_and_fit <- function(n_hidden_nodes){

  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = n_hidden_nodes, activation = 'tanh', input_shape = n_inputs) %>% 
    layer_dense(units = n_hidden_nodes, activation = 'tanh') %>%
    layer_dense(units = n_hidden_nodes, activation = 'tanh') %>%
    layer_dense(units = n_hidden_nodes, activation = 'tanh') %>%
    layer_dense(units = n_hidden_nodes, activation = 'tanh') %>%
    layer_dense(units = n_outputs)
  summary(model)
  model %>% compile(
    loss = 'mean_squared_error',
    optimizer = 'adam'
  )
  keras.time <- proc.time()
  history <- model %>% fit(
    x_train, y_train, 
    epochs = n_epochs, batch_size = N, 
    validation_split = validation_split,
    verbose = verbose
  )
  weights <- get_weights(model) #get ANN weights
  pred <- model %>% predict(x_test)
  proc.time() - keras.time #keras ann fitting time


  # plot ANN convergence
  png(filename=paste0('output/ann_convergence_influenza_4lyrs_',paste0(n_hidden_nodes,'.png')))
  plot(history)
  dev.off()
  
  loss <- history$metrics$loss
  write.csv(loss, paste0('output/ann_convergence(loss)_4lyrs_',paste0(n_hidden_nodes,'.csv')))
  val_loss <- history$metrics$val_loss
  write.csv(val_loss, paste0('output/ann_convergence(valloss)_4lyrs_',paste0(n_hiden_nodes,'.csv')))
  
  # plot ANN validation
  model %>% evaluate(x_test, y_test)
  png(filename=paste0('output/ann_validation_vs_observed_influenza_4lyrs_',paste0(n_hidden_nodes,'.png')))
  par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
  for (o in 1:n_outputs){
    plot(y_test[,o], pred[,o], xlim=c(-1,1),ylim=c(-1,1)) #main = y_names[o])
  }
  dev.off()
  # plot ANN fit (on scale)
  png(filename=paste0('output/ann_fit_vs_observed_influenza_4lyrs_',paste0(n_hidden_nodes,'.png')))
  pred_train <- model %>% predict(x_train)
  par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
  for (o in 1:n_outputs){
    plot(y_train[,o], pred_train[,o], xlim=c(-1,1),ylim=c(-1,1)) #main = y_names[o])
  }
  dev.off()
  # plot ANN fit (unscaled)
  y_train_unscaled <- outcome[train_ind,]
  pred_train <- model %>% predict(x_train)
  pred_train_unscaled <- exp(unscale_data(pred_train,ymins,ymaxs))
  png(filename=paste0('output/ann_fit_vs_observed_influenza_4lyrs_',paste0(n_hidden_nodes,'.png')))
  par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
  for (o in 1:n_outputs){
    plot(y_train_unscaled[,o], pred_train_unscaled[,o], xlim=c(0,max(y_train_unscaled[,o],pred_train_unscaled[,o])),ylim=c(0,max(y_train_unscaled[,o],pred_train_unscaled[,o]))) #main = y_names[o])
  }
  dev.off()
}

n_hidden_nodes_list <- c(50,100,200,400,1000)
for (n_hidden_nodes in n_hidden_nodes_list){
  print(n_hidden_nodes)
  compile_and_fit(n_hidden_nodes)
}

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
Xq_prior <- Xunscaled
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
    facet_wrap(~variable, scales="free")+
  geom_vline(data=b_param_t, aes(xintercept=V1), color='red')
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
