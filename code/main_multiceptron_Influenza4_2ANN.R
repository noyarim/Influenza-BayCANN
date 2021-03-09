# Code last updated on: 2021. 02.22
# Original code: main_multiceptron_Influenza4.R
# Change made: 
# (1) Run this code for 2-target simulation
# (2) This code fits PSA samples to 2 outcomes with 2 ANNs separately
# (3) BayCANN is running with 1 of 2 ANNs and the second ANN is used to reject samples with the second outcome(peak) > 0.5

#install.packages("keras")
#install_keras()
#install_keras(tensorflow = "gpu")
#install.packages("curl")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("BioStatR")
library(keras)
library(rstan)
library(ggplot2)
library(reshape2) # for melt()
library(BioStatR)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory
setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/")

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
n_chains <- 4

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

outcome3_f <- outcome3 %>%
  mutate(
    accept = (1-V1_bi)*(1-V2_bi)*(1-V3_bi)*(1-V4_bi)
  ) %>%
  filter(
    accept==1
  )

# PSA Samples
runID_list <- c("701","702","703","704","705","706","707","708","709","710")
param <- read.csv("data/PSA_2targets/EndOfGenAPMC__AllParamSets_Run#700.csv")
outcome1 <-read.csv("data/PSA_2targets/EndOfGenAPMC__AllSimnOutcome_Run#700.csv")
outcome2 <-read.csv("data/PSA_2targets/EndOfGenAPMC__AllSimnPeakAll_Run#700.csv") # 1 binary variable
outcome3 <-read.csv("data/PSA_2targets/EndOfGenAPMC__AllSimnPeakDate_Run#700.csv") # peak days in 4 seasons
for (runID in runID_list){
  print(runID)
  this_param <- read.csv(paste0("data/PSA_2targets/EndOfGenAPMC__AllParamSets_Run#",paste0(runID,".csv")))
  param <-rbind(param,this_param)
  this_outcome1 <- read.csv(paste0("data/PSA_2targets/EndOfGenAPMC__AllSimnOutcome_Run#",paste0(runID,".csv")))
  outcome1 <-rbind(outcome1,this_outcome1)
  this_outcome2 <- read.csv(paste0("data/PSA_2targets/EndOfGenAPMC__AllSimnPeakAll_Run#",paste0(runID,".csv")))
  outcome2 <-rbind(outcome2,this_outcome2)  
  this_outcome3 <- read.csv(paste0("data/PSA_2targets/EndOfGenAPMC__AllSimnPeakDate_Run#",paste0(runID,".csv")))
  outcome3 <-rbind(outcome3,this_outcome3)
}
# calculate binary variables based on peak dates in 4 seasons
outcome3<-outcome3 %>% mutate(
  V1_bi = ifelse(V1 >= 182,1,0),
  V2_bi = ifelse(V2 >= 182,1,0),
  V3_bi = ifelse(V3 >= 182,1,0),
  V4_bi = ifelse(V4 >= 182,1,0)
)
outcome <- cbind(outcome1, outcome3) #cbind(outcome1, outcome2)
colnames(param) <- c("trans.T1","trans.T2","trans.T3","trans.T4","immun1","immun2","immun3","astn12","astn13","astn14","astn15")
colnames(outcome) <- c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15","peak1","peak2","peak3","peak4","peak1_bi","peak2_bi","peak3_bi","peak4_bi") 
#c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15","peak") 
Xunscaled <- param
Yunscaled <- cbind(log(outcome[,1:ncol(outcome1)]),outcome[,(ncol(outcome1)+1):ncol(outcome)])# transformation of outcome data. outcome

# scale the PSA inputs and outputs
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

# get the true data (parameter used to generate target data)
x_true_unscaled <- read.csv(file="/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/TheBestParam.csv") 
x_true_unscaled <- x_true_unscaled[2:ncol(x_true_unscaled)]

# get the target data
targets <- read.csv("data/EndOfGenAPMC_Target_vector.csv")
#targets <- read.csv("data/SimnOutcomewithMedian.csv")
#targets <- read.csv("data/SimnOutcomewithTheBest.csv") # test with the best parameter set
#targets <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/MAP_BayCANN_seiroutcome.csv")
#targets <- t(targets$SimnRate)
targets <- t(targets$x) # test with the best parameter set
targets <- cbind(targets,0)
colnames(targets) <- c("T1.12","T2.12","T3.12","T4.12","T1.13","T2.13","T3.13","T4.13","T1.14","T2.14","T3.14","T4.14","T1.15","T2.15","T3.15","T4.15","peak")
true_targets_mean <- c(log(targets[1,1:ncol(targets)-1]),targets[1,ncol(targets)])
true_targets_se <- true_targets_mean*0.1 #TODO: standard error is arbitrary.

# =====================================
# Scale the targets and their SE
y_targets <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1
#y_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
y_targets <- t(as.matrix(y_targets))
y_targets_se <- matrix(0.05,ncol=ncol(y_targets),nrow=1) #t(as.matrix(y_targets_se)) #matrix(0.1,ncol=16,nrow=1) 

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
indcol <- 1:(num_outputs-1) # 16 instead of 17
x_train <- data.matrix(X)
y_train <- data.matrix(Y[,indcol])
x_test <- data.matrix(Xt)
y_test <- data.matrix(yt[,indcol])
n_outputs <- dim(y_test)[2]
n_inputs <- dim(x_test)[2]
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
proc.time() - keras.time #keras ann fitting time

# plot ANN convergence
png(filename='output/ann_convergence_influenza_2targ.png')
plot(history)
dev.off()
# plot ANN validation
model %>% evaluate(x_test, y_test)
weights <- get_weights(model) #get ANN weights
pred <- model %>% predict(x_test)
png(filename='output/ann_validation_vs_observed_influenza_2targ.png')
par("mar", mfrow = c(4,5), mar=c(1,1,1,1))
for (o in 1:n_outputs){
  plot(y_test[,o], pred[,o]) #main = y_names[o])
}
dev.off()

# ============== TensorFlow Keras ANN Section (for the second output) ========================
N2 <- floor(n_obs * .7)
train_ind2 <- sample(n_obs,N2)
test_ind2 <- setdiff(1:n_obs, train_ind2) 
indcol2 <- 21:ncol(outcome) # 16 instead of 17
x_train2 <- data.matrix(xscaled[train_ind2,])
x_test2 <- data.matrix(xscaled[test_ind2,])
y_train2 <- data.matrix(Yunscaled[train_ind2,indcol2])
y_test2 <- data.matrix(Yunscaled[test_ind2,indcol2])
n_outputs2 <- dim(y_test2)[2]
model2 <- keras_model_sequential() 
model2 %>% 
  layer_dense(units = n_hidden_nodes, activation = 'tanh', input_shape = n_inputs) %>% 
  layer_dense(units = n_hidden_nodes, activation = 'tanh') %>%
  layer_dense(units = n_outputs2, activation ='sigmoid')
summary(model2)
model2 %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam'
)
keras.time <- proc.time()

history2 <- model2 %>% fit(
  x_train2, y_train2, 
  epochs = n_epochs, batch_size = N2, 
  validation_split = validation_split,
  verbose = verbose
)
proc.time() - keras.time #keras ann fitting time

# plot ANN convergence
png(filename='output/ann_convergence_influenza_2targ.png')
plot(history2)
dev.off()
# plot ANN validation
model2 %>% evaluate(x_test2, y_test2)
weights2 <- get_weights(model2) #get ANN weights
pred2 <- model2 %>% predict(x_test2)
png(filename='output/ann_validation_vs_observed_influenza_2targ.png')
par("mar", mfrow = c(2,2))#, mar=c(1,1,1,1))
for (o in 1:n_outputs2){
  plot(y_test2[,o], pred2[,o]) #main = y_names[o])
}
dev.off()

hist(y_train2[,1],y_train2[,2],y_train2[,3],y_train2[,4])
y_train2_t <- melt(as.data.frame(y_train2))
ggplot(y_train2_t, aes(x=value))+
  geom_histogram()+
  facet_wrap(~variable)

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
  y_targets = y_targets,
  y_targets_se = y_targets_se,
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

s0 <- stan(file = "code/post_multi_perceptron_influenza.stan", data = stan.dat, iter = n_iter, chains = n_chains, 
           pars = c("Xq") , control = list(adapt_delta = 0.9))
s_df <- extract(s0)
Xq  <- s_df$Xq
dim(Xq) <- c(n_iter*n_chains/2, n_inputs)

# get ypred using stan output and the Keras ANN
Xqmeans <- colMeans(Xq)
pred_keras <- model %>% predict(Xq)
pred_keras_means <- colMeans(pred_keras)
plot(y_targets, pred_keras_means)

# Scale the posteriors
Xq_unscaled <- unscale_data(Xq, vec.mins = xmins, vec.maxs = xmaxs)
pred_y_unscaled<- unscale_data(pred_keras, vec.mins=ymins, vec.max=ymaxs)

# SAve the unscaled posterior samples
write.csv(Xq_unscaled, file = "output/calibrated_posteriors_influenza_NewBayCANN_2targets.csv")
write.csv(exp(pred_y_unscaled), file = "output/posterior_predicted_influenza_NewBayCANN_2targets.csv")

# Plot historgram for individual posteriors and compare to prior
priordf <- data.frame(Xunscaled)
postdf <- data.frame(Xq_unscaled)
colnames(postdf) <- x_names
priordf$type <- 'prior'
postdf$type <- paste('chain ', sort(rep(1:n_chains, n_iter/2)))
priorpost <- rbind(priordf, postdf)
melt_df <- melt(priorpost)
#line_df <- data.frame(variable = x_names, intercept = cbind(x_true_unscaled))
#colnames(line_df) <- c("variable","intercept")
ggplot(melt_df, aes(value, fill=type, colour = type)) + 
  geom_density(alpha = 0.1) + 
  facet_wrap(~variable, scales="free") 
  #geom_vline(data=bestparam_t, aes(xintercept=V1))
ggsave("output/convergence_influenza_MSE_newBayCANN.png", height = 5, width = 10)

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
