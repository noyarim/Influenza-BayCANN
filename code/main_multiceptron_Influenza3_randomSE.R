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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory
setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/")

# ==================
# Input parameters
n_obs <- 20000
n_iter <- 100000
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


# PSA Samples
runID_list <- c("602","603","607","608","609","610","611","612","613")
param <- read.csv("data/EndOfGenAPMC__AllParamSets_Run#601.csv")
outcome <-read.csv("data/EndOfGenAPMC__AllSimnOutcome_Run#601.csv")
for (runID in runID_list){
  this_param <- read.csv(paste0("data/EndOfGenAPMC__AllParamSets_Run#",paste0(runID,".csv")))
  param <-rbind(param,this_param)
  this_outcome <- read.csv(paste0("data/EndOfGenAPMC__AllSimnOutcome_Run#",paste0(runID,".csv")))
  outcome <-rbind(outcome,this_outcome)
}
colnames(param) <- c("trans.T1","trans.T2","trans.T3","trans.T4","immun1","immun2","immun3","astn12","astn13","astn14","astn15")
colnames(outcome) <- c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15") 
Xunscaled <- param
Yunscaled <- log(outcome)# transformation of outcome data. outcome

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
colnames(targets) <- c("T1.12","T2.12","T3.12","T4.12","T1.13","T2.13","T3.13","T4.13","T1.14","T2.14","T3.14","T4.14","T1.15","T2.15","T3.15","T4.15")
true_targets_mean <- log(targets[1,])
true_targets_se <- true_targets_mean*0.1 #TODO: standard error is arbitrary.

# =====================================
# Scale the targets and their SE
y_targets <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1
y_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
y_targets <- t(as.matrix(y_targets))
y_targets_se <- matrix(0.1,ncol=16,nrow=1) #t(as.matrix(y_targets_se)) 

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
indcol <- 1:16
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
png(filename='output/ann_convergence_influenza_20000_5_MSE.png')
plot(history)
dev.off()
# plot ANN validation
model %>% evaluate(x_test, y_test)
weights <- get_weights(model) #get ANN weights
pred <- model %>% predict(x_test)
png(filename='output/ann_validation_vs_observed_influenza_20000_MSE_5.png')
par("mar", mfrow = c(4,4), mar=c(1,1,1,1))
for (o in 1:n_outputs){
  plot(y_test[,o], pred[,o]) #main = y_names[o])
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

s0 <- stan(file = "code/post_multi_perceptron_influenza_randomSE.stan", data = stan.dat, iter = n_iter, chains = n_chains, 
           pars = c("Xq","y_targets_se") , control = list(adapt_delta = 0.9, max_treedepth=15))
s_df <- extract(s0)
Xq  <- s_df$Xq
dim(Xq) <- c(n_iter*n_chains/2, n_inputs)
se_post <- s_df$y_targets_se # posterior of targets_se
se_prior <- matrix(rep(abs(rnorm(nrow(se_post),0,0.2)),n_outputs),nrow=n_iter*n_chains/2,n_outputs)

# get ypred using stan output and the Keras ANN
Xqmeans <- colMeans(Xq)
pred_keras <- model %>% predict(Xq)
pred_keras_means <- colMeans(pred_keras)
plot(y_targets, pred_keras_means)

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
    facet_wrap(~variable, scales="free")
    #geom_vline(data=bestparam_t, aes(xintercept=V1))
    #ggsave("output/convergence_influenza_MSE_newBayCANN.png")  
}

plot_prpst(Xq_prior, Xq_post, 'Xq')
plot_prpst(se_prior, se_post, 'y_se')

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
