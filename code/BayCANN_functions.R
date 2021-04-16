setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/Influenza-github/")

# functions

# function to transform and scale the PSA samples
prepare_data <- function(xtrain, ytrain, xtest, ytest, log){
  xtrain <- xtrain[,2:ncol(xtrain)]
  ytrain <- ytrain[,2:ncol(ytrain)]
  xtest <- xtest[,2:ncol(xtest)]
  ytest <- ytest[,2:ncol(ytest)]
  if (log == 1){
    ytrain = log(ytrain)
    ytest = log(ytest)
  }
  y_names <- colnames(ytrain)
  x_names <- colnames(xtrain)
  n_train <- nrow(xtrain)
  n_test <- nrow(xtest)
  x <- rbind(xtrain, xtest)
  y <- rbind(ytrain, ytest)
  n <- nrow(x)
  n_inputs <- length(x_names)
  n_outputs <- length(y_names)
  # scale the PSA inputs and outputs
  xresults <- scale_data(x) 
  yresults <- scale_data(y)
  xscaled <- xresults$scaled_data 
  yscaled <- yresults$scaled_data 
  xmins <- xresults$vec.mins
  xmaxs <- xresults$vec.maxs
  ymins <- yresults$vec.mins
  ymaxs <- yresults$vec.maxs
  
  xtrain_scaled <- xscaled[1:n_train, ]
  ytrain_scaled <- yscaled[1:n_train, ]
  xtest_scaled  <- xscaled[(n_train+1):n, ]
  ytest_scaled  <- yscaled[(n_train+1):n, ]
  
  return(list(n_inputs = n_inputs,
              n_outputs = n_outputs,
              n_train = n_train,
              n_test = n_test,
              x_names = x_names, 
              y_names = y_names,
              xscaled = xscaled,
              yscaled = yscaled,
              xtrain_scaled = xtrain_scaled,
              ytrain_scaled = ytrain_scaled,
              xtest_scaled  = xtest_scaled ,
              ytest_scaled  = ytest_scaled,
              xmins = xmins,
              xmaxs = xmaxs,
              ymins = ymins,
              ymaxs = ymaxs
  ))
}

# function to scale the data
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

# function to unscale the data
unscale_data <- function(scaled_data, vec.mins, vec.maxs){
  vec.ones <- matrix(1, nrow = nrow(scaled_data), 1)
  mat.mins <- vec.ones %*% vec.mins
  mat.maxs <- vec.ones %*% vec.maxs
  unscaled_data <- (scaled_data + 1) * (mat.maxs - mat.mins) / 2 + mat.mins
}

# function to compile PSA samples and save
param <- read.csv("output/PSA_output/EndOfGenAPMC__AllParamSets_Run#621.csv")
outcome1 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnOutcome_Run#621.csv")
outcome2 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakAll_Run#621.csv") # 1 binary variable
outcome3 <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakDate_Run#621.csv") # peak days in 4 seasons
outcome <- cbind(outcome1, outcome3) #cbind(outcome1, outcome2)
colnames(param) <- c("trans.T1","trans.T2","trans.T3","trans.T4","immun1","immun2","immun3","astn12","astn13","astn14","astn15")
colnames(outcome) <- c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15","peak1","peak2","peak3","peak4") 
# remove the samples peak == 366 (max) and == 1(min)
keep_idx <- which(outcome[,17]<366&outcome[,18]<366&outcome[,19]<366&outcome[,20]<366&outcome[,17]!=1&outcome[,18]!=1&outcome[,19]!=1&outcome[,20]!=1)
outcome <- outcome[keep_idx,]
param <- param[keep_idx,]
n_obs <- nrow(outcome) # update number of psa samples
write.csv(param, "output/PSA_output/x_train.csv")
write.csv(outcome, "output/PSA_output/y_train.csv")
# for testing
param_t <- read.csv("output/PSA_output/EndOfGenAPMC__AllParamSets_Run#800.csv")
outcome1_t <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnOutcome_Run#800.csv")
outcome2_t <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakAll_Run#800.csv") # 1 binary variable
outcome3_t <-read.csv("output/PSA_output/EndOfGenAPMC__AllSimnPeakDate_Run#800.csv") # peak days in 4 seasons
outcome_t <- cbind(outcome1_t, outcome3_t) #cbind(outcome1, outcome2)
colnames(param_t) <- c("trans.T1","trans.T2","trans.T3","trans.T4","immun1","immun2","immun3","astn12","astn13","astn14","astn15")
colnames(outcome_t) <- c("rate.T1.12","rate.T2.12","rate.T3.12","rate.T4.12","rate.T1.13","rate.T2.13","rate.T3.13","rate.T4.13","rate.T1.14","rate.T2.14","rate.T3.14","rate.T4.14","rate.T1.15","rate.T2.15","rate.T3.15","rate.T4.15","peak1","peak2","peak3","peak4") 
# remove the samples peak == 366 (max) and == 1(min)
keep_idx <- which(outcome_t[,17]<366&outcome_t[,18]<366&outcome_t[,19]<366&outcome_t[,20]<366&outcome_t[,17]!=1&outcome_t[,18]!=1&outcome_t[,19]!=1&outcome_t[,20]!=1)
outcome_t <- outcome_t[keep_idx,]
param_t <- param_t[keep_idx,]
write.csv(param_t, "output/PSA_output/x_test.csv")
write.csv(outcome_t, "output/PSA_output/y_test.csv")