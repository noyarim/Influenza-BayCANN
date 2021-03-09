#################################################################
## Code updated: 2021.2.7########################################
## This code is to debug wrong convergence of BayCANN posterior##
## Run main_multiceptron_Influenza2.R before running this code###
#################################################################
library(matrixStats)
## 1. With the best parameter from the paper (MAP of paper's posterior)
bestparam <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/TheBestParam.csv")
bestparam <- bestparam[,2:ncol(bestparam)] # remove the first column for index
bestparam_scaled <- 2 * (bestparam - xmins) / (xmaxs - xmins) - 1 # scale x
# 1-1. outcome from PSA
bestpred_psa <- targets
bestpred_seir <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/bestparam_seir_outcome.csv")
# 1-2. outcome predicted by ANN
bestpred <- model %>% predict(as.matrix(bestparam_scaled))
bestpred_unscaled <- unscale_data(bestpred, vec.mins=ymins, vec.maxs=ymaxs)
bestpred_unscaled_exp <- exp(bestpred_unscaled)
# compare outcomes with the best param from SEIR and ANN
bestpred_all <- rbind(bestpred_psa,bestpred_unscaled_exp)

## 2. With a random parameter from PSA samples
rparam <- param[1,]
rparam_scaled <- 2 * (rparam - xmins) / (xmaxs - xmins) - 1 # scale x
# 2-1. outcome from PSA
routcome <-outcome[1,]
# 2-2. outcome predicted by ANN
rpred <- model %>% predict(as.matrix(rparam_scaled))
rpred_unscaled <- unscale_data(rpred, vec.mins=ymins, vec.maxs = ymaxs)
rpred_unscaled_exp <- exp(rpred_unscaled)
colnames(rpred_unscaled_exp) <- colnames(routcome)
# compare the outcomes with a random param from SEIR and ANN
rpred_all <- rbind(routcome, rpred_unscaled_exp)

## 3. With the MAP of BayCANN posterior
s0 <- stan(file = "code/post_multi_perceptron_influenza.stan", data = stan.dat, iter = n_iter, chains = n_chains, 
           pars = c("Xq") , control = list(adapt_delta = 0.99))
s_df <- extract(s0, permuted = TRUE)
Xq  <- s_df$Xq # parameters in BayCANN posterior
se <- s_df$y_targets_se
dim(Xq) <- c(n_iter*n_chains/2, n_inputs)
stan_lp <- s_df$lp__ # log density in posterior
best_idx <- which(stan_lp == max(stan_lp))
b_param <- Xq[best_idx,] # MAP of BayCANN (scaled)
b_se <- se[best_idx,]
b_se <- c(rep(b_se[1],16),rep(b_se[2],4)) # if fitting only 2 SEs
b_param_unscaled <- unscale_data(t(as.matrix(b_param)), xmins, xmaxs)
b_pred_ann <- model %>% predict(t(as.matrix(b_param)))
b_pred_ann_unscaled <- unscale_data(b_pred_ann, vec.mins = ymins, vec.maxs = ymaxs)
b_pred_ann_unscaled_exp <- exp(b_pred_ann_unscaled)
# 3-1. outcome simulated in SEIR
b_pred_seir <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/MAP_BayCANN_seiroutcome.csv")
b_pred_seir <- as.vector(b_pred_seir$SimnRate)
# 3-2. outcome predicted by ANN
write.csv(b_param_unscaled, "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/BayCANN_MAP_2targets_2.csv")
write.csv(b_pred_ann_unscaled_exp, "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/BayCANN_MAP_2targets_annoutcome.csv")
## 4. With the median of BayCANN posterior
median_X <- colMedians(Xq)
median_X_unscaled <- unscale_data(as.matrix(t(median_X)), vec.mins=xmins, vec.maxs=xmaxs)
# 4-1. outcome simulated in SEIR
m_pred_seir <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/med_BayCANN_seiroutcome.csv")
m_pred_seir <- as.vector(m_pred_seir$SimnRate)
pred_med <- model %>% predict(t(as.matrix(median_X)))
pred_med_unscaled <- unscale_data(pred_med, vec.mins = ymins, vec.maxs = ymaxs)
pred_med_unscaled_exp <-exp(pred_med_unscaled)

## Compile and compare
all <- rbind(bestpred_all,b_pred_seir, b_pred_ann_unscaled_exp,m_pred_seir,pred_med_unscaled_exp)
rownames(all) <- c("best_SEIR(target)","best_ANN","MAP_SEIR","MAP_ANN","median_SEIR","median_ANN")
all_t <- melt(all)

## Plots
# 1. check whether PSA samples cover the target
priordf <- data.frame(outcome)
postdf <- data.frame(targets)
colnames(postdf) <- colnames(priordf)
#priordf$type <- 'psa'
#postdf$type <- paste('chain ', sort(rep(1:n_chains, n_iter/2)))
#priorpost <- rbind(priordf, postdf)
melt_df <- melt(priordf)
postdf_t <- melt(postdf)
#line_df <- data.frame(variable = x_names, intercept = cbind(x_true_unscaled))
#colnames(line_df) <- c("variable","intercept")
ggplot(melt_df, aes(value)) + 
  geom_density(alpha = 0.1) + 
  facet_wrap(~variable, scales="free") + 
  geom_vline(data=postdf_t, aes(xintercept=value),color='red')+
  scale_x_log10()
ggsave("output/convergence_influenza_bestparam.png")

# 2. comparison of outcomes
ggplot(data = all_t, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(position = "dodge",stat = "identity", width = 0.6)+
  xlab("outcomes")+
  ylab("rates (per 100,000)")+
  scale_y_log10(expand = c(0,0))+
  theme(
    # backgrounds
    panel.grid.major = element_blank(), panel.grid.minor= element_blank(), panel.background = element_blank(), axis.line=element_line(colour = "black"), legend.key = element_blank(),
    # title
    plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 15),
    legend.title=element_blank(),
    legend.text=element_text(size = 10),
    legend.position = 'bottom',#c(0.1,0.8),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 30, margin = margin(t=10)),
    axis.text.y = element_text(size = 11),
    strip.text.x = element_text(size = 11))
ggsave("output/outcome_comparison.png")
# 3. posterior of Ed Hill's paper
ed_param <- read.table("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/SeasonalFluImmunityPropagation-master/Results/ParamInference/FourSeasonFitData/ParameterSets_FourSeasonFit.txt")
ed_param_t <- melt(ed_param)
bestparam_t <- melt(bestparam) # MAP of the paper's posterior
map_t <- melt(as.data.frame(b_param_unscaled))
ggplot(ed_param_t, aes(value)) + 
  geom_density(alpha = 0.1) + 
  facet_wrap(~variable, scales="free") +
  geom_vline(data=bestparam_t, aes(xintercept=value),color='red')+
  geom_vline(data=map_t, aes(xintercept=value), color='blue')
ggsave("output/paper_postdist2.png")

ggpairs(ed_param)

## 4. comparison of paper's MAP, new BayCANN MAP with the observed targets
targets <- read.csv("data/EndOfGenAPMC_Target_vector.csv")
bestpred_seir <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/bestparam_seir_outcome.csv")
newbaycan_seir <- read.csv("data/newBayCANN_best(MSE,redSE)_outcome.csv") # 1target seir outcome
newbaycan_ann <- read.csv("data/newBayCANN_best(MSE,redSE)_annoutcome.csv") # 1target ann outcome
newbaycan2_seir <- read.csv("data/BayCANN_MAP_2targets_seiroutcome.csv") # 2targets seir outcome
newbaycan2_ann <- read.csv("data/BayCANN_MAP_2targets_annoutcome.csv") # 2targets ann outcome
#med_seir <- read.csv("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/data/newBayCANN_median_outcome.csv")
all <- as.data.frame(cbind(targets$x,bestpred_seir$x,newbaycan_seir$x, newbaycan_ann$x, newbaycan2_seir$x[1:16], newbaycan2_ann$x[1:16]))
all <- as.data.frame(cbind(targets$x,bestpred_seir$x,b_pred_ann_unscaled_exp[1:16]))
colnames(all) <- c("targets","MAP_EdHill","1target_MAP_seir","1target_MAP_ann","2targets_MAP_seir","2targets_MAP_ann")
colnames(all) <- c("targets", "MAP_EdHill","MAP_BayCANN")
all$subtype <- rep(c("H1N1","H3N2","B/Y","B/V"),4)
all$season <- c(rep("2012/13",4),rep("2013/14",4),rep("2014/15",4),rep("2015/16",4))
all$subtype <- factor(all$subtype, levels=c("H1N1","H3N2","B/V","B/Y"))
all_t <- melt(all)
ggplot(data = all_t, aes(x = season, y = value, fill = variable)) +
  geom_bar(position = "dodge",stat = "identity", width = 0.6)+
  facet_wrap(~subtype)+
  xlab("outcomes")+
  ylab("rates (per 100,000)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,100))+
  theme(
    # backgrounds
    panel.grid.major = element_blank(), panel.grid.minor= element_blank(), panel.background = element_blank(), axis.line=element_line(colour = "black"), legend.key = element_blank(),
    # title
    plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 15),
    legend.title=element_blank(),
    legend.text=element_text(size = 10),
    legend.position = 'bottom',#c(0.1,0.8),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 30, margin = margin(t=10)),
    axis.text.y = element_text(size = 11),
    strip.text.x = element_text(size = 11))+
  scale_fill_manual(values=c("red","orange","green4","green1","royalblue4","royalblue1"))
ggsave("output/outcome_comparison.png")


# Draw distribution around target with MAP y_targets_se
targets <- read.csv("data/EndOfGenAPMC_Target_vector.csv")
targets <- t(targets$x) 
targets <- cbind(targets, 50,50,50,50) #cbind(targets, 182,182,182,182) # 1 is placeholder for peak date
colnames(targets) <- c("T1.12","T2.12","T3.12","T4.12","T1.13","T2.13","T3.13","T4.13","T1.14","T2.14","T3.14","T4.14","T1.15","T2.15","T3.15","T4.15","peak1","peak2","peak3","peak4" )
#true_targets_mean <- targets
true_targets_mean <- log(targets) #c(log(targets[1,1:ncol(targets)-1]),targets[1,ncol(targets)])
# scale
g_star <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1
n_ts <- 1000 
target_samples <- matrix(0,n_ts,n_outputs)
for (i in c(1:n_outputs)){
  target_samples[,i] <- rlnorm(1000,g_star[1,i],b_se[i])
}
colnames(target_samples) <- colnames(targets)
target_samples_t <- melt(as.data.frame(target_samples))
ggplot(target_samples_t, aes(value)) + 
  geom_density() + 
  facet_wrap(~variable, scales="free")
org_target_samples <- log(target_samples)
unsc_target_samples <- unscale_data(org_target_samples, ymins, ymaxs)
fin_target_samples <- exp(unsc_target_samples)
fin_target_samples_t <- melt(as.data.frame(fin_target_samples))
fin_target_samples_t$targets <- rep(targets, each=n_ts)
ggplot(fin_target_samples_t, aes(value)) + 
  geom_density() + 
  facet_wrap(~variable, scales="free")+
  geom_vline(aes(xintercept = targets),color='red')

## BayCANN with parameter_6 = 1
s0 <- stan(file = "code/post_multi_perceptron_influenza_test.stan", data = stan.dat, iter = n_iter, chains = 1, 
           pars = c("Xq_v")) #, control = list(adapt_delta = 0.99))

s_df <- extract(s0)
Xq  <- s_df$Xq
