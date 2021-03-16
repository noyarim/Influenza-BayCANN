#Purpose:
#Script to run Adaptive Population Monte Carlo Approximate Bayesian Computation
#Run serially & take arguments from command line as function inputs

#Fit multi-strain, non-age structured influenza transmission model
#(with immunity propagation) to empirical data

#Code Author: Ed Hill
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
#--------------------------------------------------------------------------
source("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Code/RunSeasonalFluModelODEs.R")
source("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Code/ExpHistUpdate.R")
source("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Code/APMC_loop_BayCANN_Peak.R")
source("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Code/SeasonalFluModel_SerialInferenceFns_BayCANN_peak(LHS).R")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
library(Matrix)
library(stats)
library(deSolve)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(matlab)
library(MASS)
library(mvtnorm)


#SetUpInferenceAPMC <- function(ARGS){

  #Take command line arguments, ARGS, assign to variable names
  #ARGS input listing
  # ARG 1: RunID
  # ARG 2: Observed data filename
  # ARG 3: SeasonsToSimulate
  # ARG 4: N_alpha
  # ARG 5: alpha
  # ARG 6: MinAcceptRate
  # ARG 7: MaxGen
  # ARG 8: PerturbVarFn
  # ARG 9: PriorFn
  # ARG 10: SummStatFn
  # ARG 11: SampleFromPriorFn
  # ARG 12: ModelSimnFn
  # ARG 13: FirstGenFromFileFlag
  
  #To convert strings to numbers, use parse
  ARGS = commandArgs(trailingOnly = TRUE)
  ##delete this lines after completing testing
  ARGS = c("611","/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/SeasonalFluImmunityPropagation-master/Data/ILIData/EmpData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv","7", "1000", "0.5", "-1", "10", "OLCMPerturbVarFn" ,"Prior_FourSeasonFit", "SummStatFun", "SampleFirstGenFn_FourSeasonFit", "RunModelSimn", "0")
  
  #--------------------------------------------------------------------------
  # Set RunID idx
  #--------------------------------------------------------------------------
  RunID = as.integer(ARGS[1])
  
  #--------------------------------------------------------------------------
  # Specify empirical data being fitted to
  #--------------------------------------------------------------------------
  ObvsDataFileName = ARGS[2]
  ObvsDataTemp = read.table(ObvsDataFileName,sep=',')
  
  #Specify number of influenza seasonsto be simulated
  #From 2009/2010 onwards: - > 7 will be up to and inc. 2015/16
  #                        - > 8 up to and inc. 2016/17 etc
  SeasonsToSimulate = as.integer(ARGS[3])
  
  #Get maximum number of seasons that could be used in inference procedure
  #In case all seasons are not used, get number of seasons omitted
  MaxSeasonNumToSimulate = nrow(ObvsDataTemp)
  SeasonsUnused =  MaxSeasonNumToSimulate - SeasonsToSimulate
  
  #Error check
  IntCheck = (SeasonsToSimulate%%1)==0
  if (IntCheck == 0 || SeasonsToSimulate <= 0 || (SeasonsToSimulate > MaxSeasonNumToSimulate)){
    stop("Incompatible SeasonsToFit value")
  }
  
  #Pick out row 4 onwards of ObvsDataTemp.
  #Row 4 corresponds to 2012/13 influenza season, row 5 to 2013/14 influenza season etc
  ObvsData = ObvsDataTemp[(4:(nrow(ObvsDataTemp)-SeasonsUnused)),]

  #--------------------------------------------------------------------------
  # SET UP APMC SCHEME RELATED PARAMETERS
  #--------------------------------------------------------------------------
  
  #Specify target number of samples and fraction kept at each iteration
  #(alpha)
  N_alpha = as.integer(ARGS[4])
  alpha = as.numeric(ARGS[5])
  
  #Set minimal acceptance rate
  MinAcceptRate = as.integer(ARGS[6])
  MaxGen = as.integer(ARGS[7])
  
  #Specify APMC related functions
  s_1 = ARGS[8] #Convert string to Symbol
  s_2 = ARGS[9] #Convert string to Symbol
  s_3 = ARGS[10] #Convert string to Symbol
  s_4 = ARGS[11] #Convert string to Symbol
  s_5 = ARGS[12] #Convert string to Symbol
  
  #Make Symbols callable functions
  PerturbVarFn = get(s_1) #Fn specifying perturbation distribution for newly proposed samples
  PriorFn = get(s_2) #Fn to check whether perturbed samples are within prior bounds
  SummStatFn = get(s_3) #Error measure, compare observed to simulated data
  SampleFromPriorFn = get(s_4) #Sampling from prior distribution
  ModelSimnFn = get(s_5) #Model simulation
  
  #Indicator variable
  # 0 - Sample first particle sets (N in total) from prior distribution
  # 1 - 1 - Read particle sets from file, the retained particles at end
  #       of a previous inference run (N_alpha in total)
  FirstGenFromFileFlag = as.integer(ARGS[13])
  
  #Run function
  RunInferenceAPMC(RunID,ObvsData,SeasonsToSimulate,
                   N_alpha,alpha,MinAcceptRate,MaxGen,PerturbVarFn,PriorFn,SummStatFn,
                   SampleFromPriorFn,ModelSimnFn,FirstGenFromFileFlag)
}

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
SetUpInferenceAPMC(ARGS)
