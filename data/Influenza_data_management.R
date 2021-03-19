# This code is to prepare for US calibration
# code generated on: 2021. 03.18

library(openxlsx)
library(tidyverse)
library(reshape2)
setwd("/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/Calibration/ANN_BUGS_Cleaned/Influenza-github/data")
## 1. Vaccine uptake
cumulvaxbymonth <- read.xlsx("USVaccUptakeByMonth.xlsx")
cumulvaxbymonth[,2:4] = cumulvaxbymonth[,2:4]/100
incvaxbymonth = cumulvaxbymonth
dailyvaxrate = cumulvaxbymonth
for (i in 2:ncol(cumulvaxbymonth)){
  inc = cumulvaxbymonth[,i] - c(0,cumulvaxbymonth[1:(nrow(cumulvaxbymonth)-1),i])
  incvaxbymonth[,i] = inc
  dailyvaxrate[,i] = -log(1-incvaxbymonth[,i])/30.5
}
dailyvaxrate <- rbind(dailyvaxrate, c('June',rep(0,ncol(dailyvaxrate)-1)))
dailyvaxrate <- rbind(dailyvaxrate, c('July',rep(0,ncol(dailyvaxrate)-1)))

# compare with UK data
path = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/SeasonalFluImmunityPropagation-master/"
ukvaxuptake <- read.xlsx(paste0(path,"Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx"),sheet="All popn.", colNames = F, rows = c(4:12), cols = c(3:368))
ukvaxuptake_p <- read.xlsx(paste0(path,"Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx"),sheet="2009_2010",colNames = F, rows =5, cols = c(3:368))

usvaxuptake <- matrix(rep(NA, (ncol(dailyvaxrate)-1)*365),nrow=365,ncol=ncol(dailyvaxrate)-1)
for (i in 2:ncol(dailyvaxrate)){
  print(i)
  start_idx=1
  end_idx = 0
  for (j in 6:17){
    if (j >12){
      j = j-12
    }
    print(j)
    if (dailyvaxrate[j,1] %in% c("Jan","Mar","May","July","August","Oct","Dec")){
      rep_num = 31
    } else if(dailyvaxrate[j,1] == 'Feb'){
      rep_num = 28
    } else{
      rep_num = 30
    }
    end_idx = end_idx + rep_num
    usvaxuptake[start_idx:end_idx,(i-1)] = rep(dailyvaxrate[j,i],rep_num)
    start_idx = end_idx + 1
  }
}

# from 2009/10 to 2017/18
usvaxuptake_s <- as.data.frame(usvaxuptake[,-c(2,11,12)])
colnames(usvaxuptake_s) <- c("2009/10","2010/11","2011/12","2012/13","2013/14","2014/15","2015/16","2016/17","2017/18")
usvaxuptake_s$day = seq(1,365)
usvaxuptake_p <- as.data.frame(usvaxuptake[,2])

ukvaxuptake <- rbind(ukvaxuptake, seq(1,365))
ukvaxuptake_s <- as.data.frame(t(ukvaxuptake))
colnames(ukvaxuptake_s) <- c("2009/10","2010/11","2011/12","2012/13","2013/14","2014/15","2015/16","2016/17","2017/18",'day')
ukvaxuptake_p <- t(ukvaxuptake_p)

usvaxuptake_s_t <- melt(usvaxuptake_s, id.vars = 'day')
usvaxuptake_s_t$country <- 'us'
usvaxuptake_s_t$value = as.numeric(usvaxuptake_s_t$value)
ukvaxuptake_s_t <- melt(ukvaxuptake_s, id.vars = 'day')
ukvaxuptake_s_t$country <- 'uk'
allvaxuptake_s <- rbind(usvaxuptake_s_t, ukvaxuptake_s_t)

ggplot( data = allvaxuptake_s)+
  geom_line(aes(x=day, y=value, color=country))+
  facet_wrap(~variable)


## 2. Vaccine efficacy
ukvaxeff = read.xlsx(paste0(path,"Data/VaccEfficacy/VaccEfficacy_AllPopn.xlsx"),sheet="MidPoint", colNames = F, rows = c(3:11), cols = c(3:6))
usvaxeff = read.xlsx("VaccEfficacy_AllPopn_US.xlsx",colNames = T)
usvaxeff = usvaxeff[-c(1),-c(1)]
colnames(ukvaxeff)<-c("H1N1","H3N2","B/Y","B/V")
colnames(usvaxeff)<-c("H1N1","H3N2","B/Y","B/V")
ukvaxeff$season <- c("2009/10","2010/11","2011/12","2012/13","2013/14","2014/15","2015/16","2016/17","2017/18")
usvaxeff$season <- c("2009/10","2010/11","2011/12","2012/13","2013/14","2014/15","2015/16","2016/17","2017/18")
ukvaxeff$country <- 'uk'
usvaxeff$country <- 'us'
allvaxeff <- rbind(melt(ukvaxeff,id.vars = c('season','country')),melt(usvaxeff,id.vars = c('season','country')))

ggplot(data = allvaxeff, aes(x=variable, y=value, fill=country))+
  geom_bar(position = "dodge",stat = "identity", width = 0.6)+
  facet_wrap(~season)
  