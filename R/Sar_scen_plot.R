#plotting the  Spawning Stock Biomass of the base OM and alternative climate informed OMs for the CJFAS paper
#the runs were produced using SSMSE and the code CJFAS2026_scenarios_runs_MICERec. R for the climate informed OMs
#and CJFAS2026_scenarios_base.R for the base case autocorrelated recruitment.
#those scripts are in the sardineMSE repository https://github.com/futureseas/SardineMSE/tree/main 
#Author: D. Tommasi
library(dplyr)
library(r4ss)
library(tidyverse) 

#Can skip this step if the files have already been extracted, see below.
#------------------EXTRACTING THE OUTPUT FROM THE SSMSE RUNS-----------------------
#setting the working directory
#setwd("C:/Users/desiree.tommasi/Documents/CAFA/Synthesis")
#when you clone the MSescenarios repository https://github.com/futureseas/MSEscenarios 
#and open the project, the working directory should automatically be MSEscenarios
# check that's the case
getwd()
#specify path of working directory
mdir <- getwd()
niter <- 100
nyrs <- 50
#start of simulation is 2020
years <- 2019:(2019+nyrs)

#create data frame to store SSB output
SSBmat <- data.frame(Year=rep(years, niter),SSB=rep((rep(1,length(years))),niter),itr=rep((rep(1,length(years))),niter))

ct=1
for (i in seq(1,dim(SSBmat)[1],length(years))){
  itr<-ct
  out_dir <- paste0("C:/Users/desiree.tommasi/Documents/CAFA/Synthesis/OM_outputs/No_Climate/",itr, "/constGrowthMidSteepNewSelex_OM_OM/", sep="")
  om_out <- SS_output(out_dir, covar = FALSE)
  #Extract spr series quantities from the SS output
  spr.om <- om_out$sprseries[which(om_out$sprseries$Yr %in% 2019:(2019+nyrs)),]
  SSB <- spr.om$SSB
  SSBmat$SSB[i:(i+length(years)-1)] <- SSB
  SSBmat$itr[i:(i+length(years)-1)] <- itr
  ct <- ct+1
}

#save output to file
write.table(SSBmat, paste0(mdir, "/OM_outputs/Out_No_climate.csv")) 
#------------------------------------------------------------------------------------------------------

#----------------------------SSB PLOTS FOR CJFAS PAPER----------------------------------------------
#read output extracted from the SSMSE simulations
#the simulations projected each OM forward from 2020 to 2068
#each OM was run 100 time with different recruitment process error
#the simulation output is available at the repository https://github.com/futureseas/MSEscenarios
SSBmatN <- read.table(paste0(mdir, "/OM_outputs/Out_No_climate.csv"))
SSBmatN$scenario = "Base OM"
SSBmatR <- read.table(paste0(mdir, "/OM_outputs/Out_R_climate.csv")) 
SSBmatR$scenario = "Recruitment"
SSBmatM <- read.table(paste0(mdir, "/OM_outputs/Out_M_climate.csv")) 
SSBmatM$scenario = "Natural Mortality"
SSBmatG <- read.table(paste0(mdir, "/OM_outputs/Out_G_climate.csv")) 
SSBmatG$scenario = "Growth Hadley"
SSBmatG2 <- read.table(paste0(mdir, "/OM_outputs/Out_Ggfd_climate.csv")) 
SSBmatG2$scenario = "Growth GFDL"
SSBmatMG <- read.table(paste0(mdir, "/OM_outputs/Out_MG_climate.csv")) 
SSBmatMG$scenario = "Natural Mortality \nand Growth Hadley"

#generate summary across the 100 iterations for each scenario
SSBmat = rbind(SSBmatN,SSBmatR,SSBmatM,SSBmatG, SSBmatG2, SSBmatMG)

SSBsmry <- SSBmat %>% group_by(Year,scenario) %>% summarize(
  meanSSB = mean(SSB),
  q5SSB = quantile(SSB, probs = 0.05),
  q95SSB = quantile(SSB, probs = 0.95)
)

#ssb time series plot
ssbplot=ggplot(SSBsmry, aes(x=Year, y=meanSSB,color=scenario, fill=scenario)) + 
  geom_ribbon(aes(x=Year, ymin=q5SSB, ymax=q95SSB),alpha=0.3)+
  geom_line(linewidth=1)+xlab("")+ylab("Spawning Stock Biomass (mt)")+
  theme_bw()+
  theme(axis.title = element_text(size = rel(1.25)),axis.text = element_text(size = rel(1.25)),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = rel(1.25)),
        legend.text = element_text(size=rel(1.20)), legend.title = element_text(size=rel(1.25)))+
  scale_fill_manual(name = "OM Scenario",values = c("Base OM"="steelblue2","Growth GFDL"="olivedrab",
                                "Growth Hadley"="green4", "Natural Mortality"="orangered3",
                               "Natural Mortality \nand Growth Hadley" = "orange",
                               "Recruitment" = "orchid3"))+
  scale_color_manual(values = c("Base OM"="steelblue2","Growth GFDL"="olivedrab",
                               "Growth Hadley"="green4", "Natural Mortality"="orangered3",
                               "Natural Mortality \nand Growth Hadley" = "orange",
                               "Recruitment" = "orchid3"))+
  guides(color = "none")

dev.new(width=10, height=6)
png(paste("ssb_tseries.png",sep=""), width = 10, height = 6, units = 'in', res = 300)
ssbplot
dev.off()


#add a new "Period" factor, this specifies the time slice of 16/17 years over which to plot each boxplot
SSBmat$Period <- "2020-2035"
SSBmat$Period[which(SSBmat$Year%in%c(2036:2052))] <- "2036-2052"
SSBmat$Period[which(SSBmat$Year%in%c(2053:2069))] <- "2053-2069"

#plot SSB distribution by time slice
ssbplot2 <- ggplot()+
  geom_boxplot(data=SSBmat%>%filter(Year>2019), aes(x=Period, y=SSB, fill = scenario), 
               outliers=FALSE)+xlab("")+
  ylab("Spawning Stock BIomass (mt)")+
  theme_bw()+theme(legend.title=element_blank())+
  scale_fill_manual(name = "OM Scenario",
                    values = c("Base OM"="cyan3","Growth GFDL"="olivedrab",
                               "Growth Hadley"="green2", 
                               "Natural Mortality"="orangered3",
                               "Natural Mortality \nand Growth Hadley" = "orange",
                               "Recruitment" = "orchid3"))
dev.new(width=5, height=4)
png("ssbboxplot.png", width = 6, height = 4, units = 'in', res = 300)
ssbplot2
dev.off()
