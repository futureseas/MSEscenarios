#Comparison of recruitment deviates around DynaMICE Stock Recruitment Curve fit from
#historical and projection period
#code modified from that of R. Wildermuth at https://github.com/futureseas/SardineMSE/blob/main/R/ExploreRecruitIndices.R
#that was used in Wildermuth et al. 2024
#author: D. Tommasi

library(tidyverse)
library(r4ss)

#when you clone the MSescenarios repository https://github.com/futureseas/MSEscenarios 
#and open the project, the working directory should automatically be MSEscenarios
#set the project directory as the main directory
dir <-getwd()

#read in the recruitment and ssb estimated from the DynaMICE. Note this is the ensemble mean
#across all nine DynaMICE ensemble members when forced by the historical reanalysis (RA),
#by each each ESM and averaged across the three ESMs.
ssbrecsMICE <- read.csv(paste0(dir, "/dat/recdevMICE2100.csv"))

#Fit the SR to the historical data only
datMICE <- ssbrecsMICE %>% filter(GCM == "RA") %>% select(ensembleSSB, ensembleRec) %>%
  rename("ssb" = "ensembleSSB", 
         "rec" = "ensembleRec") %>%
  filter(complete.cases(.))

#use Beverton-Holt function
BH_lin <- function(pars,sb) {
  a <- pars[1]
  b <- pars[2]
  recl <- log(a*sb/(b+sb))
  return(recl)
}

#negative log likelihood
nll <- function(pars, data) {
  # Values predicted by the model
  Rpred <- BH_lin(pars, data$ssb)
  Robs <- log(data$rec)
  # Negative log-likelihood 
  -sum(dnorm(x = Robs, mean = Rpred, sd = 1.25, log = TRUE))
}

#find parameters that minimize the log likelihood
fitMICE2 <- optim(par = c(1000,1e10), fn = nll, method = "L-BFGS-B",
                  lower = rep(2, 1), data = datMICE)

#plot the SR and deviations
plotdat <- ssbrecsMICE %>% filter(GCM %in% c("RA", "gcmMEAN")) %>%
  mutate(exLogRec = BH_lin(fitMICE2$par, ensembleSSB)) 
srplot <- ggplot(plotdat) +
  geom_point(aes(x = ensembleSSB, y = log(ensembleRec), color = GCM)) +
  geom_line(aes(x = ensembleSSB, y = exLogRec), linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = c("brown2", "cyan3"), 
                     labels = c("Future\nClimate", "Historical"), name = "") +
  labs(title = "Stock Recruitment Curve Fit to DynaMICE", y = "Recruits (log)", x = "SSB")+
  theme(axis.title = element_text(size = rel(1.25)),axis.text = element_text(size = rel(1.25)),
        strip.text = element_text(size = rel(1.25)),
        legend.text = element_text(size=rel(1.20)), legend.title = element_text(size=rel(1.25)),
        legend.position = "inside",legend.position.inside = c(0.75,0.25))


#load the recruitment deviations from the DynaMICE used for the future simulation
#generated with code "Synthesis_scenarios_runs_MICERec.R"
frdevs <- read.csv(paste0(dir,"/dat/FutureMICERecdevs.csv"))

#calculate the median and quantiles across the 100 iterations, 
#Note all scenarios have the same recruitment deviations, just pick one
devsmry <-frdevs %>% filter(scen=="R_climate") %>% group_by(yr) %>% summarise(
  mrecdev = median(value),
  recdev5 = quantile(value, probs= 0.05),
  recdev95= quantile(value, probs = 0.95)
)

#load historical OM deviations
hrdevs <- read.csv(paste0(dir,"/dat/HistOMRecdevs.csv"))
rdplot=ggplot(devsmry, aes(x=yr, y=mrecdev)) + 
  geom_ribbon(aes(x=yr, ymin=recdev5, ymax=recdev95),alpha=0.3)+
  geom_line(linewidth=1,color="black")+xlab("")+ylab("Recruitment Deviations")+
  theme_bw()+xlim(1994,2068)+labs(title = "OM Recruitment Deviations")+
  geom_line(data=hrdevs, aes(x=Yr,y=dev),linewidth=1,color="black")+
  theme(axis.title = element_text(size = rel(1.25)),axis.text = element_text(size = rel(1.25)),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = rel(1.25)),
        legend.text = element_text(size=rel(1.20)), legend.title = element_text(size=rel(1.25)))

#-------------------------------Figure 4--------------------------------------------------
library(gridExtra)
combo.rec <- grid.arrange(arrangeGrob(srplot,rdplot,ncol=2))
ggsave("Rplot_combo.png",plot = combo.rec, width = 8, height = 4)
