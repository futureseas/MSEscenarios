#Code to develop a climate-informed natural mortality scenario for Pacific sardine
#from output of projections of consumption rates from all predation on sardine from the 
#Atlantis ecosystem model parametrized for the California Current
#see Liu et al. 2025 for details on the Atlantis model
#to generate this projections Atlantis was forced with 2013-2100 regionally 
#downscaled projections from three Earth System Models (GFDL, IPSL, and Hadley)
#and SDMs projected from 2013-2100 using the same downscaled ESMs
#The Base run is forced by looping 2013 oceanography and SDMs from 2013 
#See Liu et al. 2025 for details on the projections and base runs

library(tidyverse)

#when this project is opened the working directory should automatically be MSEscenarios
# check that's the case
getwd()
#specify path of working directory
wpath <- getwd()

#load the Atlantis output provided by Isaac Kaplan at NWFSC
#these are projected consumption rates per year (predation M2 for sardine) and
#sardine biomass in mt for the three ESMs and base run
#and IPSL MICE recdevs computed as for Wildermuth et al. 2023 
#note the first 7 years are not included as they are part of the spin-up period
Mdat <- read.csv(paste0(wpath,"/AtlantisM.csv"))

#plot M output by 20-yr time block as projections describe mean climate over 20-30 yrs periods
Mplot <- ggplot(Mdat, aes(x=TimeBlock, y=M, fill = ESM))+
  geom_boxplot()+xlab("")+ylab("Predation Mortality (per year)")+
  stat_summary(fun = mean, geom="point", shape = 8, size=2, color="black",
               position=position_dodge2(width = 0.7, preserve="single"))+
    theme_bw()+theme(legend.title=element_blank(), 
                     legend.text = element_text(size=12),
                     axis.text = element_text(size=12))
dev.new(width=10, height=8)
png("Mplot_byesm.png", width = 10, height = 8, units = 'in', res = 300)
Mplot
dev.off()

#plot M distribution across all ESMs
Mplot2 <- ggplot(Mdat, aes(x=TimeBlock, y=M, fill = Type))+
  geom_boxplot()+xlab("")+ylab("Predation Mortality (per year)")+
  stat_summary(fun = mean, geom="point", shape = 8, size=2, color="black",
               position=position_dodge2(width = 0.7, preserve="single"))+
  theme_bw()+theme(legend.title=element_blank(), 
                   legend.text = element_text(size=12),
                   axis.text = element_text(size=12))+
  scale_fill_discrete(labels=c("Non-climate","Climate"))
dev.new(width=10, height=8)
png("Mplot.png", width = 10, height = 8, units = 'in', res = 300)
Mplot2
dev.off()

#calculate statistics over entire projection period
Mavg1 <- Mdat %>% group_by(Type) %>%summarize(meanM = mean(M),
                                                        q25=quantile(M, probs = 0.25),
                                                        q75=quantile(M, probs = 0.75),
                                                        min=min(M),
                                                        max=max(M),
                                                        sd=sd(M))


#Calculate mean and 25th and 75th quantiles by time block across ESMs
Mavg <- Mdat %>% group_by(TimeBlock, Type) %>%summarize(meanM = mean(M),
                                                        q25=quantile(M, probs = 0.25),
                                                        q75=quantile(M, probs = 0.75),
                                                        min=min(M),
                                                        max=max(M),
                                                        sd=sd(M))

#compute the proportional change in M by time block
Mchange <- Mavg %>% group_by(TimeBlock) %>%summarize(changeM = 100*(meanM[Type=="ESM"]-meanM[Type=="Base"])/meanM[Type=="Base"])
#Apply the proportional change to the assumed M in the OM
Mchange$newM <- 0.585+0.585*Mchange$changeM/100

#calculate change by ESM, also for 25th and 75th quantiles
Mavgesm <- Mdat %>% group_by(TimeBlock, Type, ESM) %>%summarize(meanM = mean(M),
                                                                q25=quantile(M, probs = 0.25),
                                                                q75=quantile(M, probs = 0.75),
                                                                sd=sd(M))
#compute the proportional change in M by time block
Mchangegfdl <- Mavgesm %>% filter(ESM %in% c("Base", "GFDL"))%>% group_by(TimeBlock) %>%summarize(
  changeM = 100*(meanM[Type=="ESM"]-meanM[Type=="Base"])/meanM[Type=="Base"])
Mchangegfdl$ESM <-"GFDL"
MchangeHad <- Mavgesm %>% filter(ESM %in% c("Base", "Hadley"))%>% group_by(TimeBlock) %>%summarize(
  changeM = 100*(meanM[Type=="ESM"]-meanM[Type=="Base"])/meanM[Type=="Base"])
MchangeHad$ESM <- "Hadley"
MchangeIpsl <- Mavgesm %>% filter(ESM %in% c("Base", "IPSL"))%>% group_by(TimeBlock) %>%summarize(
  changeM = 100*(meanM[Type=="ESM"]-meanM[Type=="Base"])/meanM[Type=="Base"])
MchangeIpsl$ESM <- "IPSL"

#join in one data frame
Mchangesm = rbind(Mchangegfdl,MchangeHad,MchangeIpsl)

#Apply the proportional change to the assumed M in the OM
Mchangesm$newM <- 0.585+0.585*Mchangesm$changeM/100

#Calculate statistics of sardine biomass and M by Type and Yr to see if relationship between
#temporal variability in M and sardine biomass
Mavgy <- Mdat %>% group_by(Year, Type) %>%summarize(meanM = mean(M),
                                                    meanBio = mean(SarBio),
                                                    meanR = mean(MICER),
                                                        q5=quantile(M, probs = 0.05),
                                                        q95=quantile(M, probs = 0.95),
                                                    bq5=quantile(SarBio, probs = 0.05),
                                                    bq95=quantile(SarBio, probs = 0.95))
#Scale sardine biomass before plotting so that it is comparable to M
Mavgy$SarBio2 <- (Mavgy$meanBio-min(Mavgy$meanBio))/(max(Mavgy$meanBio)-min(Mavgy$meanBio))

Mbio <- ggplot(Mavgy, aes(x=Year)) +
  geom_line(aes(y=SarBio2),color="blue")+geom_line(aes(y=meanM),color="red")+
  xlab("")+
  ylab("Scaled sardine biomass\nand predation mortality per year")+
  theme_bw()+facet_wrap(~Type)+
  theme(axis.title = element_text(size = rel(1.5)),axis.text = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = rel(1.5)))
dev.new(width=10, height=8)
png("Mbio.png", width = 10, height = 8, units = 'in', res = 300)
Mbio
dev.off()

#Create a time-varying mortality using the output from the Atlantis non-climate run
#the Atlantis simulation showed that predation mortality (M2) is density-dependent
#being highest at lower sardine abundance
#Since now sardine is at a low biomass, we assume the current M from the 
#assessment, 0.53, includes the highest M2 simulated by Atlantis and calculate what the
#constant, non-predation mortality (M1) is.
#The M2 will fluctuate between high and low in 20-years periods

#set M to the one from the OM
M <- 0.585
#extract the highest M2
M2high <- max(Mavg$meanM)
M1 <- M - M2high
#extract the lowest M2
M2low <- min(Mavg$meanM)
Mlow <- M2low+M1

#generate time series of fluctuating M
#since adding to OM that started in 2019, use 5 additional years, then switch to 
#low for 20, high for 20, low for 5 for a total of 50 years
Mtv = c(rep(M,5),rep(Mlow,20),rep(M,20),rep(Mlow,5))
Mtvdat = data.frame(M=Mtv, Yr = 2020:2069)

Mtvaryp <- ggplot(Mtvdat, aes(x=Yr, y=M)) +
  geom_line(color="black", linewidth=1.5)+
  xlab("")+
  ylab("Natural mortality per year")+
  theme_bw()+ylim(0.35,0.65)+
  theme(axis.text.x = element_text(angle = 90))
dev.new(width=6, height=3)
png("Mtvary.png", width = 6, height = 3, units = 'in', res = 300)
Mtvaryp
dev.off()
