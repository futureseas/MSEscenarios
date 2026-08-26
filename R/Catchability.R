#Code to develop OM catchability scenarios from SDM output
#Author:D. tommasi

library(ggplot2)
library(tidyverse)
library(zoo)

#when you clone the MSescenarios repository https://github.com/futureseas/MSEscenarios 
#and open the project, the working directory should automatically be MSEscenarios
# check that's the case
getwd()
#specify path of working directory
wpath <- getwd()

#Read projected mean monthly probability of sardine occurrence for each ESM from 1980 to 2100
#within the CalCOFI area and outside
#The .rds file was computed by Barb Muhling using an SDM trained on sardine larvae
#The SDMs uses presence/absence data

CalCOFIlar <- readRDS(paste0(wpath,"/dat/sardMeanLarvaeInOutCalcofiFootprint_1980_2100.rds"))


#change the season to reflect the OM ones
CalCOFIlar$season <- 2
CalCOFIlar$season[which(CalCOFIlar$mo %in% c(7:12))] <- 1
#add a model year (starts in July)
CalCOFIlar$yrm <- CalCOFIlar$yr
CalCOFIlar$yrm[which(CalCOFIlar$season==2)] <- CalCOFIlar$yr[which(CalCOFIlar$season==2)]-1
#choose years > 1980 to have even time slices
CalCOFIlar2 <- CalCOFIlar%>%filter(yrm>1980)
#add a time period factor over which to compute statistics
CalCOFIlar2$TimeBlock <- "1981-2000"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2001:2020))] <- "2001-2020"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2021:2040))] <- "2021-2040"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2041:2060))] <- "2041-2060"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2061:2080))] <- "2061-2080"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2081:2100))] <- "2081-2100"

#take mean for season 2 since this is the spawning season
#Compute mean seasonal prob of occurrence for each area type for months of Jan-Jun (season 2)
larmat.2 <- CalCOFIlar2 %>% filter(season==2) %>% group_by(esm,type,yrm,season,TimeBlock)%>%summarize(
  avail=mean(pred),
  p5 = quantile(pred,probs = 0.05,na.rm = TRUE),
  p95 = quantile(pred, probs = 0.95, na.rm = TRUE))

#plot a 5-yr smoother of the probability of occurrence by region
#generate moving average
po_gfdl<-larmat.2 %>% filter (type=="calcofi"&esm=="gfdl")
po_ipsl<-larmat.2 %>% filter (type=="calcofi"&esm=="ipsl")
po_hadl<-larmat.2 %>% filter (type=="calcofi"&esm=="hadl")
po_gfdl$ma <- rollmean(po_gfdl$avail, k= 5, align = "right", fill = NA)
po_ipsl$ma <- rollmean(po_ipsl$avail, k= 5, align = "right", fill = NA)
po_hadl$ma <- rollmean(po_hadl$avail, k= 5, align = "right", fill = NA)
po_gfdl$p5ma <- rollmean(po_gfdl$p5, k= 5, align = "right", fill = NA)
po_ipsl$p5ma <- rollmean(po_ipsl$p5, k= 5, align = "right", fill = NA)
po_hadl$p5ma <- rollmean(po_hadl$p5, k= 5, align = "right", fill = NA)
po_gfdl$p95ma <- rollmean(po_gfdl$p95, k= 5, align = "right", fill = NA)
po_ipsl$p95ma <- rollmean(po_ipsl$p95, k= 5, align = "right", fill = NA)
po_hadl$p95ma <- rollmean(po_hadl$p95, k= 5, align = "right", fill = NA)
po_gfdl$area <- "CalCOFI Survey Region"
po_ipsl$area <- "CalCOFI Survey Region"
po_hadl$area <- "CalCOFI Survey Region"


#plot a 5-yr smoother of the probability of occurrence by region
#generate moving average
pon_gfdl<-larmat.2 %>% filter (type=="north"&esm=="gfdl")
pon_ipsl<-larmat.2 %>% filter (type=="north"&esm=="ipsl")
pon_hadl<-larmat.2 %>% filter (type=="north"&esm=="hadl")
pon_gfdl$ma <- rollmean(pon_gfdl$avail, k= 5, align = "right", fill = NA)
pon_ipsl$ma <- rollmean(pon_ipsl$avail, k= 5, align = "right", fill = NA)
pon_hadl$ma <- rollmean(pon_hadl$avail, k= 5, align = "right", fill = NA)
pon_gfdl$p5ma <- rollmean(pon_gfdl$p5, k= 5, align = "right", fill = NA)
pon_ipsl$p5ma <- rollmean(pon_ipsl$p5, k= 5, align = "right", fill = NA)
pon_hadl$p5ma <- rollmean(pon_hadl$p5, k= 5, align = "right", fill = NA)
pon_gfdl$p95ma <- rollmean(pon_gfdl$p95, k= 5, align = "right", fill = NA)
pon_ipsl$p95ma <- rollmean(pon_ipsl$p95, k= 5, align = "right", fill = NA)
pon_hadl$p95ma <- rollmean(pon_hadl$p95, k= 5, align = "right", fill = NA)
pon_gfdl$area <- "North of CalCOFI Survey Region"
pon_ipsl$area <- "North of CalCOFI Survey Region"
pon_hadl$area <- "North of CalCOFI Survey Region"

#combine back into data frame for plotting
po_ma = rbind(po_gfdl,po_hadl,po_ipsl,pon_gfdl,pon_hadl,pon_ipsl)
names(po_ma)[1]="ESM"

t_po=ggplot(po_ma %>% filter(yrm>2000), aes(x=yrm, y=ma,color=ESM,fill=ESM)) + 
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(x=yrm, ymin=p5ma, ymax=p95ma),alpha=0.2)+
  xlab("")+ylab("Mean probability of occurence (/1)")+
  guides(color  = guide_legend(position = "inside"))+
  theme_bw()+coord_cartesian(ylim=c(0,0.6),xlim=c(2001,2100),expand=FALSE)+
  theme(legend.position.inside = c(0.85,0.865),
        legend.background = element_rect(fill = "transparent"))+
  scale_x_continuous(breaks = c(2020,2040,2060,2080))+
  facet_wrap(~area, ncol=1)

#This plot is for Fig. 8
dev.new(width=3, height=5)
png(paste("Qtimeseries.png",sep=""), width = 3, height = 5, units = 'in', res = 300)
t_po
dev.off()


#Compute mean monthly prob of occurrence for each area type for months of Jan-Jun (season 2)
larmat.mo <- CalCOFIlar2 %>% filter(season==2) %>% group_by(esm,type,mo,yrm,season,TimeBlock)%>%summarize(
  avail=mean(pred))
#compute the total prob of occurrence for the entire region
larmat_tot.mo <-larmat.mo %>% group_by(esm,mo,yrm,season,TimeBlock)%>%summarize(availT = sum(avail))
#extract prob of occurrence in calCOFI region and add it to the total dataframe
larmat_tot.mo$availC <- (larmat.mo%>%filter(type=="calcofi"))$avail
#calculate fraction in CalCOFI region by esm
Qmat <- larmat_tot.mo %>% group_by(esm,mo,yrm,season,TimeBlock)%>%summarize(fracC = availC/availT)

#plot Q output by 20-yr time block as projections describe mean climate over 20-30 yrs periods
#use as historical the 2001-2020 time block since OM starts in 2001
Qplot <- ggplot(Qmat %>% filter(yrm>2000), aes(x=TimeBlock, y=fracC, fill = esm))+
  geom_boxplot()+xlab("")+ylab("Catchability Index\n(Fraction of total probability of occurence in survey area)")+
  stat_summary(fun = mean, geom="point", shape = 8, size=2, color="black",
               position=position_dodge2(width = 0.7, preserve="single"))+
  theme_bw()+theme(legend.title=element_blank(), 
                   legend.text = element_text(size=12),
                   axis.text = element_text(size=12))
dev.new(width=10, height=8)
png("Qplot_byesm.png", width = 10, height = 8, units = 'in', res = 300)
Qplot
dev.off()

#-------------------------------Fig. 9---------------------------------------------------------
#plot Q across all ESMs
#add label for fill color
Qmat$Period <-"Historical"
Qmat$Period[Qmat$yrm>2020] <-"Climate"
Qplot2 <- ggplot(Qmat%>% filter(yrm>2000), aes(x=TimeBlock, y=fracC,fill = Period))+
  geom_boxplot()+xlab("")+ylab("Catchability Index\n(Fraction of total probability\nof occurence in survey area)")+
  #stat_summary(fun = mean, geom="point", shape = 8, size=2, color="black",
               #position=position_dodge2(width = 0.7, preserve="single"))+
  theme_bw()+theme(legend.title=element_blank(), 
                   legend.text = element_text(size=12),
                   axis.text = element_text(size=12),
                   axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_fill_manual(values = c("Historical"="cyan3","Climate"="brown2"),
                     labels = c("Historical"="Historical","Climate"="Future\nClimate"))
dev.new(width=5, height=4)
png("Qplot2.png", width = 5, height = 4, units = 'in', res = 300)
Qplot2
dev.off()

#-----------------------------------------------------------------------------------------------------

#--------------------Output for manuscript tables-----------------------------------------------------
#Calculate mean and 25th and 75th quantiles by time block across ESMs
Qavg <- Qmat %>% group_by(TimeBlock) %>%summarize(meanQ = mean(fracC),
                                                        q25=quantile(fracC, probs = 0.25),
                                                        q75=quantile(fracC, probs = 0.75),
                                                        min=min(fracC),
                                                        max=max(fracC),
                                                        sd=sd(fracC))

#compute the proportional change in Q by projection time block vs. historical (2001-2020)
histQavg <- Qavg %>% filter(TimeBlock=="2001-2020")
projQavg <- Qavg %>% filter(!TimeBlock%in% c("1981-2000","2001-2020"))
projQavg$meanQh <-histQavg$meanQ
Qchange <- projQavg %>% group_by(TimeBlock) %>%summarize(changeQ = 100*(meanQ-meanQh)/meanQh)

#calculate %change by esm
#Calculate mean and 25th and 75th quantiles by time block across ESMs
Qavg.esm <- Qmat %>% group_by(TimeBlock,esm) %>%summarize(meanQ = mean(fracC),
                                                  q25=quantile(fracC, probs = 0.25),
                                                  q75=quantile(fracC, probs = 0.75),
                                                  min=min(fracC),
                                                  max=max(fracC),
                                                  sd=sd(fracC))

#compute the proportional change in Q by projection time block vs. historical (2001-2020)
histQavg.esm <- Qavg.esm %>% filter(TimeBlock=="2001-2020")
projQavg.esm <- Qavg.esm %>% filter(!TimeBlock%in% c("1981-2000","2001-2020"))
projQavg.esm$meanQh <-rep(histQavg.esm$meanQ,4)
Qchange.esm <- projQavg.esm %>% group_by(TimeBlock,esm) %>%summarize(changeQ = 100*(meanQ-meanQh)/meanQh)

#Apply the proportional change to the assumed baseline Q in the OM
Qchange$newQ <- 0.001+0.001*Qchange$changeQ/100
Qchange.esm$newQ <- 0.001+0.001*Qchange.esm$changeQ/100


