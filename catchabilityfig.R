#Catchability figure
library(ggplot2)
library(tidyverse)
library(zoo)
library(ggpubr)

#when this project is opened the working directory should automatically be MSEscenarios
# check that's the case
getwd()
#specify path of working directory
wpath <- getwd()

#Read projected mean monthly probability of sardine occurrence for each ESM from 1980 to 2100
#within the CalCOFI area and outside
#The .rds file was computed by Barb Muhling using an SDM trained on sardine larvae
#presence/absence data
CalCOFIlar <- readRDS(paste0(wpath,"/sardMeanLarvaeInOutCalcofiFootprint_1980_2100.rds"))

#change the season to reflect the OM ones
CalCOFIlar$season <- 2
CalCOFIlar$season[which(CalCOFIlar$mo %in% c(7:12))] <- 1
#add a model year (starts in July)
CalCOFIlar$yrm <- CalCOFIlar$yr
CalCOFIlar$yrm[which(CalCOFIlar$season==2)] <- CalCOFIlar$yr[which(CalCOFIlar$season==2)]-1
#choose years > 1980 to have even periods
CalCOFIlar2 <- CalCOFIlar%>%filter(yrm>1980)
#add a time period over which to compute statistics
CalCOFIlar2$TimeBlock <- "1981-2000"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2001:2020))] <- "2001-2020"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2021:2040))] <- "2021-2040"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2041:2060))] <- "2041-2060"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2061:2080))] <- "2061-2080"
CalCOFIlar2$TimeBlock[which(CalCOFIlar2$yrm %in% c(2081:2100))] <- "2081-2100"

#take mean for season 2
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

t_po <- ggplot(po_ma %>% filter(yrm>2000), aes(x=yrm, y=ma,color=ESM,fill=ESM)) + 
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(x=yrm, ymin=p5ma, ymax=p95ma),alpha=0.2)+
  xlab("")+ylab("Mean probability of occurence (/1)")+
  guides(color  = guide_legend(position = "inside"))+
  theme_bw()+coord_cartesian(ylim=c(0,0.6),xlim=c(2001,2100),expand=FALSE)+
  theme(legend.position.inside = c(0.85,0.865),
        legend.background = element_rect(fill = "transparent"))+
  scale_x_continuous(breaks = c(2020,2040,2060,2080))+
  facet_wrap(~area, ncol=1)

#load rds to generate map
pntsToExtract <- readRDS(paste0(wpath,"/calcofiAndNorthPointsToExtract.rds"))

epo.coast <- borders("world", colour = "gray50", fill = "gray80", xlim = c(-140, -110), ylim = c(0, 60))

areaMap <- ggplot(pntsToExtract) + geom_point(aes(x = lon, y = lat, color = type), show.legend = FALSE) +
  scale_color_manual(values = c("dodgerblue", "blue")) + xlab("Longitude") + ylab("Latitude") +
  epo.coast + coord_quickmap(xlim = c(-130, -115), ylim = c(30, 45)) + theme_bw()

figurex <- ggarrange(t_po, areaMap,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

ggsave(paste0(wpath,"/Qfig.tiff"), 
       plot = figurex, width = 3200, height = 1900, units = "px", dpi = 400, compression = "lzw")
