###################################################################################################
# Calculate the "footprints" of different fishery-independent surveys: CPS, Calcofi, maybe RREAS
###################################################################################################

library(ggplot2)
library(lubridate)
library(dplyr)
library(sf)
library(sfheaders)
library(ncdf4)
library(stringr)

#get working dirctory
fdir <- getwd()
# Load Calcofi stns and distance from land grid
cal <- read.csv(paste0(fdir,"/erdCalCOFIstns_89dc_925a_54a0.csv")) # Just stns, not catch

# Load pre-calculated distance to land for ROMS grid
ll <- readRDS(paste0(fdir,"/distROMsplusAtlantis_2060N_150110W.rds"))

# Mapping
epo.coast <- borders("world", colour = "gray50", fill = "gray80", xlim = c(-140, -110), ylim = c(0, 60))

# Start with CalCOFI
# Clean up a bit
cal <- cal[2:nrow(cal),]
cal$lon <- as.numeric(as.character(cal$longitude))
cal$lat <- as.numeric(as.character(cal$latitude))
# Recent ERDDAP data upload seems to have messed up date field
# Just get year from cruise no. for now
cal$yr <- as.numeric(str_sub(cal$cruise, start = 1, end = 4))

#Just the last 20 years
cal <- subset(cal, yr >= 2003) # 9698 stations occupied
# Core stations
stns <- aggregate(cruise ~ line + station, cal, FUN = length) # 3267 line/stn combos
colnames(stns)[3] <- "timesSampled"
# Or maybe just by location?
cal$lonrd <- round(cal$lon, 1)
cal$latrd <- round(cal$lat, 1)
locns <- aggregate(cruise ~ lonrd + latrd, cal, FUN = length) # 2263
colnames(locns)[3] <- "timesSampled"
# Just include locations sampled 10+ times
common <- subset(locns, timesSampled >= 10) # 119 locations
calCommon <- left_join(common, cal, by = c("lonrd", "latrd")) # 5802 separate sampled stations
# Plot up sampling effort. Ok, that picks out the usual grid pretty well
ggplot(common) + geom_point(aes(x = lonrd, y = latrd, size = timesSampled)) +
  epo.coast + coord_quickmap(xlim = c(-135, -115), ylim = c(30, 55)) + theme_bw()

# Just a regular hull should capture this ok
dat <- as.data.frame(cbind(common$lonrd, common$latrd))
# Let's make a polygon for mapping with e.g. a 0.5 degree buffer around the actual points
dat1 <- as.data.frame(cbind(dat$lon + 0.5, dat$lat + 0.5))
dat2 <- as.data.frame(cbind(dat$lon - 0.5, dat$lat - 0.5))
dat3 <- as.data.frame(cbind(dat$lon + 0.5, dat$lat - 0.5))
dat4 <- as.data.frame(cbind(dat$lon - 0.5, dat$lat + 0.5))
dat1deg <- rbind(dat, dat1, dat2, dat3, dat4)
ch <- chull(dat1deg)
coords <- dat1deg[c(ch, ch[1]), ]  # closed polygon
colnames(coords) <- c("lon", "lat")
# Seems to work
ggplot(common) + geom_point(aes(x = lonrd, y = latrd, size = timesSampled)) +
  geom_path(data = coords, aes(x = lon, y = lat), color = "red") +
  epo.coast + coord_quickmap(xlim = c(-135, -115), ylim = c(30, 55)) + theme_bw()

# Get ROMS points inside polygon
calPoly <- sfheaders::sf_polygon(obj = coords, x = "lon", y = "lat")
# Kludge it up
calSample <- st_sample(calPoly, 100000)
buffDF <- do.call(rbind.data.frame, calSample)
colnames(buffDF) <- c("lon", "lat")
buffDF$lonrd <- round(buffDF$lon, 1)
buffDF$latrd <- round(buffDF$lat, 1)
pnts <- aggregate(lon ~ lonrd + latrd, buffDF, FUN = length)
pnts$lon <- NULL
colnames(pnts) <- c("lon", "lat")

######################################################################################################
# Also define some points north of this grid, so can show if/how habitat moves into that area
# Maybe the same distance from shore as the northernmost/offshore Calcofi line (line 60 stn 90), 
# but just extended north 
# # stn90 <- subset(calCommon, line == 60 & station == 90) #  -125.77 W, 36.61 N
# source("./albtag/calculateDistanceToLandFn.R")
# stn90 <- data.frame("lon" = -125.77, "lat" = 36.61)
# extractDistLand <- calcDistLand(myData = stn90) # 287920.6 m 
# Subset distance from land grid to just north of 37, and within 287921 m
ll <- subset(ll, lat >= 36.5 & lat <= 45 & lon <= -120 & distLand <= 287921)
# Round to 2dp
ll$lon <- round(ll$lon, 2)
ll$lat <- round(ll$lat, 2)
# Just 0.1x0.1 grid
ll$lon01 <- (ll$lon * 10) %% 1
ll$lat01 <- (ll$lat * 10) %% 1
ll <- subset(ll, lon01 == 0 & lat01 == 0)
ll$lon01 <- ll$lat01 <- NULL
# Drop overlap with the actual CalCOFI grid at ~ 37N
ll <- anti_join(ll, pnts, by = c("lon", "lat"))
# plot(ll$lon, ll$lat); points(pnts$lon, pnts$lat, col = "red")
# Join points
pnts$type <- "calcofi"
ll$type <- "north"
ll$distLand <- NULL
pntsToExtract <- rbind(pnts, ll)
# Drop the few south of 30N, as no ROMS there
pntsToExtract <- subset(pntsToExtract, lat >= 30)

# Make and save a nicer map
areaMap <- ggplot(pntsToExtract) + geom_point(aes(x = lon, y = lat, color = type), show.legend = FALSE) +
  scale_color_manual(values = c("dodgerblue", "blue")) + xlab("Longitude") + ylab("Latitude") +
  epo.coast + coord_quickmap(xlim = c(-130, -115), ylim = c(30, 45)) + theme_bw()
areaMap
ggsave(paste0(fdir,"/mapInsideOutsideCalcofi.tiff"), 
       plot = areaMap, width = 1400, height = 1600, units = "px", dpi = 400, compression = "lzw")
# Save points
saveRDS(pntsToExtract, paste0(fdir,"/calcofiAndNorthPointsToExtract.rds"))