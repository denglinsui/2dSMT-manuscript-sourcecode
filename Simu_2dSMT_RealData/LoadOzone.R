##=== Load Data
library(data.table)
Fn <- list.files("Data/Ozone",full.names = T)

##=== Re-orgnize data from 2000 to 2020
Annual_Summary <- NULL
for(fn in Fn){
  Annual <- as.data.table(read.csv(fn))
  Annual_Tmp <- Annual[Pollutant.Standard=="Ozone 8-hour 2015"&Parameter.Name == "Ozone",
                       .(Avr = mean(Arithmetic.Mean)),
                       by = .(State.Code,County.Code,Site.Num,Parameter.Code,
                              Latitude,Longitude,Datum,Pollutant.Standard,Metric.Used,
                              Year,Units.of.Measure)]
  Annual_Summary <- rbind(Annual_Summary,Annual_Tmp)
}

#--- Preparation (Lon&Lat uniquely determine an annual record.)
#-- Remove locations away from mainland
Annual_Summary <- Annual_Summary[Longitude>=-130&Latitude>=20]
#-- Remove different Datum： "NAD83" "WGS84"
Annual_Summary <- Annual_Summary[Datum=="WGS84"]
#-- Transform State.Code
Annual_Summary$State.Code <- as.numeric(Annual_Summary[,State.Code])


#=== Remove Data with different site
#-- Find Lon&Lat should be removed
Site_Length <- Annual_Summary[,.(Obs_site=length(unique(Site.Num))),
                                by = .(Latitude,Longitude)]
Site_Length <- Site_Length[Obs_site==1]
#-- Just preserve the data with single site
Annual_Summary <- Annual_Summary[paste(Longitude,Latitude)%in%paste(Site_Length$Longitude,Site_Length$Latitude)]


##=== Maintain the observations that are continuously observed from 2000 to 2021
#-- Site.Num might change, but it doesn't matter the final result
Annual_Length <- Annual_Summary[Year>=2010,.N,
                                by = .(County.Code,Site.Num,Parameter.Code,
                                       Latitude,Longitude,Datum,Metric.Used)]
Annual_Length <- Annual_Length[N==12]
#-- Just preserve the data with full record
Annual_Summary <- Annual_Summary[Year>=2010&paste(Longitude,Latitude)%in%paste(Annual_Length$Longitude,Annual_Length$Latitude)]

##=== Summary Data
#-- How many site?
nrow(Annual_Summary)/12

#-- Plot
library(ggplot2)
library(ggmap)
state <- map_data("state")
ggplot() +
  geom_map(
    data = state, map = state,
    aes(map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )+
  theme_void() +
  geom_point(aes(x=Longitude,y=Latitude), colour = "black",data=Annual_Length)
ggsave("Fig/geom_map.eps",width = 4,height=3)

#-- Mercator projection
library(mapproj)
test <- mapproject(Annual_Length$Longitude,Annual_Length$Latitude,projection="mercator")
qplot(test$x,test$y)
ggsave("Fig/mer.eps",width = 4,height=3)

#-- Save Data
#save(Annual_Summary,file="Data/Ozone.RData")
