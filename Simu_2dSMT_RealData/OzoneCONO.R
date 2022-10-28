## Extract information about (CO, NO, NO2)
##=== Load Data
library(data.table)
library(ggplot2)
library(ggmap)
library(cowplot)
state <- map_data("state")

Fn <- list.files("Data/Ozone",full.names = T)

##=== Re-orgnize data from 2000 to 2020
Annual_CO_Summary <- NULL
Annual_NO2_Summary <- NULL
for(fn in Fn){
  Annual <- as.data.table(read.csv(fn))
  Annual_CO_Tmp <- Annual[Pollutant.Standard=="CO 8-hour 1971",
                          .(Avr = mean(Arithmetic.Mean)),
                          by = .(State.Code,County.Code,Site.Num,Parameter.Code,
                                 Latitude,Longitude,Datum,Pollutant.Standard,Metric.Used,
                                 Year,Units.of.Measure)]
  Annual_CO_Summary <- rbind(Annual_CO_Summary,Annual_CO_Tmp)
  
  Annual_NO2_Tmp <- Annual[Pollutant.Standard=="NO2 1-hour 2010",
                          .(Avr = mean(Arithmetic.Mean)),
                          by = .(State.Code,County.Code,Site.Num,Parameter.Code,
                                 Latitude,Longitude,Datum,Pollutant.Standard,Metric.Used,
                                 Year,Units.of.Measure)]
  Annual_NO2_Summary <- rbind(Annual_NO2_Summary,Annual_NO2_Tmp)
}

#--- Preparation (Lon&Lat uniquely determine an annual record.)
#-- Remove locations away from mainland
Annual_CO_Summary <- Annual_CO_Summary[Longitude>=-130&Latitude>=20]
Annual_NO2_Summary <- Annual_NO2_Summary[Longitude>=-130&Latitude>=20]
#-- Remove different Datum： "NAD83" "WGS84"
Annual_CO_Summary <- Annual_CO_Summary[Datum=="WGS84"]
Annual_NO2_Summary <- Annual_NO2_Summary[Datum=="WGS84"]
#-- Transform State.Code
Annual_CO_Summary$State.Code <- as.numeric(Annual_CO_Summary[,State.Code])
Annual_NO2_Summary <- Annual_NO2_Summary[State.Code!="CC",]
Annual_NO2_Summary$State.Code <- as.numeric(Annual_NO2_Summary[,State.Code])


#=== Remove Data with different site
#-- Find Lon&Lat should be removed
Site_CO_Length <- Annual_CO_Summary[,.(Obs_site=length(unique(Site.Num))),
                                    by = .(Latitude,Longitude)]
Site_CO_Length <- Site_CO_Length[Obs_site==1]
#-- Just preserve the data with single site
Annual_CO_Summary <- Annual_CO_Summary[paste(Longitude,Latitude)%in%paste(Site_CO_Length$Longitude,
                                                                          Site_CO_Length$Latitude)]

#-- Find Lon&Lat should be removed
Site_NO2_Length <- Annual_NO2_Summary[,.(Obs_site=length(unique(Site.Num))),
                                    by = .(Latitude,Longitude)]
Site_NO2_Length <- Site_NO2_Length[Obs_site==1]
#-- Just preserve the data with single site
Annual_NO2_Summary <- Annual_NO2_Summary[paste(Longitude,Latitude)%in%paste(Site_NO2_Length$Longitude,
                                                                          Site_NO2_Length$Latitude)]

#-- Just preserve the data with full record
Annual_CO_Summary <- Annual_CO_Summary[Year>=2010]
Annual_NO2_Summary <- Annual_NO2_Summary[Year>=2010]

#--- CO increasing
Beta_CO_Stat <- Annual_CO_Summary[,
                                  .(beta_hat = beta.hat(Avr,Year),
                                    beta_hat_sd = beta.hat.sd(Avr,Year),
                                    beta_hat_norm_CO = beta.hat(Avr,Year)/beta.hat.sd(Avr,Year),
                                    State.Code=unique(State.Code),
                                    NumYear=length(Year)),
                                  by = .(Latitude,Longitude)]
Beta_CO_Stat <- Beta_CO_Stat[!is.na(beta_hat)]
Beta_NO2_Stat <- Annual_NO2_Summary[,
                                  .(beta_hat = beta.hat(Avr,Year),
                                    beta_hat_sd = beta.hat.sd(Avr,Year),
                                    beta_hat_norm_NO2 = beta.hat(Avr,Year)/beta.hat.sd(Avr,Year),
                                    State.Code=unique(State.Code),
                                    NumYear=length(Year)),
                                  by = .(Latitude,Longitude)]
Beta_NO2_Stat <- Beta_NO2_Stat[!is.na(beta_hat)]

##=== Summary Data
CO_Lat_Lon <- paste(Beta_CO_Stat$Latitude,Beta_CO_Stat$Longitude)
NO2_Lat_Lon <- paste(Beta_NO2_Stat$Latitude,Beta_NO2_Stat$Longitude)
CONO2_Lat_Lon <- unique(c(CO_Lat_Lon,NO2_Lat_Lon))
padding <- 0.15; label.size_1 <- 1;label.size_letter <- 5;p_size_l <- 3;p_size_s <- 2
#padding <- 0.1; label.size_1 <- 0.5;label.size_letter <- 4;p_size_l <- 2;p_size_s <- 1
for(hh in 2){
  for(beta0 in seq(-0.1,-0.5,by=-0.1)){
  #for(beta0 in c(-0.3)){  
  load(paste0("Fig/Ozone_h",hh,"_",beta0,".RData"))
    print(paste0("beta0: ",beta0))
    Ozone_Lat_Lon <- paste(T2_Stat$Latitude,T2_Stat$Longitude)
    Diff.ST <- paste(T2_Stat[Detect.ST%in%c("ST","2D(ST)")]$Latitude,
                     T2_Stat[Detect.ST%in%c("ST","2D(ST)")]$Longitude)
    Diff.IHW <- paste(T2_Stat[Detect.IHW%in%c("IHW","2D(IHW)")]$Latitude,
                     T2_Stat[Detect.IHW%in%c("IHW","2D(IHW)")]$Longitude)
    Diff.SA <- paste(T2_Stat[Detect.SA%in%c("SA","2D(SA)")]$Latitude,
                     T2_Stat[Detect.SA%in%c("SA","2D(SA)")]$Longitude)
    Diff.Total <- unique(c(Diff.ST,Diff.IHW,Diff.SA))
    
    #==== Find the locations need text.
    Text.ST <- Diff.ST[Diff.ST %in% CONO2_Lat_Lon]
    Text.ST <- T2_Stat %>% 
      filter(paste(Latitude,Longitude)%in%Text.ST) %>%
      select(Latitude,Longitude) %>%
      mutate(label=letters[1:length(Text.ST)],bg_CO="white",bg_NO2="white")
    Text.IHW <- Diff.IHW[Diff.IHW %in% CONO2_Lat_Lon]
    Text.IHW <- T2_Stat %>% 
      filter(paste(Latitude,Longitude)%in%Text.IHW) %>%
      select(Latitude,Longitude) %>%
      mutate(label=letters[1:length(Text.IHW)],bg_CO="white",bg_NO2="black")
    Text.SA <- Diff.SA[Diff.SA %in% CONO2_Lat_Lon]
    Text.SA <- T2_Stat %>% 
      filter(paste(Latitude,Longitude)%in%Text.SA) %>%
      select(Latitude,Longitude) %>%
      mutate(label=letters[1:length(Text.SA)],bg_CO="white",bg_NO2="black")
    
    Text.Total <- Diff.Total[Diff.Total %in% CONO2_Lat_Lon]
    Text.Total <- T2_Stat %>% 
      filter(paste(Latitude,Longitude)%in%Text.Total) %>%
      select(Latitude,Longitude) %>%
      mutate(label=letters[1:length(Text.Total)],bg_CO="white",bg_NO2="black")
    
    #==== Create results for CO and NO2
    print("CO increasing")
    CO_Res <- merge(Beta_CO_Stat[CO_Lat_Lon %in% Diff.Total],
                    T2_Stat[Ozone_Lat_Lon  %in% Diff.Total[Diff.Total%in%CO_Lat_Lon]],
                    by = "Latitude")
    print(CO_Res[,.(Latitude, Longitude.x, beta_hat_norm_CO,Detect.ST,Detect.IHW,Detect.SA)])
    print(CO_Res[,.(mean(beta_hat_norm_CO)),by="Detect.SA"])
    print(CO_Res[,.(mean(beta_hat_norm_CO)),by="Detect.ST"])
    print(CO_Res[,.(mean(beta_hat_norm_CO)),by="Detect.IHW"])
    CO_de_Lat <- as.numeric(CO_Res[order(beta_hat_norm_CO)[1],.(Latitude)])
    CO_de_Lon <- as.numeric(CO_Res[order(beta_hat_norm_CO)[1],.(Longitude.x)])
    Text.Total$bg_CO[Text.Total$Latitude==CO_de_Lat&Text.Total$Longitude==CO_de_Lon] <- "orange"
    
    print("NO2 increasing")
    NO2_Res <- merge(Beta_NO2_Stat[NO2_Lat_Lon %in% Diff.Total],
                    T2_Stat[Ozone_Lat_Lon  %in% Diff.Total[Diff.Total%in%NO2_Lat_Lon]],
                    by = "Latitude")
    print(NO2_Res[,.(Latitude, Longitude.x, beta_hat_norm_NO2,Detect.ST,Detect.IHW,Detect.SA)])
    
    print(NO2_Res[,.(mean(beta_hat_norm_NO2)),by="Detect.SA"])
    print(NO2_Res[,.(mean(beta_hat_norm_NO2)),by="Detect.ST"])
    print(NO2_Res[,.(mean(beta_hat_norm_NO2)),by="Detect.IHW"])
    
    NO2_de_Lat <- as.numeric(NO2_Res[order(beta_hat_norm_NO2)[1],.(Latitude)])
    NO2_de_Lon <- as.numeric(NO2_Res[order(beta_hat_norm_NO2)[1],.(Longitude.x)])
    Text.Total$bg_NO2[Text.Total$Latitude==NO2_de_Lat&Text.Total$Longitude==NO2_de_Lon] <- "green"
    print(Text.Total)
    #=== plot with vertex label
    p.SA <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.SA,shape=Detect.SA,size=Detect.SA),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.SA))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "SA" = "#1F78B4",
                                     "2D(SA)" = "#E31A1C",
                                     "SA & 2D(SA)" = "#6A3D9A",
                                     "black" = "black",
                                     "green" = "#33A02C"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "SA" = 16,
                                    "2D(SA)" = 17,
                                    "SA & 2D(SA)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "SA" = p_size_l,
                                   "2D(SA)" = p_size_l,
                                   "SA & 2D(SA)" = p_size_s))+
      #geom_label(aes(x=Longitude+c(0,0,1,-1,0,0,0,0,0,0),y=Latitude-1,# special case for beta0=-0.3
      geom_label(aes(x=Longitude,y=Latitude-1,
                     label=label,fill=bg_CO,col=bg_NO2),
                 data=Text.Total,size=label.size_letter,
                 label.padding = unit(padding,'cm'),
                 label.size = label.size_1,show.legend = F)+
      scale_fill_manual(values = c("white" = "white",
                                   "orange" = "#FDBF6F"))+
     # geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude+c(0,0,1,-1,0,0,0,0,0,0), yend = Latitude-0.5),# special case for beta0=-0.3
      geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude, yend = Latitude-0.5),
                   data=Text.Total,
                   arrow = arrow(length = unit(0.15, "cm")))+
      theme_void() +
      theme(legend.position="none") 
    p.SA
    p.SA.1 <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      theme_void() +
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.SA,shape=Detect.SA,size=Detect.SA),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.SA))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "SA" = "#1F78B4",
                                     "2D(SA)" = "#E31A1C",
                                     "SA & 2D(SA)" = "#6A3D9A"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "SA" = 16,
                                    "2D(SA)" = 17,
                                    "SA & 2D(SA)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "SA" = p_size_l,
                                   "2D(SA)" = p_size_l,
                                   "SA & 2D(SA)" = p_size_s))+
      theme(legend.position = "bottom",
            legend.title=element_text(size=15), 
            legend.text=element_text(size=15))
    
    # theme(legend.position = "bottom")
    #scale_colour_gradientn(colours=c("red","gray","blue"))
    #
    
    leg1 <- get_legend(p.SA.1)
    
    print(plot_grid(p.SA, leg1, nrow=2, rel_heights=c(1, .1)))
    #print(p.SA)
    
    ggsave(paste0("Fig/Ozone(SA)_h",hh,"_",beta0,".eps"),width = 7,height=5)
    
    
    p.ST <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.ST,shape=Detect.ST,size=Detect.ST),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.ST))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "ST" = "#1F78B4",
                                     "2D(ST)" = "#E31A1C",
                                     "ST & 2D(ST)" = "#6A3D9A",
                                     "black" = "black",
                                     "green" = "#33A02C"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "ST" = 16,
                                    "2D(ST)" = 17,
                                    "ST & 2D(ST)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "ST" = p_size_l,
                                   "2D(ST)" = p_size_l,
                                   "ST & 2D(ST)" = p_size_s))+
      #geom_label(aes(x=Longitude+c(0,0,1,-1,0,0,0,0,0,0),y=Latitude-1,# special case for beta0=-0.3
       geom_label(aes(x=Longitude,y=Latitude-1,
                     label=label,fill=bg_CO,col=bg_NO2),
                 data=Text.Total,size=label.size_letter,
                 label.padding = unit(padding,'cm'),
                 label.size = label.size_1,show.legend = F)+
      scale_fill_manual(values = c("white" = "white",
                                   "orange" = "#FDBF6F"))+
      #geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude+c(0,0,1,-1,0,0,0,0,0,0), yend = Latitude-0.5),# special case for beta0=-0.3
      geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude, yend = Latitude-0.5),
                   data=Text.Total,
                   arrow = arrow(length = unit(0.15, "cm")))+
      theme_void() +
      theme(legend.position="none") 
    
    p.ST.1 <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      theme_void() +
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.ST,shape=Detect.ST,size=Detect.ST),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.ST))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "ST" = "#1F78B4",
                                     "2D(ST)" = "#E31A1C",
                                     "ST & 2D(ST)" = "#6A3D9A"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "ST" = 16,
                                    "2D(ST)" = 17,
                                    "ST & 2D(ST)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "ST" = p_size_l,
                                   "2D(ST)" = p_size_l,
                                   "ST & 2D(ST)" = p_size_s))+
      theme(legend.position = "bottom",
            legend.title=element_text(size=15), 
            legend.text=element_text(size=15))
    
    # theme(legend.position = "bottom")
    #scale_colour_gradientn(colours=c("red","gray","blue"))
    #
    
    leg1 <- get_legend(p.ST.1)
    
    print(plot_grid(p.ST, leg1, nrow=2, rel_heights=c(1, .1)))
    #print(p.ST)
    ggsave(paste0("Fig/Ozone(ST)_h",hh,"_",beta0,".eps"),width = 7,height=5)
    
    
    p.IHW <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.IHW,shape=Detect.IHW,size=Detect.IHW),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.IHW))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "IHW" = "#1F78B4",
                                     "2D(IHW)" = "#E31A1C",
                                     "IHW & 2D(IHW)" = "#6A3D9A",
                                     "black" = "black",
                                     "green" = "#33A02C"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "IHW" = 16,
                                    "2D(IHW)" = 17,
                                    "IHW & 2D(IHW)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "IHW" = p_size_l,
                                   "2D(IHW)" = p_size_l,
                                   "IHW & 2D(IHW)" = p_size_s))+
      #geom_label(aes(x=Longitude+c(0,0,1,-1,0,0,0,0,0,0),y=Latitude-1,# special case for beta0=-0.3
      geom_label(aes(x=Longitude,y=Latitude-1,
                     label=label,fill=bg_CO,col=bg_NO2),
                 data=Text.Total,size=label.size_letter,
                 label.padding = unit(padding,'cm'),
                 label.size = label.size_1,show.legend = F)+
      scale_fill_manual(values = c("white" = "white",
                                   "orange" = "#FDBF6F"))+
      #geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude+c(0,0,1,-1,0,0,0,0,0,0), yend = Latitude-0.5),# special case for beta0=-0.3
      geom_segment(aes(x = Longitude, y = Latitude, xend = Longitude, yend = Latitude-0.5),
                   data=Text.Total,
                   arrow = arrow(length = unit(0.15, "cm")))+
      theme_void() +
      theme(legend.position="none") 
    
    p.IHW.1 <- ggplot() +
      geom_map(
        data = state, map = state,
        aes(map_id = region),
        color = "black", fill = "lightgray", size = 0.1
      )+
      theme_void() +
      geom_point(aes(x=Longitude,y=Latitude,
                     col=Detect.IHW,shape=Detect.IHW,size=Detect.IHW),
                 #size=1,
                 data=T2_Stat %>% arrange(Order.IHW))+
      scale_colour_manual(values = c("Null" = "#999999",
                                     "IHW" = "#1F78B4",
                                     "2D(IHW)" = "#E31A1C",
                                     "IHW & 2D(IHW)" = "#6A3D9A"))+
      scale_shape_manual(values = c("Null" = 3,
                                    "IHW" = 16,
                                    "2D(IHW)" = 17,
                                    "IHW & 2D(IHW)" = 2))+
      scale_size_manual(values = c("Null" = p_size_s,
                                   "IHW" = p_size_l,
                                   "2D(IHW)" = p_size_l,
                                   "IHW & 2D(IHW)" = p_size_s))+
      theme(legend.position = "bottom",
            legend.title=element_text(size=15), 
            legend.text=element_text(size=15))
    
    # theme(legend.position = "bottom")
    #scale_colour_gradientn(colours=c("red","gray","blue"))
    #
    
    leg1 <- get_legend(p.IHW.1)
    
    print(plot_grid(p.IHW, leg1, nrow=2, rel_heights=c(1, .1)))
    #print(p.IHW)
    
    ggsave(paste0("Fig/Ozone(IHW)_h",hh,"_",beta0,".eps"),width = 7,height=5)
    
  }
}

Avr_Ozone_Stat <- Annual_Summary[,.(avr_hat = mean(Avr),State.Code=unique(State.Code)),
                                 by = .(Latitude,Longitude)]
Avr_Ozone_Lat_Lon <- paste(Avr_Ozone_Stat$Latitude,Avr_Ozone_Stat$Longitude)

for(hh in 2){
  for(beta0 in seq(-0.1,-0.5,by=-0.1)){
    #file_num <- file_num+1
    load(paste0("Fig/Ozone_h",hh,"_",beta0,".RData"))
    print(paste0("beta0: ",beta0))
    Ozone_Lat_Lon <- paste(T2_Stat$Latitude,T2_Stat$Longitude)
    Diff.ST <- paste(T2_Stat[Detect.ST%in%c("ST","2D(ST)")]$Latitude,
                     T2_Stat[Detect.ST%in%c("ST","2D(ST)")]$Longitude)
    Diff.IHW <- paste(T2_Stat[Detect.IHW%in%c("IHW","2D(IHW)")]$Latitude,
                      T2_Stat[Detect.IHW%in%c("IHW","2D(IHW)")]$Longitude)
    Diff.SA <- paste(T2_Stat[Detect.SA%in%c("SA","2D(SA)")]$Latitude,
                     T2_Stat[Detect.SA%in%c("SA","2D(SA)")]$Longitude)
    Diff.Total <- unique(c(Diff.ST,Diff.IHW,Diff.SA))
    print("Avr Ozone increasing")
    Avr_Ozone_Res <- merge(Avr_Ozone_Stat[Avr_Ozone_Lat_Lon %in% Diff.Total],
                           T2_Stat[Ozone_Lat_Lon %in% Diff.Total[Diff.Total%in%Avr_Ozone_Lat_Lon]],
                           by = "Latitude")
    Avr_Ozone_All <- merge(Avr_Ozone_Stat,
                           T2_Stat,by = "Latitude")
    if(F){
      print(Avr_Ozone_Res[,.(mean(avr_hat)),by="Detect.SA"])
      print(Avr_Ozone_Res[,.(mean(avr_hat)),by="Detect.ST"])
      print(Avr_Ozone_Res[,.(mean(avr_hat)),by="Detect.IHW"])}
    if(T){
      print( Avr_Ozone_All[,.(mean(avr_hat)),by="Detect.SA"])
      print( Avr_Ozone_All[,.(mean(avr_hat)),by="Detect.ST"])
      print( Avr_Ozone_All[,.(mean(avr_hat)),by="Detect.IHW"])}
    #print(Avr_Ozone_Res[,.(Latitude, Longitude.x,avr_hat,beta_hat,Detect.ST,Detect.IHW,Detect.SA)])
   }
}
