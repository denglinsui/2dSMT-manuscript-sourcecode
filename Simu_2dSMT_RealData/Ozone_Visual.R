#=== Visualization
library(ggplot2)
library(ggmap)

#=== Original Visualization
state <- map_data("state")
ggplot() +
  geom_map(
    data = state, map = state,
    aes(map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )+
  theme_void() +
  geom_point(aes(x=Longitude,y=Latitude,color=beta_hat),data=T2_Stat)+
  scale_colour_gradientn(colours=c("red","gray","blue"))+
  theme(legend.position = "bottom")

ggsave("Fig/test_stat.eps",width = 4,height=3)


ggplot() +
  geom_map(
    data = state, map = state,
    aes(map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )+
  theme_void() +
  geom_point(aes(x=Longitude,y=Latitude,color=T2),data=T2_Stat)+
  scale_colour_gradientn(colours=c("red","gray","blue"))+
  theme(legend.position = "bottom")

ggsave("Fig/beta_hat.eps",width = 4,height=3)


ggplot() +
  geom_map(
    data = state, map = state,
    aes(map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )+
  theme_void() +
  geom_point(aes(x=Longitude,y=Latitude,color=beta_hat/beta_sd),data=T2_Stat)+
  scale_colour_gradientn(colours=c("red","gray","blue"))+
  theme(legend.position = "bottom")

ggsave("Fig/norm_beta_hat.eps",width = 4,height=3)

#=== Visualization for result
#T2_Stat$Detect <- sample(c(TRUE,FALSE),nrow(T2_Stat),replace  = TRUE)

ggplot() +
  geom_map(
    data = state, map = state,
    aes(map_id = region),
    color = "black", fill = "lightgray", size = 0.1
  )+
  theme_void() +
  geom_point(aes(x=Longitude,y=Latitude,color=Detect,shape=Detect),
             data=T2_Stat)+
  theme(legend.position = "bottom")
  #scale_colour_gradientn(colours=c("red","gray","blue"))+
  #
