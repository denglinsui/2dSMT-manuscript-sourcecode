library(ggplot2)
library(latex2exp)
source("Spatial_Detection.R")
set.seed(10)
n <- 10
k <- 2
x <- sample(n);y <- sample(n)
T1 <- x; T2 <- y
cutoffs <- cutoff.gen(x,y)
pt <- data.frame(T1=x, T2=y)
threshod <- data.frame(T1=cutoffs[-1,1], T2=cutoffs[-1,2])

# Total Set
all_threshod <- data.frame(T1=rep(T1,times = length(T2)), 
                           T2=rep(T1,each = length(T2)))
p_all <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(data = all_threshod,aes(x=T1,y=T2),
             shape = 0,size =2,
             color = "blue")+
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_all
ggsave("Tools/Reduced_Set/ALL_Threshod.pdf", 
       plot = last_plot(),
       width = 4, height = 4)

# Total Set
all_threshod_0 <- data.frame(T1=rep(k,times = length(T2)), 
                           T2=T2)
p_k_thre_0 <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(data = all_threshod_0,aes(x=T1,y=T2),
             shape = 0,size =2,
             color = "blue")+
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_k_thre_0
ggsave("Tools/Reduced_Set/k_Threshod_0.pdf", 
       plot = last_plot(),
       width = 4, height = 4)

# Step_k_0to1
p_k_thre_0to1 <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(aes(x=c(2),y=c(5)),
             shape = 0,size =2,
             color = "blue") + 
  geom_point(aes(x=c(2),y=c(6)),
             size = 2,
             color = "purple",
             shape = 2) +   
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_k_thre_0to1
ggsave("Tools/Reduced_Set/k_Threshod_0to1.pdf", 
       plot = last_plot(),
       width = 4, height = 4)

# Step_k_1
k4_T2 <- T2[T1>=k]
k4_threshod_1 <- data.frame(T1=rep(k,times = length(k4_T2)), 
                            T2=k4_T2)
p_k_thre_1 <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(data = k4_threshod_1,aes(x=T1,y=T2),
             size = 2,
             color = "purple",
             shape = 2)+
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_k_thre_1
ggsave("Tools/Reduced_Set/k_Threshod_1.pdf", 
       plot = last_plot(),
       width = 4, height = 4)

# Step_k_1to2
p_k_thre_1to2 <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(aes(x=c(2),y=c(4)),
             size = 2,
             color = "purple",
             shape = 2)+
  geom_point(aes(x=c(4),y=c(4)),
             size = 5,
             color = "red",
             shape = "o") +     
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_k_thre_1to2
ggsave("Tools/Reduced_Set/k_Threshod_1to2.pdf", 
       plot = last_plot(),
       width = 4, height = 4)

# Step_k_3
T2_0 <- T2[T1==k]
k4_T2_2 <- T2[T1>=k & T2<=T2_0]
k4_threshod_2 <- data.frame(T1=rep(k,times = length(k4_T2_2)), 
                            T2=k4_T2_2)
p_k_thre_2 <- ggplot(data=pt,aes(x=T1,y=T2))+
  geom_point()+
  geom_point(data = k4_threshod_2,aes(x=T1,y=T2),
             shape = 'o',size =5,
             color = "red")+
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_k_thre_2
ggsave("Tools/Reduced_Set/k_Threshod_2.pdf", 
       plot = last_plot(),
       width = 4, height = 4)
# Final Set
p_final <- ggplot()+
  geom_point(data=pt,aes(x=T1,y=T2))+
  geom_point(data = threshod,aes(x=T1,y=T2),
             shape = 'o',size =5,
             color = "red")+
  xlab(TeX("$\\widehat{T}_2(s)$"))+
  ylab(TeX("$\\widehat{T}_1(s)$"))
p_final
ggsave("Tools/Reduced_Set/Final_Threshod.pdf", 
       plot = last_plot(),
       width = 4, height = 4)


# Combine
library(ggpubr)
p_arr_1 <- ggarrange(p_all,
          p_final+ 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank()),
          labels = "AUTO",
          nrow = 1, ncol = 2)
p_arr_1
ggsave("Tools/Reduced_Set/ALL_Final.pdf", 
       plot = last_plot(),
       width = 8, height = 4)

ggsave("Tools/Reduced_Set/ALL_Final.eps", 
       plot = last_plot(),
       width = 8, height = 4)


p_arr_2 <- ggarrange(p_k_thre_0,
                     p_k_thre_0to1+ 
                       theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank()),
                     p_k_thre_1+ 
                       theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank()),
                     p_k_thre_1to2+ 
                       theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank()),
                     p_k_thre_2+ 
                       theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank()),
          labels = "AUTO",
          nrow = 1, ncol = 5)
p_arr_2
ggsave("Tools/Reduced_Set/K_Threshold.pdf", 
       plot = last_plot(),
       width = 20, height = 4)

ggarrange(p_arr_1,p_arr_2,
          nrow = 2,ncol=1)

