### Describe the simulation setting

library(latex2exp)
library(ggplot2)
source("Tools/Fix_GenData.R")
source("Tools/Basic.R")
#Initialize point
###================== One Dimensional: mu (Spline) ====================
m <- 900 # point size

point <- matrix(seq(0,30,length.out=m), 
                m, 1, byrow = F) 
magnitude <- 1
method <- "uc.spline"


mu1 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=1,single_bump_prop=1/20)
mu2 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=2,single_bump_prop=1/20)
mu3 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=4,single_bump_prop=1/20)

data <- data.frame(mu=c(mu1$mu,mu2$mu,mu3$mu),
                   x=rep(point,times=3),
                   Type=factor(rep(c("Sparse", "Medium", "Dense"),each=m),
                               levels =c("Sparse", "Medium", "Dense"))
)

p <- 
  ggplot(data,aes(x=x,y=mu))+
  geom_line()+
  facet_grid(.~Type)+ theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  xlab("s")+ylab(TeX("\\mu(s)"))

cairo_ps("Figure/1D_mu.eps",width = 13,height=5)
p
dev.off()



#ggsave(file="Figure/1D_mu.eps",width = 9,height=3)

###======================================

#Initialize point
###================== One Dimensional: mu(Unif) ====================
m <- 900 # point size

point <- matrix(seq(0,30,length.out=m), 
                m, 1, byrow = F) 
magnitude <- 1
method <- "uc.unif"

set.seed(10)
mu1 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=1,single_bump_prop=1/20)

set.seed(10)
mu2 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=2,single_bump_prop=1/20)

set.seed(10)
mu3 <- mu.fix.gen.1D(point, method=method, magnitude,num_bump=4,single_bump_prop=1/20)

data <- data.frame(mu=c(mu1$mu,mu2$mu,mu3$mu),
                   x=rep(point,times=3),
                   Type=factor(rep(c("Sparse", "Medium", "Dense"),each=m),
                               levels =c("Sparse", "Medium", "Dense"))
)

p <- 
  ggplot(data,aes(x=x,y=mu))+
  geom_point(shape=1,size=1,alpha=1)+
  facet_grid(.~Type)+ theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  xlab("s")+ylab(TeX("\\mu(s)"))


cairo_ps("Figure/1DUnif_mu.eps",width = 13,height=5)
p
dev.off()



#ggsave(file="Figure/1DUnif_mu.eps",width = 9,height=3)

###======================================

###======================================

#Initialize point
###================== One Dimensional: mu(mvnorm) ====================
m <- 900 # point size

point <- matrix(seq(0,30,length.out=m), 
                m, 1, byrow = F) 
magnitude <- 1
method <- "mvnorm"

Dist.p <- as.matrix(dist(point))
Sigma.mu <- K.mu.d(Dist.p, k =1, rho.mu = 0.3, sig2.mu=3)

set.seed(10)
mu1 <- mu.fix.gen.1D(point, method=method, magnitude,mu.mean = -2.5,Sigma.mu = Sigma.mu)

set.seed(10)
mu2 <- mu.fix.gen.1D(point, method=method, magnitude,mu.mean = -2,Sigma.mu = Sigma.mu)

set.seed(10)
mu3 <- mu.fix.gen.1D(point, method=method, magnitude,mu.mean = -1,Sigma.mu = Sigma.mu)

data <- data.frame(mu=c(mu1$mu,mu2$mu,mu3$mu),
                   x=rep(point,times=3),
                   Type=factor(rep(c("Sparse", "Medium", "Dense"),each=m),
                               levels =c("Sparse", "Medium", "Dense"))
)

p <- 
  ggplot(data,aes(x=x,y=mu))+
  geom_point(shape=1,size=1,alpha=1)+
  facet_grid(.~Type)+ theme_bw() +
  xlab("s")+ylab(TeX("\\mu(s)"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15))

cairo_ps("Figure/1Dmvnorm_mu.eps",width = 13,height=5)
p
dev.off()

#ggsave(file="Figure/1Dmvnorm_mu.eps",width = 9,height=3)

###======================================

###================== One Dimensional + Two Dimensional: Cov ====================

r <- c(0.1,0.3,0.3)
k <- c(1,1,2)
rho.eps <- c(0.1,0.1,0.3)

r <- c(0.5,0.8,0.6)
k <- c(1,1,2)
rho.eps <- c(0.05,0.1,0.2)

for(ii in 1:3){
  assign(paste("Corr",ii,sep = ""),
         sapply(1:m, 
                function(i){
                  K.eps(point[1, ], point[i, ], r=r[ii], k=k[ii], rho.eps = rho.eps[ii])
                }
         )
  )
}
#Corr1 <- sapply(1:m, function(i){K.eps(point[1, ], point[i, ], r=r[1], k=k[1], rho.eps = rho.eps[1])})
#Corr2 <- sapply(1:m, function(i){K.eps(point[1, ], point[i, ], r=r[2], k=k[2], rho.eps = rho.eps[2])})
#Corr3 <- sapply(1:m, function(i){K.eps(point[1, ], point[i, ], r=r[3], k=k[3], rho.eps = rho.eps[3])})


data <- data.frame(Corr=c(Corr1,Corr2,Corr3),
                   distance=rep(point,times=3),
                   Type=factor(rep(paste("r=",r," k=",k," rho.eps=",rho.eps,sep=""),each=m))
)

#levels(data$Type) <- c("r=0.5 k=1 rho.eps=0.05" = TeX("$r=0.1$ $k=1$ $\\rho_{\\epsilon}=0.1$"), 
#                       "r=0.8 k=1 rho.eps=0.1" = TeX("$r=0.3$ $k=1$ $\\rho_{\\epsilon}=0.1$"), 
#                       "r=0.6 k=2 rho.eps=0.2" = TeX("$r=0.5$ $k=2$ $\\rho_{\\epsilon}=0.3$"))

levels(data$Type) <- c("r=0.5 k=1 rho.eps=0.05" = "Weak", 
                       "r=0.6 k=2 rho.eps=0.2" = "Strong", 
                       "r=0.8 k=1 rho.eps=0.1" = "Medium")
data$Type <- factor(data$Type, levels=c("Weak","Medium","Strong"))

data1 <- 
  data %>%
  filter(distance<0.8)
p <- 
  ggplot(data1,aes(x=distance,y=Corr))+
  geom_point()+
  facet_grid(
    cols=vars(Type),
    labeller=label_parsed)+ theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  xlab("||s-t||")+ylab(TeX("$cov(\\epsilon(s),\\epsilon(t))$"))


cairo_ps("Figure/1D_Cov.eps",width = 13,height=5)
p
dev.off()


#ggsave(file="Figure/1D_Cov.eps",width = 9,height=3)

#============================================================


#Initialize point
###================== Two Dimensional: mu ====================
m <- 900 # point size
m_x <- 30
m_y <- 30
point_x <- seq(0,5,length.out=m_x)
point_y <- seq(0,5,length.out=m_x)
point <- cbind(rep(point_x,times=length(point_y)), 
               rep(point_y,each=length(point_x))) 

magnitude <- 1
mu1 <- mu.fix.gen.2D(point, method="uc.square", magnitude)

set.seed(10)
mu2 <- mu.fix.gen.2D(point, method="uc.circle", magnitude)

set.seed(10)
mu3 <- mu.fix.gen.2D(point, method="uc.mix", magnitude)

data <- data.frame(mu=c(mu1$mu,mu2$mu,mu3$mu),
                   x=rep(point[,1],times=3),
                   y=rep(point[,2],times=3),
                   Type=factor(rep(c("Sparse", "Medium", "Dense"),each=m),
                               levels =c("Sparse", "Medium", "Dense"))
)

#library(grDevices) #绘图颜色相关x
#library(RColorBrewer)#绘图颜色相关
#rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
#colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)

p <- 
  ggplot(data,aes(x=x,y=y,z=mu))+
  geom_tile(aes(fill=mu))+
  facet_grid(.~Type)+
  xlab(as.expression(TeX('$s_1,$')))+
  ylab(as.expression(TeX('$s_2$')))+
  scale_fill_gradient(as.expression(TeX('$\\mu(s_1,s_2)$')),
                      low="white", high="red") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15))

cairo_ps("Figure/2D_mu.eps",width = 13,height=5)
p
dev.off()



#ggsave(file="Figure/2D_mu.eps",width = 9,height=3)

###======================================
