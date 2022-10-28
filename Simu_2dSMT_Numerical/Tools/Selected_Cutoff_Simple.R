#' Simulation for selected cutoff

##=== Start of problem independent section
rm(list=ls())
set.seed(10)
options(digits=5)
require(np)
require(latex2exp)
require(ggplot2)
setwd('~/project/2DSpatial')
#setwd("D:/RUC/project/multiple\ testing/2DSpatial")
source('Spatial_Detection.R')
source('Tools/Basic.R')
source('Tools/Fix_GenData.R')

#=== Initialize Parameters
q <- 0.1
n <- 1 # sample size
m <- 300 # point size
data.type <- "SquareSpline"
k <- 1
rho.eps <- 0.1
cor.type <- "cor"
tau <- 0.5

#=== parameters for 2D
detect.m <- "top.k"
hh <- 5
const <- 0
pis <- rep(1,m)

#=== Setting
magnitude <- 2
r <- 0

point <- matrix(seq(0,30,length.out=m), 
                m, 1, byrow = F) 

#=== Initial Dist and Sigma
Dist.p <- Dist(m)
Sigma.eps.p <- Sigma.eps(m, r, k, rho.eps, point)
Sigma.eps.est <- Sigma.eps.p

#=== Generate signal
mu.res <- mu.fix.gen.1D(point, "uc.spline", magnitude)
mu <- mu.res$mu
pluszero <- mu>0
sgm <- sqrt(Sigma.eps.est[1,1])

#=== Generate Data
X <- MASS::mvrnorm(n = n, mu = mu, 
                   Sigma = Sigma.eps.p)
X <- matrix(X, nrow = n)

#==== Perform Algorithm
dig <- 7
T1 <- apply(X,2,function(x){sum(x)/sgm/sqrt(n)})
T1 <- round(T1,dig)

#================================
#==== Without spatial Info: One D
#================================
res.1D <- OneD_Detect(T1, q,
                      const = const)
T1.star <- res.1D$t1.min 
max.rej <- res.1D$max.rej

#==== Set pre-chosen T1, T2
T1.star <- round(T1.star,dig)
T2.star <- Inf

#=== Detect Neighbor
#hh <- 2
hh <- 5
hh.seq <- rep(hh,m)
Neigh_Detect_res <- Neigh_Detect(hh = hh.seq,
                                 X = X, 
                                 Dist = Dist.p,
                                 Sigma.eps = Sigma.eps.est,
                                 detect.m = detect.m)
T2 <- Neigh_Detect_res$T2
V2 <- Neigh_Detect_res$V2
V1V2.cov <- Neigh_Detect_res$V1V2.cov
ind <- Neigh_Detect_res$ind

#== round numerical values to ensures the correctness of taking maximum
T1 <- round(T1,dig)
T2 <- round(T2,dig)
V2 <- round(V2,dig)
V1V2.cov <- round(V1V2.cov,dig)

  
eta <- T2
#=== Perform Non-parametric Emprical Bayesian
mm <- REBayes::GLmix(x = eta[ind])
normalized.prob <- mm$y / sum(mm$y)
  
  #=== Estimate false discovery
  fdr.est <- function (t1, t2, NP.max, etype){
    p <- length(T1)
    L <- numeric(p)
    NP <- sum(T1 >= t1 & T2 >= t2)
    NP <- ifelse(is.na(NP), 0, NP)
    
    if (NP == 0){
      FDP <- 0
      FD <- 0
    }else if(NP<NP.max){
      FDP <- NA
      FD <- NA
    }else{
      #group_by (T2,V1V2.cov)
      grp.val <- unique(cbind(V2,V1V2.cov))
      for(j in 1:nrow(grp.val)){
        vV2 <- grp.val[j,1]
        vV1V2.cov <- grp.val[j,2]
        row.ind <- which(V2==vV2 & V1V2.cov==vV1V2.cov)
        L[row.ind] <- L.cal(t1,t2,
                            mm,normalized.prob,
                            vV2,vV1V2.cov)
      }
      
      if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        FD <- sum(pis * L)
        FDP <- FD / NP
        #print(paste(FDP,p*(1-B2)/NP))
      }
    } 
    return(list(FDP = FDP, 
                NP = NP,
                FD = FD))	
  }
  
#=== Initialize the searching grid
cutoff.rec <- cutoff.gen.rec(T1,T2,T1.star,T2.star)
cutoff <- cutoff.rec$cutoff
ind.T1 <- cutoff.rec$ind.T1
  
  FDP <- NULL
  NP <- NULL
  #t1t2 <- NULL
  t1.cand.set <- NULL
  t2.cand.set <- NULL
  cutoff_search <- NULL
  NP.max <- max.rej
  
  index.T1 <- max.rej+1
  
  while(index.T1 <= m+1){
    i.up <- ind.T1[index.T1+1]-1
    cur.ind.T1 <- ind.T1[index.T1]
    i.down <- cur.ind.T1
    
    # We only consider the rej > current max rej
    t1.down <- cutoff[i.down,1]
    t2.down <- cutoff[i.down,2]
    
    NP.down <- sum(T1 >= t1.down & T2 >= t2.down)
    i.down <- i.down + max(0, NP.max - NP.down)
    
    while(i.up >= i.down){
      
      #=== Search Down
      # If not stand, we change i.down to accelerate
      t1.down <- cutoff[i.down,1]
      t2.down <- cutoff[i.down,2]
      cutoff_search <- rbind(cutoff_search, c(t1.down, t2.down))
      obj <- fdr.est(t1.down, t2.down, NP.max, etype="FDR")
      
      if(is.na(obj$FDP)){
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        t1.cand.set <- c(t1.cand.set, t1.down)
        t2.cand.set <- c(t2.cand.set, t2.down)
        NP.max <- obj$NP
        t10 <- t1.down
        t20 <- t2.down
        
        i.down <- i.down +1
      } else{
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        
        i.down <- i.down + min.REJ - NP.down
      }
    }
    index.T1 <- index.T1 + 1
  }
  
data <- data.frame(x = c(cutoff[,1],cutoff_search[,1],t1.cand.set),
                   y = c(cutoff[,2],cutoff_search[,2],t2.cand.set),
                   type = c(rep("Initialized Candidate Set",nrow(cutoff)),
                            rep("Pruned Candidate Set",nrow(cutoff_search)+length(t1.cand.set))),
                   color = c(rep("cutoff",nrow(cutoff)+nrow(cutoff_search)),
                             rep("active cutoff",length(t1.cand.set))))
data_item <- data.frame(x=T1,
                        y=T2,
                        type="Scatter plot of $T_1$ vs $T_2$",
                        color = ifelse(mu>0,"NonNull","Null"))
data_merge <- rbind(data,data_item)
#data <- data[data$x>0,]
levels(data_merge$type) <- c(
                             "Initialized Candidate Set"=
                               TeX("Initialized Candidate Set$\\ $"),
                             "Pruned Candidate Set"=
                               TeX("Pruned Candidate Set$\\ $"),
                             "Scatter plot of $T_1$ vs $T_2$"=
                               TeX("Scatter plot of $T_1$ vs $T_2$"))
data_merge$type <- factor(data_merge$type,
                                  levels = c( TeX("Scatter plot of $T_1$ vs $T_2$"),
                                              TeX("Initialized Candidate Set$\\ $"),
                                              TeX("Pruned Candidate Set$\\ $")
                                              ))
typeplot <- data.frame(type=c("Scatter plot of $T_1$ vs $T_2$"))
levels(typeplot$type) <- c("Scatter plot of $T_1$ vs $T_2$"=
                             TeX("Scatter plot of $T_1$ vs $T_2$"))

typeplot <- data.frame(type=c("Scatter plot of $T_1$ vs $T_2$",
                              "Pruned Candidate Set",
                              "Initialized Candidate Set"))


levels(typeplot$type) <- c("Initialized Candidate Set"=
                             TeX("Initialized Candidate Set$\\ $"),
                           "Pruned Candidate Set"=
                             TeX("Pruned Candidate Set$\\ $"),
                           "Scatter plot of $T_1$ vs $T_2$"=
                             TeX("Scatter plot of $T_1$ vs $T_2$"))


t1 <- t1.cand.set[which.max(NP)]
t2 <- t2.cand.set[which.max(NP)]
t1.1D <- min(T1[p.adjust(1-pnorm(T1),method="BH") <=q])
test <- data.frame(Tm=T1,Ta=T2,mu=mu,hypo=mu>0)
fdp(which(T1>=t1.1D),mu);Pow(which(T1>=t1.1D),mu)
fdp(which(T1>=t1&T2>=t2),mu);Pow(which(T1>=t1&T2>=t2),mu)

Cutoff_Set <- 
  ggplot(data = data_merge, aes(x=x,y=y,color=color,shape=color,size=color))+
  geom_point()+
  #facet_grid(. ~ type)+
  facet_grid(
    #rows=vars(type),
    cols=vars(type),
    #scales = "free_x",
    labeller=label_parsed
  )+
  scale_color_manual(values  = c("cutoff"="#4DAF4A",
                                 "active cutoff"= "#FF7F00",
                                 "Null"="#1F78B4",
                                 "NonNull"="#E31A1C"))+
  scale_shape_manual(values  = c(4,16,2,1))+
  scale_size_manual(values  = c(2,1,2,2))+
  xlab(TeX("$t_2$"))+
  ylab(TeX("$t_1$"))+
  theme_bw()+
  annotate("rect", xmin=t1, xmax=Inf, ymin=t2, ymax=Inf, fill = "#FB9A99", alpha = 0.1)+
  geom_rect(data=typeplot, 
            aes(xmin=t1, xmax=Inf, ymin=t2, ymax=Inf), 
            inherit.aes = FALSE,
            fill = "#FB9A99", alpha = 0.2)+
  geom_rect(data=typeplot, 
            aes(xmin=t1.1D, xmax=Inf,ymin=-Inf, ymax=Inf), 
            inherit.aes = FALSE,
            fill = "#A6CEE3", alpha = 0.2)+
  geom_vline(data=typeplot, 
             aes(xintercept = t1.1D),colour= "#A6CEE3",size=1,linetype= "dashed")+
  geom_vline(data=typeplot, 
             aes(xintercept = t1),colour= "#FB9A99",size=1)+
  geom_hline(data=typeplot, 
             aes(yintercept = t2),colour= "#FB9A99",size=1)+
  annotate("rect", xmin=t1.1D, xmax=Inf, 
           ymin=-Inf, ymax=Inf, fill = "#A6CEE3", alpha = 0.1)+
  theme(legend.position = "none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        strip.text.x = element_text(size = 15))

Cutoff_Set  



ggsave(file="Figure/SelectCutoff3.eps", 
       plot = last_plot(),
       width = 12, height = 4)
    

ggsave(file="Figure/SelectCutoff3.pdf", 
       plot = Cutoff_Set,
       width = 9, height = 4)

t1 <- t1.cand.set[which.max(NP)]
t2 <- t2.cand.set[which.max(NP)]
t1.1D <- min(T1[p.adjust(1-pnorm(T1),method="BH") <=q])
test <- data.frame(Tm=T1,Ta=T2,mu=mu,hypo=mu>0)

p_2Dillu <- ggplot(data=test,aes(Tm,Ta,color=hypo,shape=hypo))+
  geom_point()+
  annotate("rect", xmin=t1, xmax=Inf, ymin=t2, ymax=Inf, fill = "#FB9A99", alpha = 0.1)+
  annotate("rect", xmin=t1.1D, xmax=Inf, 
           ymin=-Inf, ymax=Inf, fill = "#A6CEE3", alpha = 0.1)+
  scale_color_manual(values = c("#1F78B4","#E31A1C"))+
  geom_vline(xintercept = t1.1D,colour= "#A6CEE3",size=1)+
  geom_vline(xintercept = t1,colour= "#FB9A99",size=1)+
  geom_hline(yintercept = t2,colour= "#FB9A99",size=1)+
  theme_bw()+
  ylab("T1")+xlab("T2")+
  scale_shape_manual(values = c(1,2))+
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title=element_text(size=15),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))

ggarrange(p_2Dillu,Cutoff_Set)
nrow(cutoff)
nrow(cutoff_search)

