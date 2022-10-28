#==================================================
#' Estimating adaptive pis via kernel regression
#==================================================

#====================
# 1 dimensional
#====================
pis_1D.func.ker.reg <- function(x, point, tau=0.1, h=10)
{
  ## pis_1D.func.ker.reg calculates the conditional proportions pis via kernel regression
  ## Arguments
  # x: a matrix of z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: conditional proportions
  m <- length(x)
  x.vec<-c(t(x))
  pv.vec<-1-pnorm(x.vec) #Change to one-side test
  scr.idx<-which(pv.vec>=tau)
  identify<-as.numeric(pv.vec>=tau)
  
  bw <- npregbw(formula=identify~point[,1],
                bandwidth.compute=TRUE)
  ker.reg <- npreg(bws = bw)
  
  p.plus.est <- ker.reg$mean
  p.est <- sapply( p.plus.est/(1-tau), function(x){min(x,1-1e-5)})
  
  pis.est<-1-c(p.est)
  return(pis.est)
}

#====================
# 2 dimensional
#====================
pis_2D.func.ker.reg <- function(x, point, tau=0.1, h=10)
{
  ## pis_2D.func calculates the conditional proportions pis via kernal regression
  ## Arguments
  # x: a matrix of z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: conditional proportions
  m <- length(x)
  x.vec<-c(t(x))
  pv.vec<-1-pnorm(x.vec) #Change to one-side test
  scr.idx<-which(pv.vec>=tau)
  identify<-as.numeric(pv.vec>=tau)
  
  bw <- npregbw(formula=identify~point[,1]+point[,2],
                bandwidth.compute=TRUE)
  ker.reg <- npreg(bws = bw)
  
  p.plus.est <- ker.reg$mean
  p.est <- sapply( p.plus.est/(1-tau), function(x){min(x,1-1e-5)})
  
  pis.est<-1-c(p.est)
  return(pis.est)
}