## Functions used here.
bh.func<-function(pv, q)
{ 
  # the input 
  # pv: the p-values
  # q: the FDR level
  # the output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)   
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}

disvec.func<-function(dims, s)
{
  # disvec computes the distances of all points on a m1 times m2 spatial domain to a point s
  ## Arguments:
  # dims=c(d1, d2): the dimensions
  # s=c(s1, s2): a spatial point
  ## Values:
  # a vector of distances
  
  m<-dims[1]*dims[2]
  dis.vec<-rep(0, m)
  for(i in 1:dims[1])
  {
    dis.vec[((i-1)*dims[2]+1):(i*dims[2])]<-sqrt((i-s[1])^2+(1:dims[2]-s[2])^2)
  }
  return(dis.vec) 
}


pis_2D.func <- function(x, tau=0.1, h=10)
{
  ## pis_2D.func calculates the conditional proportions pis
  ## Arguments
  # x: a matrix of z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: conditional proportions
  
  dims<-dim(x)
  m<-dims[1]*dims[2]
  x.vec<-c(t(x))
  pv.vec<-1-pnorm(x.vec) #Change to one-side test
  scr.idx<-which(pv.vec>=tau)
  p.est<-matrix(rep(0, m), dims[1], dims[2])  
  
  for (i in 1:dims[1]) 
  {
    for (j in 1:dims[2]) 
    {
      s<-c(i, j)
      dis.vec<-disvec.func(dims, s)
      kht<-dnorm(dis.vec, 0, h)
      p.est[i,j]<-min(1-1e-5, sum(kht[scr.idx])/((1-tau)*sum(kht)))
    }
  }
  pis.est<-1-c(p.est)
  return(pis.est)
}

law.func<-function(pvs, pis, q)
{
  ## implementing "spatial multiple testing by locally adaptive weighting"
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values
  
  m<-length(pvs)
  nu<-10e-5
  pis[which(pis<nu)]<-nu # stabilization
  pis[which(pis>1-nu)]<-1-nu # stabilization
  ws<-pis/(1-pis)
  pws<-pvs/ws
  st.pws<-sort(pws)
  fdps<-sum(pis)*st.pws/(1:m)
  de<-rep(0, m)
  if(sum(fdps<=q)==0)
  {
    k<-0
    pwk<-1
  }
  else
  {
    k<-max(which(fdps<=q))
    pwk<-st.pws[k]
    de[which(pws<=pwk)]<-1
  }
  y<-list(nr=k, th=pwk, de=de)
  return (y)
}	 


sab.func<-function(pvs, pis, q)
{
  ## implementing "SABHA" by Li and Barber
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values
  
  m<-length(pvs)
  nu<-10e-5
  pis[which(pis>1-nu)]<-1-nu # stabilization
  pws<-pvs*(1-pis)
  st.pws<-sort(pws)
  
  pwi<-st.pws/1:m
  de<-rep(0, m)
  if (sum(pwi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pwi<=(q/m)))
    pk<-st.pws[k]
    de[which(pws<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}	 

pis_1D.func<- function(x, tau=0.1, h=50)
{
  ## pis_est.func calculates the conditional proportions pis
  ## Arguments
  # x: z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: the conditional proportions
  
  m <- length(x)
  s <- 1:m # auxiliary variable
  pval <- 2*pnorm(-abs(x))
  p.est <-rep(0, m)
  for (i in 1:m) { 
    kht<-dnorm(s-i, 0, h)
    p.est[i]<-sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <-1
  return(1-p.est)
}
