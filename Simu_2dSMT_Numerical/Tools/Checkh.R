# Simulating spatial data

# S: 2D spatial domain with m1*m2 points


library('foreach')
library('doParallel')

setwd('~/project/Multiple-Testing-Replication/TDFDR_X/Spatial/Tools')
source('Basic.R')
source('laws_used.R')

#numCores = detectCores()
registerDoParallel(32)


dims<-c(200, 200)
m<-dims[1]*dims[2]
thetaS<-matrix(rep(0, m), dims[1], dims[2]) # true states of nature
pis<-matrix(0, dims[1], dims[2])

# A and B: 2D region with signals
# A is a triangle
# B is a rectangle 

pisA.func<-function(i, j)
{
  # decide whether (i, j) is in triangle A
  if((j <= 1.25*i+75) && (j >= 100) && j <= -1.25*i+225)
  {y <- 0.9}
  else if((j <= 1.25*i-25) && (j >= 100) && j <= -1.25*i+325)
  {y = 0.6}
  else
  {y = 0.01}
  return(y)
}

pisB.func<-function(i, j)
{
  # decide whether (i, j) is in the rectangle B
  if(i>=20 && i<=60 && j>=40 && j<=90)
  {y <- 0.9}
  else if(i>60 && i<=100 && j>=40 && j<=90)
  {y = 0.6}
  else
  {y<-0.01}
  return(y)
}

q<-0.1


mu0.vec<-seq(from=2.5, to=4, by=0.5)
np<-length(mu0.vec)
nrep<-10
locs<-1:m
loop_matrix = expand.grid(1:np, 1:nrep)

# Methods

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)

law.dd.fdr<-rep(0, np)
law.dd.etp<-rep(0, np)

tansey.fdr<-rep(0, np)
tansey.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.or.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)

tansey.fdp<-matrix(rep(0, nrep*np), np, nrep)
tansey.ntp<-matrix(rep(0, nrep*np), np, nrep)


system.time(for_out1 <- foreach(iter = 1:(np*nrep), .combine='rbind') %dopar%
              {
                i<-as.numeric(loop_matrix[iter,][1])
                j<-as.numeric(loop_matrix[iter,][2])
                
                cat("\n", "iteration i = ", i, "\n", "iteration j = ", j, "\n")
                mu0<-mu0.vec[i]
                
                set.seed(iter+1000)
                for(i1 in 1:dims[1])
                {
                  for (j1 in 1:dims[2])
                  {
                    pis[i1,j1] = pisB.func(i1,j1)
                    thetaS[i1,j1]<-rbinom(1, 1, pis[i1,j1])
                  }
                }
                
                thetaS.vec<-c(thetaS)
                pis.vec = c(pis)
                
                pii<-sum(thetaS.vec)/m
                x0<-rnorm(m, mean=0, sd=1)
                x1<-rnorm(m, mean=mu0, sd=1)
                x.vec<-(1-thetaS.vec)*x0+thetaS.vec*x1
                
                x<-matrix(x.vec, dims[1], dims[2])
                
                pv.vec<-2*pnorm(-abs(x.vec), 0, 1)
                bh.th<-bh.func(pv.vec, 0.9)$th
                pis.hat<-pis_2D.func(x, tau=bh.th, h=10)
                
                
                bh.res<-bh.func(pv.vec, q)
                bh.de<-bh.res$de
                bh.fdp[i, j]<-sum((1-thetaS.vec)*bh.de)/max(sum(bh.de), 1)
                bh.ntp[i, j]<-sum(thetaS.vec*bh.de)/sum(thetaS.vec)
                
                law.or.res<-law.func(pvs=pv.vec, pis.vec, q)
                law.or.de<-law.or.res$de
                law.or.fdp[i, j]<-sum((1-thetaS.vec)*law.or.de)/max(sum(law.or.de), 1)
                law.or.ntp[i, j]<-sum(thetaS.vec*law.or.de)/sum(thetaS.vec)
                
                law.dd.res<-law.func(pvs=pv.vec, pis.hat, q)
                law.dd.de<-law.dd.res$de
                law.dd.fdp[i, j]<-sum((1-thetaS.vec)*law.dd.de)/max(sum(law.dd.de), 1)
                law.dd.ntp[i, j]<-sum(thetaS.vec*law.dd.de)/sum(thetaS.vec)
                
                
                # use python package smoothfdr
                #np_x = npy$array(x)
                #results = smooth_fdr(np_x, q, verbose=0, missing_val=0)
                #tansey.fdp[i, j] = sum((1-thetaS.vec)*results$discoveries)/max(sum(results$discoveries), 1)
                #tansey.ntp[i, j] = sum(thetaS.vec*results$discoveries)/sum(thetaS.vec)
                
                tansey.fdp[i, j] = 0
                tansey.ntp[i, j] = 0
                
                cbind(i,j,bh.fdp[i, j],law.or.fdp[i, j],law.dd.fdp[i, j],tansey.fdp[i, j],bh.ntp[i, j],law.or.ntp[i, j],law.dd.ntp[i, j],tansey.ntp[i, j])
              })

mean_out1<-aggregate(for_out1,by=list(for_out1[,1]),FUN='mean')
fdr1<-mean_out1[,4:7]
etp1<-mean_out1[,8:11]

save(mu0.vec,fdr1,etp1,for_out1,file="rect_mu.RData")


################################################################################
##### different sparsity levels
################################################################################

isA.func<-function(i, j)
{
  # decide whether (i, j) is in triangle A
  if((j <= 1.25*i+75) && (j >= 100) && j <= -1.25*i+225)
  {y<-1}
  else if((j <= 1.25*i-25) && (j >= 100) && j <= -1.25*i+325)
  {y<-1}
  else
  {y<-0}
  return(y)
}

isB.func<-function(i, j)
{
  # decide whether (i, j) is in the rectangle B
  if(i>=20 && i<=100 && j>=40 && j<=90)
    y<-1
  else
    y<-0
  return(y)
}

pi0.vec<-seq(from=0.6, to=0.9, by=0.05)
pis<-matrix(0, dims[1], dims[2])
np<-length(pi0.vec)
mu0<-3
loop_matrix<-expand.grid(1:np, 1:nrep)

# Methods

bh.fdr2<-rep(0, np)
bh.etp2<-rep(0, np)

law.or.fdr2<-rep(0, np)
law.or.etp2<-rep(0, np)

law.dd.fdr2<-rep(0, np)
law.dd.etp2<-rep(0, np)

tansey.fdr2<-rep(0, np)
tansey.etp2<-rep(0, np)

bh.fdp2<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp2<-matrix(rep(0, nrep*np), np, nrep)

law.or.fdp2<-matrix(rep(0, nrep*np), np, nrep)
law.or.ntp2<-matrix(rep(0, nrep*np), np, nrep)

law.dd.fdp2<-matrix(rep(0, nrep*np), np, nrep)
law.dd.ntp2<-matrix(rep(0, nrep*np), np, nrep)

tansey.fdp2<-matrix(rep(0, nrep*np), np, nrep)
tansey.ntp2<-matrix(rep(0, nrep*np), np, nrep)


system.time(for_out2 <- foreach(iter = 1:(np*nrep), .combine='rbind') %dopar%
              {
                i<-as.numeric(loop_matrix[iter,][1])
                j<-as.numeric(loop_matrix[iter,][2])
                
                cat("\n", "iteration i = ", i, "\n", "iteration j = ", j, "\n")
                pi0<-pi0.vec[i]
                
                set.seed(iter+2000)
                for(i1 in 1:dims[1])
                {
                  for (j1 in 1:dims[2])
                  {
                    if (isB.func(i1,j1))
                    {pis[i1,j1] = pi0
                    thetaS[i1,j1]<-rbinom(1, 1, pi0)}
                    else
                    {pis[i1,j1] = 0.01
                    thetaS[i1,j1]<-rbinom(1, 1, 0.01)}
                  }
                }
                
                thetaS.vec<-c(thetaS)
                pis.vec = c(pis)
                
                pii<-sum(thetaS.vec)/m
                x0<-rnorm(m, mean=0, sd=1)
                x1<-rnorm(m, mean=mu0, sd=1)
                x.vec<-(1-thetaS.vec)*x0+thetaS.vec*x1
                
                x<-matrix(x.vec, dims[1], dims[2])
                
                pv.vec<-2*pnorm(-abs(x.vec), 0, 1)
                bh.th<-bh.func(pv.vec, 0.9)$th
                pis.hat<-pis_2D.func(x, tau=bh.th, h=10)
                
                
                bh.res<-bh.func(pv.vec, q)
                bh.de<-bh.res$de
                bh.fdp2[i, j]<-sum((1-thetaS.vec)*bh.de)/max(sum(bh.de), 1)
                bh.ntp2[i, j]<-sum(thetaS.vec*bh.de)/sum(thetaS.vec)
                
                law.or.res<-law.func(pvs=pv.vec, pis.vec, q)
                law.or.de<-law.or.res$de
                law.or.fdp2[i, j]<-sum((1-thetaS.vec)*law.or.de)/max(sum(law.or.de), 1)
                law.or.ntp2[i, j]<-sum(thetaS.vec*law.or.de)/sum(thetaS.vec)
                
                law.dd.res<-law.func(pvs=pv.vec, pis.hat, q)
                law.dd.de<-law.dd.res$de
                law.dd.fdp2[i, j]<-sum((1-thetaS.vec)*law.dd.de)/max(sum(law.dd.de), 1)
                law.dd.ntp2[i, j]<-sum(thetaS.vec*law.dd.de)/sum(thetaS.vec)
                
                
                # use python package smoothfdr
                #np_x = npy$array(x)
                #results = smooth_fdr(np_x, q, verbose=0, missing_val=0)
                #tansey.fdp2[i, j] = sum((1-thetaS.vec)*results$discoveries)/max(sum(results$discoveries), 1)
                #tansey.ntp2[i, j] = sum(thetaS.vec*results$discoveries)/sum(thetaS.vec)
                
                tansey.fdp2[i, j] = 0
                tansey.ntp2[i, j] = 0
                
                cbind(i,j,bh.fdp2[i, j],law.or.fdp2[i, j],law.dd.fdp2[i, j],tansey.fdp2[i, j],bh.ntp2[i, j],law.or.ntp2[i, j],law.dd.ntp2[i, j],tansey.ntp2[i, j])
              })

mean_out2<-aggregate(for_out2,by=list(for_out2[,1]),FUN='mean')
fdr2<-mean_out2[,4:7]
etp2<-mean_out2[,8:11]

save(pi0.vec,fdr2,etp2,for_out2,file="rect_spa.RData")

stopImplicitCluster()

