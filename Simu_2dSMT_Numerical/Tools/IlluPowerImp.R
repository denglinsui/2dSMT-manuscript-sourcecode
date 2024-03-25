library(Rsolnp)
library(latex2exp)
library(pbivnorm)
set.seed(10)
rho.seq =  seq(0,0.5,length.out= 10)
kappa.seq = 2:7
pi0.seq =  seq(0.9,0.99,length.out= 4)
mu1.seq  = seq(2,5,length.out = 100)

para.df = expand.grid(rho.seq,kappa.seq,pi0.seq,mu1.seq)
pow.df = NULL
for(ind in 1:dim(para.df)[1]){
  print(ind)
  para.seq = as.numeric(para.df[ind,])
  tmp =pow_imp(rho = para.seq[1],
               kappa = para.seq[2],
               pi0=para.seq[3],
              mu1=para.seq[4])
  pow.df = rbind(pow.df,tmp)
}
big.df = cbind(para.df,pow.df)
names(big.df) = c("rho","kappa","pi0","mu1","pow1d","pow2d")
big.df$pi0anno = paste0("$\\pi_0=$",big.df$pi0)
big.df$rhoanno = paste0("$\\rho_X=$",big.df$rho)
big.df = as.data.frame(big.df)
library(ggplot2)
library(dplyr)
ggplot(big.df%>%filter(kappa == 4),
       aes(x = pow1d,y=(pow2d-pow1d)/pow1d,
           color=rho))+
 # facet_wrap(~pi0,nrow = 1,scales = "free")+
  facet_wrap(.~TeX(pi0anno, output = "character"), 
             nrow = 1,scales = "free",
             labeller = label_parsed)+
  scale_colour_gradient(high = "#1F78B4",
                        low ="#E31A1C")+
  #geom_point()+
  geom_line(aes(group = rho))+
  guides(color = guide_colourbar(barwidth = 1, 
                                 barheight = 10,
                                 reverse=TRUE))+
  theme_bw()+
  labs(x = TeX("$PTD^{1d}$"),
       y = TeX("$(PTD^{2d}-PTD^{1d})/PTD^{1d}$"),
       color= TeX("$\\rho_X$"))+
  theme_here()
ggsave("Figure/kappa4.pdf",width = 13,height=4)

ggplot(big.df%>%filter(rho == 0.5|rho == 0),
       aes(x = pow1d,y=(pow2d-pow1d)/pow1d,
           color=kappa))+
  # facet_wrap(~pi0anno,nrow = 1,scales = "free")+
  facet_wrap(factor(TeX(rhoanno, output = "character"),
                    levels = c("rho[X] * {phantom() == phantom()} * '0'",
                               "rho[X] * {phantom() == phantom()} * '0.5'"))~TeX(pi0anno, output = "character"), 
             nrow = 2,
             scales = "free",
             labeller = label_parsed)+
  #geom_point()+
  geom_line(aes(group = kappa))+
  #geom_line(aes(group ))+
  scale_colour_gradient(low = "#1F78B4",
                        high ="#E31A1C")+
  guides(color = guide_colourbar(barwidth = 1, 
                                 barheight = 10))+
  theme_bw()+
  labs(x = TeX("$PTD^{1d}$"),
       y = TeX("$(PTD^{2d}-PTD^{1d})/PTD^{1d}$"),
       color= TeX("$\\kappa$"))+
  theme_here()

ggsave("Figure/rho0_05.pdf",width = 13,height=7)

#ggsave("Figure/rho0.pdf",width = 13,height=4)


ggplot(big.df%>%filter(kappa == 4 & pow1d>0.2&rho%in%c(0,0.5)),
       aes(x = pow1d,y=(pow2d-pow1d)/pow1d,
           color=pi0))+
  facet_wrap(.~factor(TeX(rhoanno, output = "character"),
                      levels = c("rho[X] * {phantom() == phantom()} * '0'",
                                 "rho[X] * {phantom() == phantom()} * '0.5'")), 
             nrow = 1,scales = "free",
             labeller = label_parsed)+
  
  geom_line(aes(group = pi0))+
  scale_colour_gradient(low = "#1F78B4",
                        high ="#E31A1C")+
  guides(color = guide_colourbar(barwidth = 1, 
                                 barheight = 10))+
  theme_bw()+
  labs(x = TeX("$PTD^{1d}$"),
       y = TeX("$(PTD^{2d}-PTD^{1d})/PTD^{1d}$"),
       color= TeX("$\\pi_0$"))+
  theme_here()

ggsave("Figure/kappa4pi0.pdf",width = 9,height=4)
pow_imp = function(rho,kappa,pi0,mu1){
Kval_2d = function(x){
  x20=x[1]
  y20=x[2]
  taus = sqrt(kappa +(kappa^2-kappa)*rho)
  mu_a = kappa*mu1/taus
  rho_s = kappa * rho /taus
  x2 = x20
  y2 = y20
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0 = (1-B2 -C2 + A3)
  x2 = x20 - mu1
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K1 = (1-B2 -C2 + A3)
  
  x2 = x20 
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0est2 = (1-B2 -C2 + A3)
  finalval =(pi0*K0+(1-pi0)*K1)
  return(-finalval)
}


K1val_2d = function(x){
  x20=x[1]
  y20=x[2]
  taus = sqrt(kappa +(kappa^2-kappa)*rho)
  mu_a = kappa*mu1/taus
  rho_s = kappa * rho /taus
  x2 = x20
  y2 = y20
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0 = (1-B2 -C2 + A3)
  x2 = x20 - mu1
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K1 = (1-B2 -C2 + A3)
  
  x2 = x20 
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0est2 = (1-B2 -C2 + A3)
  finalval =K1
  return(finalval)
}

FDPval_2d = function(x){
  x20=x[1]
  y20=x[2]
  taus = sqrt(kappa +(kappa^2-kappa)*rho)
  mu_a = kappa*mu1/taus
  rho_s = kappa * rho /taus
  x2 = x20
  y2 = y20
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0 = (1-B2 -C2 + A3)
  x2 = x20 - mu1
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K1 = (1-B2 -C2 + A3)
  
  x2 = x20 
  y2 = y20 - mu_a
  A3 <- pbivnorm(x = x2, y = y2, rho = rho_s)
  B2 <- pnorm(x2)
  C2 <- pnorm(y2)
  K0est2 = (1-B2 -C2 + A3)
  if(pi0*(pi0*K0 +(1-pi0)*K0est2)==0){
    return(0)
  }else{
  return((pi0*(pi0*K0 +(1-pi0)*K0est2))/(pi0*K0+(1-pi0)*K1))
  }
}



Kval_1d = function(x){
  x20=x[1]
  
  K1 = 1-pnorm(x20,mean = mu1)
  
  K0 = K0est2 =  1-pnorm(x20,mean = 0)
  finalval =(pi0*K0+(1-pi0)*K1)
  return(-finalval)
}


K1val_1d = function(x){
  x20=x[1]
  
  K1 = 1-pnorm(x20,mean = mu1)
  
  return(K1)
}

FDPval_1d = function(x){
  x20=x[1]
  
  K1 = 1-pnorm(x20,mean = mu1)
  
  K0 = K0est2 =  1-pnorm(x20)
  if(pi0*(pi0*K0 +(1-pi0)*K0est2)==0){
    return(0)}else{
      
      return((pi0*(pi0*K0 +(1-pi0)*K0est2))/(pi0*K0+(1-pi0)*K1))
    }
}

res = uniroot(function(x){FDPval_1d(x)-0.1}, lower = -5, upper = 20)
#FDPval_1d(res$root)
#Kval_1d(res$root)
pow1d = K1val_1d(res$root)

res = solnp(c(1,2),
      Kval_2d, 
      eqfun=FDPval_2d,
      eqB=0.1,
      LB=c(-10,-10), 
      UB=c(20,20)) 
#res$pars
#FDPval_2d(res$pars)
#K01val_2d(res$pars)
pow2d = K1val_2d(res$pars)


return(c(pow1d,pow2d))
}

theme_here = function(){
  
  return(
    # guides(colour=guide_legend(title = TeX("$\\rho_X$")))+
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      axis.title=element_text(size=15),
      axis.text=element_text(size=12),
      legend.text =  element_text(size = 15),
      legend.title =  element_text(size = 15)#,
      #legend.position = "bottom"
      ))
}

