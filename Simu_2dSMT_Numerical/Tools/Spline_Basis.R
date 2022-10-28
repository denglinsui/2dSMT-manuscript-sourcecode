##' Spline Generation
library(fda)
library(dplyr)
library(plotly)

#=== Init point
point <- matrix(c(rep(seq(0,1,length.out=m.row),each=m.row),
                  rep(seq(0,1,length.out=m.row),times=m.row)), 
                m, 2, byrow = F) 
colnames(point) <- c("x","y")

#=== Init Basis
n.break <- 3
n.order <- 4
n.end <- n.order + n.break -2
sp.basis <- create.bspline.basis(c(0.25,0.75), dropind=c(1,n.end),
                                 breaks = seq(0.25, 0.75, length.out = n.break), 
                                 norder = n.order)

point.basis <- point %>%  as.data.frame() %>% mutate(mu.x = 0) %>% mutate(mu.y = 0)
colnames(point.basis) <- c("x","y","mu.x","mu.y")

#=== Basis product
point.basis[point.basis$x>=0.25 & point.basis$x<=0.75,] <- 
  point.basis %>%
  filter(x>=0.25 & x<=0.75) %>% 
  mutate(mu.x = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(x,basisobj = sp.basis)))

point.basis[point.basis$y>=0.25 & point.basis$y<=0.75,] <- 
  point.basis %>%
  filter(y>=0.25 & y<=0.75) %>% 
  mutate(mu.y = sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(y,basisobj = sp.basis)))
#sqrt(n.end)/sqrt(n.end-2) is for normalization
point.basis <- point.basis %>% mutate(mu = mu.x*mu.y) %>% mutate(mu = mu)

p <- ggplot(point.basis,aes(x=x, y = y, z = mu, color =stat(level)))+
  geom_contour()+
  scale_colour_distiller(palette = "YlGn", direction = 1)
ggplotly(p)

#=== Exploration
point.basis %>% 
  dplyr::select(mu) %>% 
  dplyr::filter(mu>0) %>% 
  unlist %>%
  mean()



if(F){ # Generate Spline by fitting
  sp1 <- bs(point[,1],df = 20, intercept = F, Boundary.knots = c(0.25,0.75))
  sp2 <- bs(point[,2],df = 20, intercept = F, Boundary.knots = c(0.25,0.75))
  fit <- lm(mu ~ sp1 + sp2 -1)
  
  data <- data.frame(x=point[,1], y=point[,2], mu = mu,
                     X=X, pluszero = pluszero, mu.smooth = fit$fitted.values)
  
  p <- ggplot()+
    geom_point(data = data,aes(x=x, y = X, color = pluszero))+
    geom_line(data=data,aes(x=x, y=mu.smooth))
  
  p <- ggplot(data,aes(x=x, y = y, z = mu.smooth, color =stat(level)))+
    geom_contour()+
    scale_colour_distiller(palette = "YlGn", direction = 1)
  ggplotly(p)
  
}