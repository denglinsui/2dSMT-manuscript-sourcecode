library(reshape2)
library(dplyr)
library(ggplot2)
Files = 
  c("Result/Simulation2D_GMLE/GMLE_Gm_tGm/mag_0.5 mu_Dense Cov_Strong 2024-03-20 19-23-15.RData",
    "Result/Simulation2D_GMLE/GMLE_Gm_tGm/mag_1 mu_Dense Cov_Strong 2024-03-20 19-00-33.RData",
    "Result/Simulation2D_GMLE/GMLE_Gm_tGm/mag_1.5 mu_Dense Cov_Strong 2024-03-20 19-02-36.RData")
Err_all <- NULL
for(files in Files){
  load(files)
  Err_df =  melt(Err_res,id = c("mSub","m"))
  Err_all = rbind(Err_all,
                  data.frame(Err_df,
                             mag=magnitude))
  # Err_df = Err_df%>% group_by(m,mSub,variable) %>%
  #   summarise(Val =  mean(value))
  # p = ggplot(data=Err_df)+
  #   geom_smooth(aes(x=m,y=value,color = variable,fill=variable),method = "loess")+
  #   labs(title = paste0(Cov_type,mu_type,magnitude))
  # p = ggplot(data=Err_df)+
  #   geom_line(aes(x=m,y=Val,color = variable))+
  #   labs(title = paste0(Cov_type,mu_type,magnitude))
  #print(p)
}
library(latex2exp)
Err_all$mag = paste0("$\\gamma=",Err_all$mag,"$")
Err_all <- Err_all %>% 
  filter(variable != "EmtG.GMLEsub") 
Err_all$Type =  as.character(Err_all$variable)
Err_all$Type[Err_all$Type=="EmG.EmtG"]=
  r"($d_H(f_{G_{\tilde{m}}})$)"
Err_all$Type[Err_all$Type=="EmG.GMLEfull"]=
  r"($d_H(f_{G_m},f_{\hat{G}_m})$)"
Err_all$Type[Err_all$Type=="EmpG.GMLEsub"]=
  r"($d_H(f_{G_m},f_{\hat{G}_{\tilde{m}}})$)"
Err_all.sub = Err_all[c(1:10,1301:1310,3601:3610),]
ggplot(data=Err_all.sub,
       aes(x=m,y=value,color = Type))+
  scale_color_manual(values=c( r"($d_H(f_{G_{\tilde{m}}})$)"="blue",
                              r"($d_H(f_{G_m},f_{\hat{G}_{\tilde{m}}})$)"="red",
                              r"($d_H(f_{G_m},f_{\hat{G}_m})$)"="black"),
                     labels = TeX(Err_all.sub$Type)) +
                       #TeX(Err_all.sub$Type)) +
  geom_point() 
p = ggplot(data=Err_all)+
  geom_smooth(aes(x=m,y=value,color = Type,fill=Type),
              method = "loess")+
  #facet_grid(~mag,parser=T)+
  facet_grid(~ TeX(mag, output = "character"), 
             labeller = label_parsed)+
  theme_bw()+
  scale_color_manual(values=c( r"($d_H(f_{G_{\tilde{m}}})$)"="#1F78B4",
                              r"($d_H(f_{G_m},f_{\hat{G}_{\tilde{m}}})$)"="#33A02C",
                              r"($d_H(f_{G_m},f_{\hat{G}_m})$)"="#E31A1C"),
                     labels = TeX(Err_all$Type)) +
  scale_fill_manual(values=c(  r"($d_H(f_{G_{\tilde{m}}})$)"="#1F78B4",
                              r"($d_H(f_{G_m},f_{\hat{G}_{\tilde{m}}})$)"="#33A02C",
                              r"($d_H(f_{G_m},f_{\hat{G}_m})$)"="#E31A1C"),
                     labels = TeX(Err_all$Type)) +
  theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.text = element_text(size=12),
    legend.position = "bottom")+
  ylab("Hellinger Distance")

ggsave("Figure/Est_NPEB.pdf",plot = p, width = 8,height =4)

Err_df = Err_all%>% group_by(m,mSub,variable) %>%
  summarise(Val =  mean(value))
ggplot(Err_df)+
  geom_line(data = Err_df %>% filter(variable=="EmG.GMLEfull"),
            aes(x=m,y=Val))+
  geom_line(data = Err_df %>% filter(variable=="EmpG.GMLEsub"),
            aes(x=mSub,y=Val))
  
