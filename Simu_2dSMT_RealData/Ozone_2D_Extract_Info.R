library(data.table)
library(dplyr)

DetectNum_ST <- NULL
DetectNum_IHW <- NULL
DetectNum_SA <- NULL
file_num <- 0
for(hh in 2){
  for(beta0 in seq(-0.1,-0.5,by=-0.1)){
    file_num <- file_num+1
    load(paste0("Fig/Ozone_h",hh,"_",beta0,".RData"))
    Tmp_Num_ST <- T2_Stat[order(State.Code),.(Null = sum(Detect.ST=="Null"),
                                              ST = sum(Detect.ST%in%c("ST","ST & 2D(ST)")),
                                              ST_2D = sum(Detect.ST%in%c("2D(ST)","ST & 2D(ST)")),
                                              NonEqual = any(Detect.ST== "ST"|Detect.ST=="2D(ST)"),
                                              Total = length(Detect.ST)),
                          by = .(State.Code)]
    colnames(Tmp_Num_ST) <- paste(paste0("h",hh,"_",beta0),c("State.Code","Null","ST","ST_2D","NonEqual","Total"))
    DetectNum_ST <- cbind(DetectNum_ST,Tmp_Num_ST)
    Tmp_Num_IHW <- T2_Stat[order(State.Code),.(Null = sum(Detect.IHW=="Null"),
                                               ST = sum(Detect.IHW%in%c("IHW","IHW & 2D(IHW)")),
                                               ST_2D = sum(Detect.IHW%in%c("2D(IHW)","IHW & 2D(IHW)")),
                                               NonEqual = any(Detect.IHW == "IHW"|Detect.IHW =="2D(IHW)"),
                                               Total = length(Detect.IHW)),
                           by = .(State.Code)]
    colnames(Tmp_Num_IHW) <- paste(paste0("h",hh,"_",beta0),c("State.Code","Null","IHW","IHW_2D","NonEqual","Total"))
    DetectNum_IHW <- cbind(DetectNum_IHW,Tmp_Num_IHW)
    
    Tmp_Num_SA <- T2_Stat[order(State.Code),.(Null = sum(Detect.SA=="Null"),
                                              SA = sum(Detect.SA%in%c("SA","SA & 2D(SA)")),
                                              SA_2D = sum(Detect.SA%in%c("2D(SA)","SA & 2D(SA)")),
                                              NonEqual = any(Detect.SA== "SA"|Detect.SA=="2D(SA)"),
                                              Total = length(Detect.SA)),
                          by = .(State.Code)]
    colnames(Tmp_Num_SA) <- paste(paste0("h",hh,"_",beta0),c("SAate.Code","Null","SA","SA_2D","NonEqual","Total"))
    DetectNum_SA <- cbind(DetectNum_SA,Tmp_Num_SA)
  }
}

OneD_col <- 1:file_num*6-3
TwoD_col <- 1:file_num*6-2
NonEqu_col <- 1:file_num*6-1
MD_col <- rep(1:file_num*6,each=2)-c(3,2)
diff_row.ST <- apply(DetectNum_ST[,..NonEqu_col],1,any)
DetectNum_ST[diff_row.ST,..MD_col]


diff_row.IHW <- apply(DetectNum_IHW[,..NonEqu_col],1,any)
DetectNum_IHW[diff_row.IHW,]
DetectNum_IHW[diff_row.IHW,..MD_col]

diff_row.SA <- apply(DetectNum_SA[,..NonEqu_col],1,any)
DetectNum_SA[diff_row.SA,..MD_col]



DetectNum_ST[diff_row.ST|diff_row.IHW|diff_row.SA,]
DetectNum_ST[diff_row.ST|diff_row.IHW||diff_row.SA,..MD_col]
DetectNum_IHW[diff_row.ST|diff_row.IHW|diff_row.SA,,..MD_col]
