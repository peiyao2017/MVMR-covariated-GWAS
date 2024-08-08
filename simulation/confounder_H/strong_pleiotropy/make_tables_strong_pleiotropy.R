setwd("D:/paper2/simulation/strong_pleiotropy/")
filename1=numeric()
filename2=numeric()
for(i in 1:20){
  filename1[i]=paste("result",i,"confounder_Hdim=2_directional_correlated_pleiotropy30.RData",sep="")
  filename2[i]=paste("result",i,"confounder_Hdim=2_directional_correlated_pleiotropy50.RData",sep="")
}
for(i in 1:20){
  load(filename1[i])
  if(i==1){
    data1=output[[1]]
  }
  if(i>1){
    data1=rbind(data1,output[[1]])
  }
  load(filename2[i])
  if(i==1){
    data2=output[[1]]
  }
  if(i>1){
    data2=rbind(data2,output[[1]])
  }
}
 
tableU=data.frame(
  theta1=rep(0,times=10),
  sd1=rep(0,times=10),
  se1=rep(0,times=10),
  theta2=rep(0,times=10),
  sd2=rep(0,times=10),
  se2=rep(0,times=10)
)
tableU$theta1[1]=mean(data1$theta_uv_cml_X_and_Y_adj)
tableU$sd1[1]=sd(data1$theta_uv_cml_X_and_Y_adj)
tableU$se1[1]=mean(data1$se_theta_uv_cml_X_and_Y_adj)
tableU$theta1[2]=mean(data1$theta_mv_cml_X_and_Y_adj)
tableU$sd1[2]=sd(data1$theta_mv_cml_X_and_Y_adj)
tableU$se1[2]=mean(data1$se_theta_mv_cml_X_and_Y_adj)
tableU$theta1[3]=mean(data1$theta_uv_lasso_X_and_Y_adj)
tableU$sd1[3]=sd(data1$theta_uv_lasso_X_and_Y_adj)
tableU$se1[3]=mean(data1$se_theta_uv_lasso_X_and_Y_adj)
tableU$theta1[4]=mean(data1$theta_mv_lasso_X_and_Y_adj)
tableU$sd1[4]=sd(data1$theta_mv_lasso_X_and_Y_adj)
tableU$se1[4]=mean(data1$se_theta_mv_lasso_X_and_Y_adj)
tableU$theta1[5]=mean(data1$theta_uv_median_X_and_Y_adj)
tableU$sd1[5]=sd(data1$theta_uv_median_X_and_Y_adj)
tableU$se1[5]=mean(data1$se_theta_uv_median_X_and_Y_adj)
tableU$theta1[6]=mean(data1$theta_mv_median_X_and_Y_adj)
tableU$sd1[6]=sd(data1$theta_mv_median_X_and_Y_adj)
tableU$se1[6]=mean(data1$se_theta_mv_median_X_and_Y_adj)
tableU$theta1[7]=mean(data1$theta_uv_egger_X_and_Y_adj)
tableU$sd1[7]=sd(data1$theta_uv_egger_X_and_Y_adj)
tableU$se1[7]=mean(data1$se_theta_uv_egger_X_and_Y_adj)
tableU$theta1[8]=mean(data1$theta_mv_egger_X_and_Y_adj)
tableU$sd1[8]=sd(data1$theta_mv_egger_X_and_Y_adj)
tableU$se1[8]=mean(data1$se_theta_mv_egger_X_and_Y_adj)
tableU$theta1[9]=mean(data1$theta_uv_ivw_X_and_Y_adj)
tableU$sd1[9]=sd(data1$theta_uv_ivw_X_and_Y_adj)
tableU$se1[9]=mean(data1$se_theta_uv_ivw_X_and_Y_adj)
tableU$theta1[10]=mean(data1$theta_mv_ivw_X_and_Y_adj)
tableU$sd1[10]=sd(data1$theta_mv_ivw_X_and_Y_adj)
tableU$se1[10]=mean(data1$se_theta_mv_ivw_X_and_Y_adj)
tableU$theta2[1]=mean(data2$theta_uv_cml_X_and_Y_adj)
tableU$sd2[1]=sd(data2$theta_uv_cml_X_and_Y_adj)
tableU$se2[1]=mean(data2$se_theta_uv_cml_X_and_Y_adj)
tableU$theta2[2]=mean(data2$theta_mv_cml_X_and_Y_adj)
tableU$sd2[2]=sd(data2$theta_mv_cml_X_and_Y_adj)
tableU$se2[2]=mean(data2$se_theta_mv_cml_X_and_Y_adj)
tableU$theta2[3]=mean(data2$theta_uv_lasso_X_and_Y_adj)
tableU$sd2[3]=sd(data2$theta_uv_lasso_X_and_Y_adj)
tableU$se2[3]=mean(data2$se_theta_uv_lasso_X_and_Y_adj)
tableU$theta2[4]=mean(data2$theta_mv_lasso_X_and_Y_adj)
tableU$sd2[4]=sd(data2$theta_mv_lasso_X_and_Y_adj)
tableU$se2[4]=mean(data2$se_theta_mv_lasso_X_and_Y_adj)
tableU$theta2[5]=mean(data2$theta_uv_median_X_and_Y_adj)
tableU$sd2[5]=sd(data2$theta_uv_median_X_and_Y_adj)
tableU$se2[5]=mean(data2$se_theta_uv_median_X_and_Y_adj)
tableU$theta2[6]=mean(data2$theta_mv_median_X_and_Y_adj)
tableU$sd2[6]=sd(data2$theta_mv_median_X_and_Y_adj)
tableU$se2[6]=mean(data2$se_theta_mv_median_X_and_Y_adj)
tableU$theta2[7]=mean(data2$theta_uv_egger_X_and_Y_adj)
tableU$sd2[7]=sd(data2$theta_uv_egger_X_and_Y_adj)
tableU$se2[7]=mean(data2$se_theta_uv_egger_X_and_Y_adj)
tableU$theta2[8]=mean(data2$theta_mv_egger_X_and_Y_adj)
tableU$sd2[8]=sd(data2$theta_mv_egger_X_and_Y_adj)
tableU$se2[8]=mean(data2$se_theta_mv_egger_X_and_Y_adj)
tableU$theta2[9]=mean(data2$theta_uv_ivw_X_and_Y_adj)
tableU$sd2[9]=sd(data2$theta_uv_ivw_X_and_Y_adj)
tableU$se2[9]=mean(data2$se_theta_uv_ivw_X_and_Y_adj)
tableU$theta2[10]=mean(data2$theta_mv_ivw_X_and_Y_adj)
tableU$sd2[10]=sd(data2$theta_mv_ivw_X_and_Y_adj)
tableU$se2[10]=mean(data2$se_theta_mv_ivw_X_and_Y_adj)


a=matrix(c("UVMR-cML&","MVMR-cML&","UVMR-Lasso&","MVMR-Lasso&","UVMR-Median&","MVMR-Median&","UVMR-Egger&","MVMR-Egger&","UVMR-IVW&","MVMR-IVW&"),ncol=1)

for(i in 1:nrow(a)){
  for(j in 1:ncol(tableU)){
    if(j %in%c(2,3,5,6)&j<ncol(tableU)){
      a[i,]=paste(a[i,],"$(",formatC(tableU[i,j],format = "f",digits = 2),")$&")
    }
    if(j %in%c(2,3,5,6)&j==ncol(tableU)){
      a[i,]=paste(a[i,],"$(",formatC(tableU[i,j],format = "f",digits = 2),")$\\")
    }
    if((!j %in%c(2,3,5,6))&j<ncol(tableU)){
      a[i,]=paste(a[i,],"$",formatC(tableU[i,j],format = "f",digits = 2),"$&")
    }
    if((!j %in%c(2,3,5,6))&j==ncol(tableU)){
      a[i,]=paste(a[i,],"$",formatC(tableU[i,j],format = "f",digits = 2),"$\\")
    }
  }
}


prmatrix(a,rowlab = rep(" ",times=100),collab = rep(" ",times=100),quote=FALSE)

tableV=data.frame(
   no1=rep(0,times=5),
   sdno1=rep(0,times=5),
   seno1=rep(0,times=5),
   cml1=rep(0,times=5),
   sdcml1=rep(0,times=5),
   secml1=rep(0,times=5),
   ivw1=rep(0,times=5),
   sdivw1=rep(0,times=5),
   seivw1=rep(0,times=5),
   no2=rep(0,times=5),
   sdno2=rep(0,times=5),
   seno2=rep(0,times=5),
   cml2=rep(0,times=5),
   sdcml2=rep(0,times=5),
   secml2=rep(0,times=5),
   ivw2=rep(0,times=5),
   sdivw2=rep(0,times=5),
   seivw2=rep(0,times=5)
)
tableV$no1[1]=mean(data1$theta_uv_cml_X_and_Y_adj)
tableV$sdno1[1]=sd(data1$theta_uv_cml_X_and_Y_adj)
tableV$seno1[1]=mean(data1$se_theta_uv_cml_X_and_Y_adj)
tableV$no1[2]=mean(data1$theta_uv_lasso_X_and_Y_adj)
tableV$sdno1[2]=sd(data1$theta_uv_lasso_X_and_Y_adj)
tableV$seno1[2]=mean(data1$se_theta_uv_lasso_X_and_Y_adj)
tableV$no1[3]=mean(data1$theta_uv_median_X_and_Y_adj)
tableV$sdno1[3]=sd(data1$theta_uv_median_X_and_Y_adj)
tableV$seno1[3]=mean(data1$se_theta_uv_median_X_and_Y_adj)
tableV$no1[4]=mean(data1$theta_uv_egger_X_and_Y_adj)
tableV$sdno1[4]=sd(data1$theta_uv_egger_X_and_Y_adj)
tableV$seno1[4]=mean(data1$se_theta_uv_egger_X_and_Y_adj)
tableV$no1[5]=mean(data1$theta_uv_ivw_X_and_Y_adj)
tableV$sdno1[5]=sd(data1$theta_uv_ivw_X_and_Y_adj)
tableV$seno1[5]=mean(data1$se_theta_uv_ivw_X_and_Y_adj)
tableV$cml1[1]=mean(data1$theta_corrected_cml_cml)
tableV$sdcml1[1]=sd(data1$theta_corrected_cml_cml)
tableV$secml1[1]=mean(data1$theta_corrected_cml_cml_se)
tableV$cml1[2]=mean(data1$theta_corrected_lasso_cml)
tableV$sdcml1[2]=sd(data1$theta_corrected_lasso_cml)
tableV$secml1[2]=mean(data1$theta_corrected_lasso_cml_se)
tableV$cml1[3]=mean(data1$theta_corrected_median_cml)
tableV$sdcml1[3]=sd(data1$theta_corrected_median_cml)
tableV$secml1[3]=mean(data1$theta_corrected_median_cml_se)
tableV$cml1[4]=mean(data1$theta_corrected_egger_cml)
tableV$sdcml1[4]=sd(data1$theta_corrected_egger_cml)
tableV$secml1[4]=mean(data1$theta_corrected_egger_cml_se)
tableV$cml1[5]=mean(data1$theta_corrected_ivw_cml)
tableV$sdcml1[5]=sd(data1$theta_corrected_ivw_cml)
tableV$secml1[5]=mean(data1$theta_corrected_ivw_cml_se)
tableV$ivw1[1]=mean(data1$theta_corrected_cml_ivw)
tableV$sdivw1[1]=sd(data1$theta_corrected_cml_ivw)
tableV$seivw1[1]=mean(data1$theta_corrected_cml_ivw_se)
tableV$ivw1[2]=mean(data1$theta_corrected_lasso_ivw)
tableV$sdivw1[2]=sd(data1$theta_corrected_lasso_ivw)
tableV$seivw1[2]=mean(data1$theta_corrected_lasso_ivw_se)
tableV$ivw1[3]=mean(data1$theta_corrected_median_ivw)
tableV$sdivw1[3]=sd(data1$theta_corrected_median_ivw)
tableV$seivw1[3]=mean(data1$theta_corrected_median_ivw_se)
tableV$ivw1[4]=mean(data1$theta_corrected_egger_ivw)
tableV$sdivw1[4]=sd(data1$theta_corrected_egger_ivw)
tableV$seivw1[4]=mean(data1$theta_corrected_egger_ivw_se)
tableV$ivw1[5]=mean(data1$theta_corrected_ivw_ivw)
tableV$sdivw1[5]=sd(data1$theta_corrected_ivw_ivw)
tableV$seivw1[5]=mean(data1$theta_corrected_ivw_ivw_se)
tableV$no2[1]=mean(data2$theta_uv_cml_X_and_Y_adj)
tableV$sdno2[1]=sd(data2$theta_uv_cml_X_and_Y_adj)
tableV$seno2[1]=mean(data2$se_theta_uv_cml_X_and_Y_adj)
tableV$no2[2]=mean(data2$theta_uv_lasso_X_and_Y_adj)
tableV$sdno2[2]=sd(data2$theta_uv_lasso_X_and_Y_adj)
tableV$seno2[2]=mean(data2$se_theta_uv_lasso_X_and_Y_adj)
tableV$no2[3]=mean(data2$theta_uv_median_X_and_Y_adj)
tableV$sdno2[3]=sd(data2$theta_uv_median_X_and_Y_adj)
tableV$seno2[3]=mean(data2$se_theta_uv_median_X_and_Y_adj)
tableV$no2[4]=mean(data2$theta_uv_egger_X_and_Y_adj)
tableV$sdno2[4]=sd(data2$theta_uv_egger_X_and_Y_adj)
tableV$seno2[4]=mean(data2$se_theta_uv_egger_X_and_Y_adj)
tableV$no2[5]=mean(data2$theta_uv_ivw_X_and_Y_adj)
tableV$sdno2[5]=sd(data2$theta_uv_ivw_X_and_Y_adj)
tableV$seno2[5]=mean(data2$se_theta_uv_ivw_X_and_Y_adj)
tableV$cml2[1]=mean(data2$theta_corrected_cml_cml)
tableV$sdcml2[1]=sd(data2$theta_corrected_cml_cml)
tableV$secml2[1]=mean(data2$theta_corrected_cml_cml_se)
tableV$cml2[2]=mean(data2$theta_corrected_lasso_cml)
tableV$sdcml2[2]=sd(data2$theta_corrected_lasso_cml)
tableV$secml2[2]=mean(data2$theta_corrected_lasso_cml_se)
tableV$cml2[3]=mean(data2$theta_corrected_median_cml)
tableV$sdcml2[3]=sd(data2$theta_corrected_median_cml)
tableV$secml2[3]=mean(data2$theta_corrected_median_cml_se)
tableV$cml2[4]=mean(data2$theta_corrected_egger_cml)
tableV$sdcml2[4]=sd(data2$theta_corrected_egger_cml)
tableV$secml2[4]=mean(data2$theta_corrected_egger_cml_se)
tableV$cml2[5]=mean(data2$theta_corrected_ivw_cml)
tableV$sdcml2[5]=sd(data2$theta_corrected_ivw_cml)
tableV$secml2[5]=mean(data2$theta_corrected_ivw_cml_se)
tableV$ivw2[1]=mean(data2$theta_corrected_cml_ivw)
tableV$sdivw2[1]=sd(data2$theta_corrected_cml_ivw)
tableV$seivw2[1]=mean(data2$theta_corrected_cml_ivw_se)
tableV$ivw2[2]=mean(data2$theta_corrected_lasso_ivw)
tableV$sdivw2[2]=sd(data2$theta_corrected_lasso_ivw)
tableV$seivw2[2]=mean(data2$theta_corrected_lasso_ivw_se)
tableV$ivw2[3]=mean(data2$theta_corrected_median_ivw)
tableV$sdivw2[3]=sd(data2$theta_corrected_median_ivw)
tableV$seivw2[3]=mean(data2$theta_corrected_median_ivw_se)
tableV$ivw2[4]=mean(data2$theta_corrected_egger_ivw)
tableV$sdivw2[4]=sd(data2$theta_corrected_egger_ivw)
tableV$seivw2[4]=mean(data2$theta_corrected_egger_ivw_se)
tableV$ivw2[5]=mean(data2$theta_corrected_ivw_ivw)
tableV$sdivw2[5]=sd(data2$theta_corrected_ivw_ivw)
tableV$seivw2[5]=mean(data2$theta_corrected_ivw_ivw_se)

a1=matrix(c("UVMR-cML&","UVMR-Lasso&","UVMR-Median&","UVMR-Egger&","UVMR-IVW&"),ncol=1)
b1=c(2,3,5,6,8,9,11,12,14,15,17,18)
for(i in 1:nrow(a1)){
  for(j in 1:ncol(tableV)){
    if(j %in%b1&j<ncol(tableV)){
      a1[i,]=paste(a1[i,],"$(",formatC(tableV[i,j],format = "f",digits = 2),")$&")
    }
    if(j %in%b1&j==ncol(tableV)){
      a1[i,]=paste(a1[i,],"$(",formatC(tableV[i,j],format = "f",digits = 2),")$\\")
    }
    if((!j %in%b1)&j<ncol(tableV)){
      a1[i,]=paste(a1[i,],"$",formatC(tableV[i,j],format = "f",digits = 2),"$&")
    }
    if((!j %in%b1)&j==ncol(tableV)){
      a1[i,]=paste(a1[i,],"$",formatC(tableV[i,j],format = "f",digits = 2),"$\\")
    }
  }
}


prmatrix(a1,rowlab = rep(" ",times=100),collab = rep(" ",times=100),quote=FALSE)