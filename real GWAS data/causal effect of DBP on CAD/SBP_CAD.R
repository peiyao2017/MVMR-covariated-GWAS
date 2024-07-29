setwd("D:/paper2/real_data/SBP_CAD")
DP=100
library(data.table)
library(MendelianRandomization)
library(MVMR)
library(MASS)
library(glmnet)
cv.mvmr_lasso = function(bx, by, seby){
  p = dim(bx)[1]
  k = dim(bx)[2]
  S = diag(seby^-2)
  b = S^(1/2) %*% bx
  Pb = b %*% solve(t(b) %*% b, t(b))
  xlas = (diag(p) - Pb) %*% S^(1/2)
  ylas = (diag(p) - Pb) %*% S^(1/2) %*% by
  alas = glmnet(xlas, ylas, intercept = FALSE,lambda = seq(from=0,to=10,by=0.1))
  av=numeric()
  for(i in 1:length(alas$lambda)){
    av[i]=sum(alas$beta[, i] == 0)
  }
  lambda_total=alas$lambda[av>0]
  alas = glmnet(xlas, ylas, intercept = FALSE,lambda = lambda_total)
  lamseq = sort(alas$lambda)
  lamlen = length(lamseq)
  rse = sapply(1:lamlen, function(j){
    av = which(alas$beta[, (lamlen - j + 1)] == 0)
    
    mod = lm.fit(as.matrix(S[av, av]^(1/2) %*% bx[av, ]), S[av, av]^(1/2) %*% by[av])
    c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
    
  })
  rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
  het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
  if (length(het) == 0){
    lam_pos = 1
  } else {
    lam_pos = min(het)
  }
  num_valid = rev(sapply(1:lamlen, function(j){sum(alas$beta[, j]==0)}))
  min_lam_pos = min(which(num_valid > k))
  if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
  return(list(fit = alas$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos]))
}


mvmr_lasso = function(bx, by, seby){
  p = dim(as.matrix(bx))[1]
  k = dim(as.matrix(bx))[2]
  S = diag(seby^-2)
  sn = sign(bx[, 1])
  bx_or = bx * sn
  by_or = by * sn
  cv.alas = cv.mvmr_lasso(bx_or, by_or, seby)
  a1 = cv.alas$fit
  e = by_or - a1
  thest = solve(t(bx_or) %*% S %*% bx_or, t(bx_or) %*% S %*% e)
  v = which(a1==0)
  mvmr_mod = mr_mvivw(mr_mvinput(bx = as.matrix(bx_or[v, ],nrow=sum(a1==0),ncol=ncol(bx_or)), bxse = as.matrix(bx_or[v, ],nrow=sum(a1==0),ncol=ncol(bx_or)),
                                 by = by_or[v], byse = seby[v]))
  th_post = mvmr_mod$Estimate
  se_post = mvmr_mod$StdError
  return(list(thest = thest, a = a1, lambda = cv.alas$lambda,
              th_post = th_post, se_post = se_post))
}    


BMI=data.table(fread("D:/paper2/real_data/SBP_CAD/GWAS_BMI_harmo.txt",sep="\t",header = TRUE))
CAD=data.table(fread("D:/paper2/real_data/SBP_CAD/GWAS_CAD_harmo.txt",sep="\t",header = TRUE))
SBP=data.table(fread("D:/paper2/real_data/SBP_CAD/GWAS_SBP_harmo.txt",sep="\t",header = TRUE))
DBP=data.table(fread("D:/paper2/real_data/SBP_CAD/GWAS_DBP_harmo.txt",sep="\t",header = TRUE))
load("D:/paper2/real_data/SBP_CAD/ind_snps.RData")
ind_snp=ind_snps[[1]]
 
IV_SBP=SBP$SNP[SBP$SNP%in%ind_snp&SBP$p<5e-8]
IV_SBP=IV_SBP[!CAD$p[CAD$SNP%in%IV_SBP]<5e-8]
IV_BMI=BMI$SNP[BMI$SNP%in%ind_snp&BMI$p<5e-8]
IV_BMI=IV_BMI[!CAD$p[CAD$SNP%in%IV_BMI]<5e-8]
IV_DBP=DBP$SNP[DBP$SNP%in%ind_snp&DBP$p<5e-8]
IV_DBP=IV_DBP[!CAD$p[CAD$SNP%in%IV_DBP]<5e-8]
IV=unique(c(IV_DBP,IV_BMI,IV_SBP))
NIV=length(IV)
IV_beta_exp=cbind(DBP$beta[DBP$SNP%in%IV],BMI$beta[BMI$SNP%in%IV],SBP$beta[SBP$SNP%in%IV])
IV_se_exp=cbind(DBP$se[DBP$SNP%in%IV],BMI$se[BMI$SNP%in%IV],SBP$se[SBP$SNP%in%IV])
IV_beta_out=CAD$beta[CAD$SNP%in%IV]
IV_se_out=CAD$se[CAD$SNP%in%IV]
null=(SBP$p>0.1&BMI$p>0.1&DBP$p>0.1&CAD$p>0.1)&(SBP$SNP%in%ind_snp&BMI$SNP%in%ind_snp&DBP$SNP%in%ind_snp&CAD$SNP%in%ind_snp)
null_z=cbind(DBP$z[null],BMI$z[null],SBP$z[null],CAD$z[null])
correlation=matrix(0,nrow=4,ncol=4)
for(i in 1:4){
  for(j in 1:4){
    correlation[i,j]=cor(null_z[,i],null_z[,j])
  }
}
SIG=list()
SIG_inv=list()
SIG_exp=list()
for(i in 1:NIV){
  SIG[[i]]=diag(c(IV_se_exp[i,],IV_se_out[i]),nrow=4,ncol=4)%*%correlation%*%diag(c(IV_se_exp[i,],IV_se_out[i]),nrow=4,ncol=4)
  SIG_inv[[i]]=solve(SIG[[i]])
  SIG_exp[[i]]=SIG[[i]][1:3,1:3]
}
theta_cml_DP=matrix(0,ncol=3,nrow=DP)
theta_egger_DP=matrix(0,ncol=3,nrow=DP)
theta_lasso_DP=matrix(0,ncol=3,nrow=DP)
theta_median_DP=matrix(0,ncol=3,nrow=DP)
theta_ivw_DP=matrix(0,ncol=3,nrow=DP)
F_DP=matrix(0,nrow=DP,ncol=3)


for(i in 1:DP){
  beta_exp_DP=matrix(0,nrow=NIV,ncol=3)
  beta_out_DP=numeric()
  for(j in 1:NIV){
    beta_DP=mvrnorm(n=1,mu=c(IV_beta_exp[j,],IV_beta_out[j]),Sigma = SIG[[j]])
    beta_exp_DP[j,]=beta_DP[1:3]
    beta_out_DP[j]=beta_DP[4]
    if(beta_exp_DP[j,1]<0){
      beta_exp_DP[j,]=-beta_exp_DP[j,]
      beta_out_DP[j]=-beta_out_DP[j]
    }
  }
  object1=mr_mvinput(bx=beta_exp_DP,bxse = IV_se_exp,by=beta_out_DP,byse = IV_se_out)
  object2=format_mvmr(BXGs =beta_exp_DP,BYG=beta_out_DP,seBYG = IV_se_out,seBXGs = IV_se_exp,RSID=c(1:NIV) )
  
  
  cml_out=MVmr_cML(b_exp = beta_exp_DP,b_out = as.matrix(beta_out_DP,nrow=NIV,ncol=1),se_bx=IV_se_exp,Sig_inv_l=SIG_inv,n=339224,K_vec=c(1:(NIV-4)),min_theta_range=0,max_theta_range=0,thres=0.001,maxit=500  )
  ivw_out=mr_mvivw(object = object1)
  lasso_out=mvmr_lasso(bx=beta_exp_DP,by=beta_out_DP,seby = IV_se_out)
  egger_out=mr_mvegger(object = object1)
  median_out=mr_mvmedian(object = object1,iterations = 100)
  theta_cml_DP[i,]=cml_out$BIC_theta 
  theta_ivw_DP[i,]=ivw_out$Estimate 
  theta_lasso_DP[i,]=lasso_out$th_post 
  theta_egger_DP[i,]=egger_out$Estimate  
  theta_median_DP[i,]=median_out$Estimate 
  weak=strength_mvmr(r_input = object2,gencov = SIG_exp)[1,]
  for(j in 1:3){
    F_DP[i,j]=weak[1,j]
  }
  
}
DBP_CAD_all=list(data.frame(
  theta_cml=colMeans(theta_cml_DP),
  se_theta_cml=sqrt(diag(cov(theta_cml_DP))),
  theta_egger=colMeans(theta_egger_DP),
  se_theta_egger=sqrt(diag(cov(theta_egger_DP))),
  theta_lasso=colMeans(theta_lasso_DP),
  se_theta_lasso=sqrt(diag(cov(theta_lasso_DP))),
  theta_ivw=colMeans(theta_ivw_DP),
  se_theta_ivw=sqrt(diag(cov(theta_ivw_DP))),
  theta_median=colMeans(theta_median_DP),
  se_theta_median=sqrt(diag(cov(theta_median_DP)))),condF=colMeans(F_DP))
save(DBP_CAD_all,file="DBP_CAD_adjBMI_SBP.RData")


IV_mv=unique(c(IV_DBP,IV_BMI))
NIV_mv=length(IV_mv)
IV_beta_exp_mv=cbind(DBP$beta[DBP$SNP%in%IV_mv],BMI$beta[BMI$SNP%in%IV_mv] )
IV_se_exp_mv=cbind(DBP$se[DBP$SNP%in%IV_mv],BMI$se[BMI$SNP%in%IV_mv] )
IV_beta_out_mv=CAD$beta[CAD$SNP%in%IV_mv]
IV_se_out_mv=CAD$se[CAD$SNP%in%IV_mv]
null_mv=( BMI$p>0.1&DBP$p>0.1&CAD$p>0.1)&( BMI$SNP%in%ind_snp&DBP$SNP%in%ind_snp&CAD$SNP%in%ind_snp)
null_z_mv=cbind(DBP$z[null_mv],BMI$z[null_mv] ,CAD$z[null_mv])
correlation_mv=matrix(0,nrow=3,ncol=3)
for(i in 1:3){
  for(j in 1:3){
    correlation_mv[i,j]=cor(null_z_mv[,i],null_z_mv[,j])
  }
}
SIG_mv=list()
SIG_inv_mv=list()
SIG_exp_mv=list()
for(i in 1:NIV_mv){
  SIG_mv[[i]]=diag(c(IV_se_exp_mv[i,],IV_se_out_mv[i]),nrow=3,ncol=3)%*%correlation_mv%*%diag(c(IV_se_exp_mv[i,],IV_se_out_mv[i]),nrow=3,ncol=3)
  SIG_inv_mv[[i]]=solve(SIG_mv[[i]])
  SIG_exp_mv[[i]]=SIG_mv[[i]][1:2,1:2]
}
theta_cml_DP_mv=matrix(0,ncol=2,nrow=DP)
theta_egger_DP_mv=matrix(0,ncol=2,nrow=DP)
theta_lasso_DP_mv=matrix(0,ncol=2,nrow=DP)
theta_median_DP_mv=matrix(0,ncol=2,nrow=DP)
theta_ivw_DP_mv=matrix(0,ncol=2,nrow=DP)
F_DP_mv=matrix(0,nrow=DP,ncol=2)


for(i in 1:DP){
  beta_exp_DP_mv=matrix(0,nrow=NIV_mv,ncol=2)
  beta_out_DP_mv=numeric()
  for(j in 1:NIV_mv){
    beta_DP_mv=mvrnorm(n=1,mu=c(IV_beta_exp_mv[j,],IV_beta_out_mv[j]),Sigma = SIG_mv[[j]])
    beta_exp_DP_mv[j,]=beta_DP_mv[1:2]
    beta_out_DP_mv[j]=beta_DP_mv[3]
    if(beta_exp_DP_mv[j,1]<0){
      beta_exp_DP_mv[j,]=-beta_exp_DP_mv[j,]
      beta_out_DP_mv[j]=-beta_out_DP_mv[j]
    }
  }
  object1=mr_mvinput(bx=beta_exp_DP_mv,bxse = IV_se_exp_mv,by=beta_out_DP_mv,byse = IV_se_out_mv)
  object2=format_mvmr(BXGs =beta_exp_DP_mv,BYG=beta_out_DP_mv,seBYG = IV_se_out_mv,seBXGs = IV_se_exp_mv,RSID=c(1:NIV_mv) )
  
  
  cml_out=MVmr_cML(b_exp = beta_exp_DP_mv,b_out = as.matrix(beta_out_DP_mv,nrow=NIV_mv,ncol=1),se_bx=IV_se_exp_mv,Sig_inv_l=SIG_inv_mv,n=339224,K_vec=c(1:(NIV_mv-3)),min_theta_range=0.01,max_theta_range=0.01,thres=0.001,maxit=500  )
  ivw_out=mr_mvivw(object = object1)
  lasso_out=mvmr_lasso(bx=beta_exp_DP_mv,by=beta_out_DP_mv,seby = IV_se_out_mv)
  egger_out=mr_mvegger(object = object1)
  median_out=mr_mvmedian(object = object1,iterations = 100)
  theta_cml_DP_mv[i,]=cml_out$BIC_theta 
  theta_ivw_DP_mv[i,]=ivw_out$Estimate 
  theta_lasso_DP_mv[i,]=lasso_out$th_post 
  theta_egger_DP_mv[i,]=egger_out$Estimate  
  theta_median_DP_mv[i,]=median_out$Estimate 
  weak=strength_mvmr(r_input = object2,gencov = SIG_exp_mv)[1,]
  for(j in 1:2){
    F_DP_mv[i,j]=weak[1,j]
  }
  
}
DBP_CADadjBMI=list(data.frame(
  theta_cml=colMeans(theta_cml_DP_mv),
  se_theta_cml=sqrt(diag(cov(theta_cml_DP_mv))),
  theta_egger=colMeans(theta_egger_DP_mv),
  se_theta_egger=sqrt(diag(cov(theta_egger_DP_mv))),
  theta_lasso=colMeans(theta_lasso_DP_mv),
  se_theta_lasso=sqrt(diag(cov(theta_lasso_DP_mv))),
  theta_ivw=colMeans(theta_ivw_DP_mv),
  se_theta_ivw=sqrt(diag(cov(theta_ivw_DP_mv))),
  theta_median=colMeans(theta_median_DP_mv),
  se_theta_median=sqrt(diag(cov(theta_median_DP_mv)))),condF_mv=colMeans(F_DP_mv))
save(DBP_CADadjBMI,file="DBP_CAD_adjBMI.RData")


IV_uv=unique(c(IV_DBP))
NIV_uv=length(IV_uv)
IV_beta_exp_uv=as.matrix(DBP$beta[DBP$SNP%in%IV_uv],nrow=NIV_uv,ncol=1 )
IV_se_exp_uv=as.matrix(DBP$se[DBP$SNP%in%IV_uv] ,nrow=NIV_uv,ncol=1  )
IV_beta_out_uv=CAD$beta[CAD$SNP%in%IV_uv]
IV_se_out_uv=CAD$se[CAD$SNP%in%IV_uv]
null_uv=( DBP$p>0.1&CAD$p>0.1)&( DBP$SNP%in%ind_snp&CAD$SNP%in%ind_snp)
null_z_uv=cbind(DBP$z[null_uv]  ,CAD$z[null_uv])
correlation_uv=matrix(0,nrow=2,ncol=2)
for(i in 1:2){
  for(j in 1:2){
    correlation_uv[i,j]=cor(null_z_uv[,i],null_z_uv[,j])
  }
}
SIG_uv=list()
SIG_inv_uv=list()

for(i in 1:NIV_uv){
  SIG_uv[[i]]=diag(c(IV_se_exp_uv[i,],IV_se_out_uv[i]),nrow=2,ncol=2)%*%correlation_uv%*%diag(c(IV_se_exp_uv[i,],IV_se_out_uv[i]),nrow=2,ncol=2)
  SIG_inv_uv[[i]]=solve(SIG_uv[[i]])
  
}
theta_cml_DP_uv=matrix(0,ncol=1,nrow=DP)
theta_egger_DP_uv=matrix(0,ncol=1,nrow=DP)
theta_lasso_DP_uv=matrix(0,ncol=1,nrow=DP)
theta_median_DP_uv=matrix(0,ncol=1,nrow=DP)
theta_ivw_DP_uv=matrix(0,ncol=1,nrow=DP)

unloadNamespace("MVMR")
library(MendelianRandomization)
for(i in 1:DP){
  beta_exp_DP_uv=matrix(0,nrow=NIV_uv,ncol=1)
  beta_out_DP_uv=numeric()
  for(j in 1:NIV_uv){
    beta_DP_uv=mvrnorm(n=1,mu=c(IV_beta_exp_uv[j,],IV_beta_out_uv[j]),Sigma = SIG_uv[[j]])
    beta_exp_DP_uv[j,]=beta_DP_uv[1]
    beta_out_DP_uv[j]=beta_DP_uv[2]
    if(beta_exp_DP_uv[j,1]<0){
      beta_exp_DP_uv[j,]=-beta_exp_DP_uv[j,]
      beta_out_DP_uv[j]=-beta_out_DP_uv[j]
    }
  }
  object1=mr_input(bx=as.vector(beta_exp_DP_uv),bxse = as.vector(IV_se_exp_uv),by=beta_out_DP_uv,byse = IV_se_out_uv)
  object2=mr_mvinput(bx= beta_exp_DP_uv ,bxse =  IV_se_exp_uv ,by=beta_out_DP_uv,byse = IV_se_out_uv)
  
  cml_out=MVmr_cML(b_exp = beta_exp_DP_uv,b_out = as.matrix(beta_out_DP_uv,nrow=NIV_uv,ncol=1),se_bx=IV_se_exp_uv,Sig_inv_l=SIG_inv_uv,n=547261,K_vec=c(1:(NIV_uv-2)),min_theta_range=0,max_theta_range=0,thres=0.001,maxit=500  )
  ivw_out=mr_mvivw(object = object2)
  lasso_out=mvmr_lasso(bx=beta_exp_DP_uv,by=beta_out_DP_uv,seby = IV_se_out_uv)
  egger_out=mr_mvegger(object = object2)
  median_out=MendelianRandomization::mr_median(object = object1,iterations = 100)
  theta_cml_DP_uv[i,]=cml_out$BIC_theta 
  theta_ivw_DP_uv[i,]=ivw_out$Estimate 
  theta_lasso_DP_uv[i,]=lasso_out$th_post 
  theta_egger_DP_uv[i,]=egger_out$Estimate  
  theta_median_DP_uv[i,]=median_out$Estimate 
  
  
}
DBP_CAD_UVMR=list(data.frame(
  theta_cml=colMeans(theta_cml_DP_uv),
  se_theta_cml=sqrt(diag(cov(theta_cml_DP_uv))),
  theta_egger=colMeans(theta_egger_DP_uv),
  se_theta_egger=sqrt(diag(cov(theta_egger_DP_uv))),
  theta_lasso=colMeans(theta_lasso_DP_uv),
  se_theta_lasso=sqrt(diag(cov(theta_lasso_DP_uv))),
  theta_ivw=colMeans(theta_ivw_DP_uv),
  se_theta_ivw=sqrt(diag(cov(theta_ivw_DP_uv))),
  theta_median=colMeans(theta_median_DP_uv),
  se_theta_median=sqrt(diag(cov(theta_median_DP_uv)))))
save(DBP_CAD_UVMR,file="DBP_CAD_UVMR.RData")