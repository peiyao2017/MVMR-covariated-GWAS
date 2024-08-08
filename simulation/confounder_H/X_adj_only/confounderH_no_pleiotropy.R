setwd("/panfs/jay/groups/20/panwei/wan01299/paper2/simulation/confounder_H/collider_bias_X_adj_corrected/")
N_non_null=40
N_H_only=40
N_null=3000
DP=100
Nsnps=N_non_null+N_null+N_H_only

sample_size=20000
snp_names=numeric()
for(i in 1:Nsnps){
  snp_names[i]=paste("snp",i)
}
library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-12, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
input_file=c(
  "confounder_Hdim=1_no_pleiotropy.RData",
  "confounder_Hdim=2_no_pleiotropy.RData",
  "confounder_Hdim=4_no_pleiotropy.RData" 
)


repeat_parallel=20
repetition=25
output_file=matrix(0,nrow=length(input_file),ncol=repeat_parallel)
for(i in 1:length(input_file)){
  for(j in 1:repeat_parallel){
    output_file[i,j]=paste("result",j,input_file[i],sep="")
  }
}
for(i1 in 1:length(input_file)){
  load(input_file[i1])
  betaGU=as.matrix(effects_confounder_H$betaGU,nrow=Nsnps) 
  betaGH=as.matrix(effects_confounder_H$betaGH,nrow=Nsnps) 
  betaGX=as.matrix(effects_confounder_H$betaGX,nrow=Nsnps)
  betaHX=as.matrix(effects_confounder_H$betaHX,ncol=1) 
  betaGY=as.matrix(effects_confounder_H$betaGY,nrow=Nsnps) 
  betaHY=as.matrix(effects_confounder_H$betaHY,ncol=1) 
  betaUX=effects_confounder_H$betaUX 
  betaUY=effects_confounder_H$betaUY 
  d=nrow(betaHY)
  betaUH=matrix(effects_confounder_H$betaUH,ncol=d,nrow=1) 
  theta=effects_confounder_H$theta
  
  for(i2 in 1:repeat_parallel){
    total_estimate=foreach(i3=1:repetition,.combine = "rbind")%dopar%{
           
 library(MVMR)

library(MASS)
library(glmnet)
library(quantreg)
library(robustbase)
library(Matrix)
library(MendelianRandomization)
library(mvtnorm)


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
Gmatrix=matrix(0,nrow=sample_size,ncol=Nsnps)
for(i in 1:ncol(Gmatrix)){
  Gmatrix[,i]=rbinom(n=nrow(Gmatrix),size = 2,prob = 0.3)
  Gmatrix[,i]=scale(Gmatrix[,i],center = TRUE,scale = FALSE)
}
U=Gmatrix%*%betaGU+rnorm(n=sample_size,mean=0,sd=1)
H_known=Gmatrix%*%betaGH+U%*%betaUH
VarEH=numeric()
 
EH=rmvnorm(n=sample_size,mean=rep(0,times=d),sigma = diag(1,nrow = d,ncol = d))
H=H_known+EH
X_known=Gmatrix%*%betaGX+U%*%betaUX+H%*%betaHX
VarEX=(0.1/0.9)*var(X_known)
X=X_known+rnorm(n=sample_size,mean=0,sd=sqrt(1))
Y_known=Gmatrix%*%betaGY+H%*%betaHY+U%*%betaUY+theta*X
VarEY=(0.1/0.9)*var(Y_known)
Y=Y_known+rnorm(n=sample_size,mean=0,sd=sqrt(1))

gwas_H=list()
for(i in 1:d){
  temp=data.frame(snp=snp_names,beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  for(j in 1:Nsnps){
    dataH=data.frame(Htemp=H[,i],Gtemp=Gmatrix[,j])
    m=lm(Htemp~.-1,data=dataH)
    obj=summary(m)
    temp[j,2]=m$coefficients[length(m$coefficients)]
    temp[j,3]=obj[[4]][,2][length(obj[[4]][,2])]
    temp[j,4]=obj[[4]][,3][length(obj[[4]][,3])]
    temp[j,5]=obj[[4]][,4][length(obj[[4]][,4])]
    print(j)
  }
  gwas_H[[i]]=temp
}

gwas_X_adjusted=data.frame(snp=snp_names,beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
gwas_Y_unadjusted=data.frame(snp=snp_names,beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
for(i in 1:Nsnps){
  dataX_adjusted=data.frame(X=X,H,G=Gmatrix[,i])
  dataY_unadjusted=data.frame(Y=Y,G=Gmatrix[,i])
  m_X_adjusted=lm(X~.-1,data=dataX_adjusted)
  m_Y_unadjusted=lm(Y~.-1,data=dataY_unadjusted)
  obj_X_adjusted=summary(m_X_adjusted)
  
  
  obj_Y_unadjusted=summary(m_Y_unadjusted)
  gwas_X_adjusted[i,2]=m_X_adjusted$coefficients[length(m_X_adjusted$coefficients)]
  gwas_X_adjusted[i,3]=obj_X_adjusted[[4]][,2][length(obj_X_adjusted[[4]][,2])]
  gwas_X_adjusted[i,4]=obj_X_adjusted[[4]][,3][length(obj_X_adjusted[[4]][,3])]
  gwas_X_adjusted[i,5]=obj_X_adjusted[[4]][,4][length(obj_X_adjusted[[4]][,4])]
   
  
  gwas_Y_unadjusted[i,2]=m_Y_unadjusted$coefficients[length(m_Y_unadjusted$coefficients)]
  gwas_Y_unadjusted[i,3]=obj_Y_unadjusted[[4]][,2][length(obj_Y_unadjusted[[4]][,2])]
  gwas_Y_unadjusted[i,4]=obj_Y_unadjusted[[4]][,3][length(obj_Y_unadjusted[[4]][,3])]
  gwas_Y_unadjusted[i,5]=obj_Y_unadjusted[[4]][,4][length(obj_Y_unadjusted[[4]][,4])]
  print(i)
}

rm(Gmatrix)

 
correlation_X_adj=diag(1,nrow=d+2,ncol=d+2)
correlation_uv_X_adj=diag(1,nrow=2,ncol=2)
 
for(i in 1:d){
  null_z=gwas_X_adjusted$p>0.1&
    gwas_Y_unadjusted$p>0.1&
    gwas_H[[i]]$p>0.1
}
null_z_H=matrix(0,nrow = sum(null_z),ncol=d)
for(i in 1:d){
  null_z_H[,i]=gwas_H[[i]]$z[null_z]
}
 
null_z_X_adj=cbind(gwas_X_adjusted$z[null_z],null_z_H,gwas_Y_unadjusted$z[null_z])
null_z_uv_X_adj=cbind(gwas_X_adjusted$z[null_z],gwas_Y_unadjusted$z[null_z])
 
for(i in 1:nrow(correlation_X_adj)){
  for(j in 1:ncol(correlation_X_adj)){
     
    correlation_X_adj[i,j]=cor(null_z_X_adj[,i],null_z_X_adj[,j])
     
  }
}
for(i in 1:nrow(correlation_uv_X_adj)){
  for(j in 1:ncol(correlation_uv_X_adj)){
     
    correlation_uv_X_adj[i,j]=cor(null_z_uv_X_adj[,i],null_z_uv_X_adj[,j])
     
  }
}

beta_Y_no_adj=matrix(gwas_Y_unadjusted$beta[1:N_non_null],nrow = N_non_null,ncol=1)
beta_X_adj=matrix(gwas_X_adjusted$beta[1:N_non_null],nrow = N_non_null,ncol=1)
se_Y_no_adj=matrix(gwas_Y_unadjusted$se[1:N_non_null],nrow = N_non_null,ncol=1)
se_X_adj=matrix(gwas_X_adjusted$se[1:N_non_null],nrow = N_non_null,ncol=1)
beta_exp_adj=matrix(0,nrow=N_non_null,ncol = d+1)
beta_exp_adj[,1]=gwas_X_adjusted$beta[1:N_non_null]
se_exp_adj=matrix(0,nrow=N_non_null,ncol = d+1)
se_exp_adj[,1]=gwas_X_adjusted$se[1:N_non_null]

for(i in 1:d){
  beta_exp_adj[,i+1]=gwas_H[[i]]$beta[1:N_non_null]
  
  se_exp_adj[,i+1]=gwas_H[[i]]$se[1:N_non_null]
  
}
 
SIG_X_adjust=list()
SIG_uv_X_adjust=list()
SIG_inv_X_adjust=list()
SIG_inv_uv_X_adjust=list()
 
for(i in 1:N_non_null){
  a1=diag(c(se_exp_adj[i,],se_Y_no_adj[i]))%*%correlation_X_adj%*%diag(c(se_exp_adj[i,],se_Y_no_adj[i]))
  a2=diag(c(se_X_adj[i,],se_Y_no_adj[i]))%*%correlation_uv_X_adj%*%diag(c(se_X_adj[i,],se_Y_no_adj[i]))
  SIG_X_adjust[[i]]=a1
  SIG_uv_X_adjust[[i]]=a2
  SIG_inv_X_adjust[[i]]=solve(a1)
  SIG_inv_uv_X_adjust[[i]]=solve(a2)
  
}


if(d==1){
  library(SlopeHunter)
  library(indexevent)
  
  H_exp=matrix(0,ncol=d,nrow=N_H_only)
  H_exp_se=matrix(0,ncol=d,nrow=N_H_only)
  for(i in 1:d){
    H_exp[,i]=gwas_H[[i]]$beta[(N_non_null+1):(N_non_null+N_H_only)]
    H_exp_se[,i]=gwas_H[[i]]$se[(N_non_null+1):(N_non_null+N_H_only)]
  }
  X_out=gwas_X_adjusted$beta[(N_non_null+1):(N_non_null+N_H_only)]
  X_out_se=gwas_X_adjusted$se[(N_non_null+1):(N_non_null+N_H_only)]
   
  
  meanHX=list()
  covHX=list()
  
  cov_invHX=list()
   
   
   
  correlationHX=diag(x=1,nrow=d+1,ncol=d+1)
  
  null_zHX=cbind(null_z_X_adj[,2:(2+d-1)],null_z_X_adj[,1])
  
  
  for(i in 1:(d+1)){
    for(j in 1:(d+1)){
      correlationHX[i,j]=cor(null_zHX[,i],null_zHX[,j])
      
    }
  }
  for(i in 1:N_H_only){
    meanHX[[i]]=c(H_exp[i,],X_out[i])
    covHX[[i]]=diag(c(H_exp_se[i,],X_out_se[i]),nrow=d+1,ncol=d+1)%*%correlationHX%*%diag(c(H_exp_se[i,],X_out_se[i]),nrow=d+1,ncol=d+1)
    cov_invHX[[i]]=l=solve(covHX[[i]])
    
  }
  b_DP_X_cml=matrix(0,nrow=DP,ncol = d)
  b_DP_X_ivw=matrix(0,nrow=DP,ncol = d)
  b_DP_X_dho=matrix(0,nrow=DP,ncol = d)
  b_DP_X_sh=matrix(0,nrow=DP,ncol = d)
   
  for(i in 1:DP){
    H_exp_DP=matrix(0,nrow=N_H_only,ncol = d)
    X_out_DP=numeric()
    Y_out_DP=numeric()
    for(j in 1:N_H_only){
      gwas_DP=rmvnorm(n=1,mean=meanXHY[[j]],sigma = covXHY[[j]])
      H_exp_DP[j,]=gwas_DP[,2:(ncol(gwas_DP)-1)]
      X_out_DP[j]=gwas_DP[,1]
      Y_out_DP[j]=gwas_DP[,ncol(gwas_DP)]
    }
    object_DP_X=mr_mvinput(bx=H_exp_DP,bxse=H_exp_se,by=X_out_DP,byse = X_out_se)
    object_DP_Y=mr_mvinput(bx=H_exp_DP,bxse=H_exp_se,by=Y_out_DP,byse = Y_out_se)
    hunterX=data.frame(snp=1:N_H_only,xbeta=H_exp_DP,xse=H_exp_se,ybeta=X_out_DP,yse = X_out_se)
    hunterY=data.frame(snp=1:N_H_only,xbeta=H_exp_DP,xse=H_exp_se,ybeta=Y_out_DP,yse = Y_out_se)
    b_DP_X_cml[i,]=MVmr_cML(b_exp =H_exp_DP,b_out = as.matrix(X_out_DP),se_bx = H_exp_se,Sig_inv_l = cov_invHX,n=sample_size,random_start = 1,min_theta_range = -0.2,max_theta_range = 0.2,maxit = 200,thres = 0.01 )$BIC_theta
    b_DP_X_ivw[i,]=mr_mvivw(object =object_DP_X )$Estimate
    b_DP_X_dho[i,]=indexevent(xbeta = as.vector(H_exp_DP),xse=as.vector(H_exp_se),ybeta = X_out_DP,yse=X_out_se,method = "Hedges-Olkin")$b
    b_DP_X_sh[i,]=hunt(dat=hunterX,Bootstrapping = FALSE,snp_col ="snp",xbeta_col = "xbeta",xse_col = "xse",ybeta_col = "ybeta",yse_col = "yse",xp_thresh=1 ,Plot = FALSE)$b
    print(i)
  }
  b_X_cml=colMeans(b_DP_X_cml) 
  b_X_ivw=colMeans(b_DP_X_ivw) 
  b_X_dho=colMeans(b_DP_X_dho)  
  b_X_sh=colMeans(b_DP_X_sh) 
   
  covb_X_cml=cov(b_DP_X_cml) 
  covb_X_ivw=cov(b_DP_X_ivw) 
  covb_X_dho=cov(b_DP_X_dho)  
  covb_X_sh=cov(b_DP_X_sh) 
   
  namesX_b_cml=numeric()
  namesX_b_ivw=numeric()
  namesX_b_dho=numeric()
  namesX_b_sh=numeric()
  namesX_b_cml_se=numeric()
  namesX_b_ivw_se=numeric()
  namesX_b_dho_se=numeric()
  namesX_b_sh_se=numeric()
   
  for(i in 1:d){
    namesX_b_cml[i]=paste("b_X_cml", i ,sep="")
    namesX_b_ivw[i]=paste("b_X_ivw", i ,sep="")
    namesX_b_dho[i]=paste("b_X_dho", i ,sep="")
    namesX_b_sh[i]=paste("b_X_sh", i ,sep="")
    namesX_b_cml_se[i]=paste("seb_X_cml", i ,sep="")
    namesX_b_ivw_se[i]=paste("seb_X_ivw", i ,sep="")
    namesX_b_dho_se[i]=paste("seb_X_dho", i ,sep="")
    namesX_b_sh_se[i]=paste("seb_X_sh", i ,sep="")
 
  }
  namesb=c(namesX_b_cml,namesX_b_ivw,namesX_b_dho,namesX_b_sh,
           namesX_b_cml_se,namesX_b_ivw_se,namesX_b_dho_se,namesX_b_sh_se)
  btotal=matrix(c(b_X_cml=b_X_cml,b_X_ivw=b_X_ivw,b_X_dho=b_X_dho,b_X_sh=b_X_sh,
                  se_b_X_cml=diag(covb_X_cml),se_b_X_ivw=diag(covb_X_ivw),se_b_X_dho=diag(covb_X_dho),se_b_X_sh=diag(covb_X_sh)),nrow=1,ncol=8*d)
  colnames(btotal)=namesb           
  
 
  gwas_X_corrected_cml=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  gwas_X_corrected_ivw=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  gwas_X_corrected_dho=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  gwas_X_corrected_sh=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  H_se=matrix(0,ncol=d,nrow=Nsnps)
  H_beta=matrix(0,ncol=d,nrow=Nsnps)
  for(i in 1:d){
    H_se[,i]=gwas_H[[i]]$se
    H_beta[,i]=gwas_H[[i]]$beta
  }
  for(i in 1:Nsnps){
    
    gwas_X_corrected_cml$beta[i]=gwas_X_adjusted$beta[i]-b_X_cml%*%H_beta[i,]
    gwas_X_corrected_cml$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_cml)%*%H_se[i,]+t(b_X_cml)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_cml+t(H_beta[i,])%*%covb_X_cml%*%H_beta[i,])
    gwas_X_corrected_cml$z[i]=gwas_X_corrected_cml$beta[i]/gwas_X_corrected_cml$se[i]
    gwas_X_corrected_cml$p[i]=2*(1-pnorm(abs(gwas_X_corrected_cml$z[i])))
    gwas_X_corrected_ivw$beta[i]=gwas_X_adjusted$beta[i]-b_X_ivw%*%H_beta[i,]
    gwas_X_corrected_ivw$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_ivw)%*%H_se[i,]+t(b_X_ivw)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_ivw+t(H_beta[i,])%*%covb_X_ivw%*%H_beta[i,])
    gwas_X_corrected_ivw$z[i]=gwas_X_corrected_ivw$beta[i]/gwas_X_corrected_ivw$se[i]
    gwas_X_corrected_ivw$p[i]=2*(1-pnorm(abs(gwas_X_corrected_ivw$z[i])))
    gwas_X_corrected_dho$beta[i]=gwas_X_adjusted$beta[i]-b_X_dho%*%H_beta[i,]
    gwas_X_corrected_dho$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_dho)%*%H_se[i,]+t(b_X_dho)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_dho+t(H_beta[i,])%*%covb_X_dho%*%H_beta[i,])
    gwas_X_corrected_dho$z[i]=gwas_X_corrected_dho$beta[i]/gwas_X_corrected_dho$se[i]
    gwas_X_corrected_dho$p[i]=2*(1-pnorm(abs(gwas_X_corrected_dho$z[i])))
    gwas_X_corrected_sh$beta[i]=gwas_X_adjusted$beta[i]-b_X_sh%*%H_beta[i,]
    gwas_X_corrected_sh$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_sh)%*%H_se[i,]+t(b_X_sh)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_sh+t(H_beta[i,])%*%covb_X_sh%*%H_beta[i,])
    gwas_X_corrected_sh$z[i]=gwas_X_corrected_sh$beta[i]/gwas_X_corrected_sh$se[i]
    gwas_X_corrected_sh$p[i]=2*(1-pnorm(abs(gwas_X_corrected_sh$z[i])))
  }
  null_z_cml=cbind(gwas_X_corrected_cml$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_cml$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                   gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_cml$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  null_z_ivw=cbind(gwas_X_corrected_ivw$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_ivw$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                   gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_ivw$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  null_z_dho=cbind(gwas_X_corrected_dho$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_dho$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                   gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_dho$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  null_z_sh=cbind(gwas_X_corrected_sh$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_sh$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                  gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_sh$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  correlation_cml=diag(1,nrow=2,ncol=2)
  correlation_ivw=diag(1,nrow=2,ncol=2)
  correlation_dho=diag(1,nrow=2,ncol=2)
  correlation_sh=diag(1,nrow=2,ncol=2)
  for(i in 1:2){
    for(j in 1:2){
      correlation_cml[i,j]=cor(null_z_cml[,i],null_z_cml[,j])
      correlation_ivw[i,j]=cor(null_z_ivw[,i],null_z_ivw[,j])
      correlation_dho[i,j]=cor(null_z_dho[,i],null_z_dho[,j])
      correlation_sh[i,j]=cor(null_z_sh[,i],null_z_sh[,j])
    }
  }
  SIG_cml=list()
  SIG_ivw=list()
  SIG_dho=list()
  SIG_sh=list()
  
  SIGinv_cml=list()
  SIGinv_ivw=list()
  SIGinv_dho=list()
  SIGinv_sh=list()
  for(i in 1:N_non_null){
    a1=diag(c(gwas_X_corrected_cml$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_cml%*%diag(c(gwas_X_corrected_cml$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    a2=diag(c(gwas_X_corrected_ivw$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_ivw%*%diag(c(gwas_X_corrected_ivw$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    a3=diag(c(gwas_X_corrected_dho$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_dho%*%diag(c(gwas_X_corrected_dho$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    a4=diag(c(gwas_X_corrected_sh$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_sh%*%diag(c(gwas_X_corrected_sh$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    
    SIGinv_cml[[i]]=solve(a1)
    SIGinv_ivw[[i]]=solve(a2)
    SIGinv_dho[[i]]=solve(a3)
    SIGinv_sh[[i]]=solve(a4)
    
    SIG_cml[[i]]=a1 
    SIG_ivw[[i]]=a2 
    SIG_dho[[i]]=a3 
    SIG_sh[[i]]=a4 
  }
  
  theta_corrected_cml_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_cml_DP=matrix(0,nrow=DP,ncol = 1)
  X_cml=gwas_X_corrected_cml$beta[1:N_non_null]
  X_ivw=gwas_X_corrected_ivw$beta[1:N_non_null]
  Y_cml=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_ivw=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_cml=gwas_X_corrected_cml$se[1:N_non_null]
  seX_ivw=gwas_X_corrected_ivw$se[1:N_non_null]
  seY_cml=gwas_Y_unadjusted$se[1:N_non_null]
  seY_ivw=gwas_Y_unadjusted$se[1:N_non_null]
  
  X_dho=gwas_X_corrected_dho$beta[1:N_non_null]
  X_sh=gwas_X_corrected_sh$beta[1:N_non_null]
  Y_dho=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_sh=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_dho=gwas_X_corrected_dho$se[1:N_non_null]
  seX_sh=gwas_X_corrected_sh$se[1:N_non_null]
  seY_dho=gwas_Y_unadjusted$se[1:N_non_null]
  seY_sh=gwas_Y_unadjusted$se[1:N_non_null]
  obj_uv_corrected_cml=mr_mvinput(bx=matrix(X_cml),bxse =matrix(seX_cml),by=Y_cml,byse = seY_cml )
  obj_uv_corrected_ivw=mr_mvinput(bx=matrix(X_ivw),bxse =matrix(seX_ivw),by=Y_ivw,byse = seY_ivw )
  obj_uv_corrected_dho=mr_mvinput(bx=matrix(X_dho),bxse =matrix(seX_dho),by=Y_dho,byse = seY_dho )
  obj_uv_corrected_sh=mr_mvinput(bx=matrix(X_sh),bxse =matrix(seX_sh),by=Y_dho,byse = seY_sh )
  obj_uv_corrected_cml1=mr_input(bx=X_cml,bxse =seX_cml,by=Y_cml,byse = seY_cml )
  obj_uv_corrected_ivw1=mr_input(bx=X_ivw,bxse =seX_ivw,by=Y_ivw,byse = seY_ivw )
  obj_uv_corrected_dho1=mr_input(bx=X_dho,bxse =seX_dho,by=Y_dho,byse = seY_dho )
  obj_uv_corrected_sh1=mr_input(bx=X_sh,bxse =seX_sh,by=Y_sh,byse = seY_sh )
  est=MVmr_cML(b_exp =matrix(X_cml,ncol=1),b_out = matrix(Y_cml,ncol = 1),se_bx = matrix(seX_cml),
               Sig_inv_l = SIGinv_cml,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_cml,ncol=1),b_out = matrix(Y_cml,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_cml)
  theta_corrected_uv_cml_cml=theta_all[1]
  theta_corrected_uv_cml_cml_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_cml1 )
  theta_corrected_uv_ivw_cml=est$Estimate[1]
  theta_corrected_uv_ivw_cml_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_cml),seby = seY_cml,by=Y_cml)
  theta_corrected_uv_lasso_cml=est$th_post[1]
  theta_corrected_uv_lasso_cml_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_cml,orientate = 1)
  theta_corrected_uv_egger_cml=est$Estimate[1]
  theta_corrected_uv_egger_cml_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_cml1,iterations = 100)
  theta_corrected_uv_median_cml=est$Estimate[1]
  theta_corrected_uv_median_cml_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =matrix(X_ivw,ncol=1),b_out = matrix(Y_ivw,ncol = 1),se_bx = matrix(seX_ivw),
               Sig_inv_l = SIGinv_ivw,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_ivw,ncol=1),b_out = matrix(Y_ivw,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_ivw)
  theta_corrected_uv_cml_ivw=theta_all[1]
  theta_corrected_uv_cml_ivw_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_ivw1 )
  theta_corrected_uv_ivw_ivw=est$Estimate[1]
  theta_corrected_uv_ivw_ivw_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_ivw),seby = seY_ivw,by=Y_ivw)
  theta_corrected_uv_lasso_ivw=est$th_post[1]
  theta_corrected_uv_lasso_ivw_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_ivw,orientate = 1)
  theta_corrected_uv_egger_ivw=est$Estimate[1]
  theta_corrected_uv_egger_ivw_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_ivw1,iterations = 100)
  theta_corrected_uv_median_ivw=est$Estimate[1]
  theta_corrected_uv_median_ivw_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =matrix(X_dho,ncol=1),b_out = matrix(Y_dho,ncol = 1),se_bx = matrix(seX_dho),
               Sig_inv_l = SIGinv_dho,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_dho,ncol=1),b_out = matrix(Y_dho,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_dho)
  theta_corrected_uv_cml_dho=theta_all[1]
  theta_corrected_uv_cml_dho_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_dho1 )
  theta_corrected_uv_ivw_dho=est$Estimate[1]
  theta_corrected_uv_ivw_dho_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_dho),seby = seY_dho,by=Y_dho)
  theta_corrected_uv_lasso_dho=est$th_post[1]
  theta_corrected_uv_lasso_dho_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_dho,orientate = 1)
  theta_corrected_uv_egger_dho=est$Estimate[1]
  theta_corrected_uv_egger_dho_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_dho1,iterations = 100)
  theta_corrected_uv_median_dho=est$Estimate[1]
  theta_corrected_uv_median_dho_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =matrix(X_sh,ncol=1),b_out = matrix(Y_sh,ncol = 1),se_bx = matrix(seX_sh),
               Sig_inv_l = SIGinv_sh,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start =4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_sh,ncol=1),b_out = matrix(Y_sh,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_sh)
  theta_corrected_uv_cml_sh=theta_all[1]
  theta_corrected_uv_cml_sh_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_sh1 )
  theta_corrected_uv_ivw_sh=est$Estimate[1]
  theta_corrected_uv_ivw_sh_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_sh),seby = seY_sh,by=Y_sh)
  theta_corrected_uv_lasso_sh=est$th_post[1]
  theta_corrected_uv_lasso_sh_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_sh,orientate = 1)
  theta_corrected_uv_egger_sh=est$Estimate[1]
  theta_corrected_uv_egger_sh_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_sh1,iterations = 100)
  theta_corrected_uv_median_sh=est$Estimate[1]
  theta_corrected_uv_median_sh_se_raw=est$StdError[1]
  
  
  theta_corrected_cml_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_cml_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_cml=numeric()
    Y_DP_cml=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_cml[j],Y_cml[j]),Sigma = SIG_cml[[j]])
      X_DP_cml[j]=DP_IV[1]
      Y_DP_cml[j]=DP_IV[2]
    }
    obj_corrected_cml=mr_mvinput(bx=as.matrix(X_DP_cml),bxse = as.matrix(seX_cml),by=Y_DP_cml,byse = seY_cml)
    obj_corrected_cml1=mr_input(bx=X_DP_cml,bxse = seX_cml,by=Y_DP_cml,byse = seY_cml)
    theta_corrected_ivw_cml_DP[i,]=mr_ivw(object = obj_corrected_cml1)$Estimate[1]
    theta_corrected_median_cml_DP[i,]=mr_median(object = obj_corrected_cml1,iterations = 200)$Estimate[1]
    theta_corrected_egger_cml_DP[i,]=mr_mvegger(object = obj_corrected_cml,orientate = 1)$Estimate[1]
    theta_corrected_lasso_cml_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_cml),by=Y_DP_cml,seby = seY_cml)$th_post[1]
    theta_corrected_cml_cml_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_cml),b_out =as.matrix(Y_DP_cml),se_bx = as.matrix(seX_cml),Sig_inv_l = SIGinv_cml,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
   
  theta_corrected_uv_cml_cml_se=sd(theta_corrected_cml_cml_DP)
  theta_corrected_uv_ivw_cml_se=sd(theta_corrected_ivw_cml_DP)
  theta_corrected_uv_lasso_cml_se=sd(theta_corrected_lasso_cml_DP)
  theta_corrected_uv_egger_cml_se=sd(theta_corrected_egger_cml_DP)
  theta_corrected_uv_median_cml_se=sd(theta_corrected_median_cml_DP)
  
  theta_corrected_cml_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_ivw=numeric()
    Y_DP_ivw=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_ivw[j],Y_ivw[j]),Sigma = SIG_ivw[[j]])
      X_DP_ivw[j]=DP_IV[1]
      Y_DP_ivw[j]=DP_IV[2]
    }
    obj_corrected_ivw=mr_mvinput(bx=as.matrix(X_DP_ivw),bxse = as.matrix(seX_ivw),by=Y_DP_ivw,byse = seY_ivw)
    obj_corrected_ivw1=mr_input(bx=X_DP_ivw,bxse = seX_ivw,by=Y_DP_ivw,byse = seY_ivw)
    theta_corrected_ivw_ivw_DP[i,]=mr_ivw(object = obj_corrected_ivw1)$Estimate[1]
    theta_corrected_median_ivw_DP[i,]=mr_median(object = obj_corrected_ivw1,iterations = 200)$Estimate[1]
    theta_corrected_egger_ivw_DP[i,]=mr_mvegger(object = obj_corrected_ivw,orientate = 1)$Estimate[1]
    theta_corrected_lasso_ivw_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_ivw),by=Y_DP_ivw,seby = seY_ivw)$th_post[1]
    theta_corrected_cml_ivw_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_ivw),b_out =as.matrix(Y_DP_ivw),se_bx = as.matrix(seX_ivw),Sig_inv_l = SIGinv_ivw,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  
  
   
  theta_corrected_uv_cml_ivw_se=sd(theta_corrected_cml_ivw_DP)
  theta_corrected_uv_ivw_ivw_se=sd(theta_corrected_ivw_ivw_DP)
  theta_corrected_uv_lasso_ivw_se=sd(theta_corrected_lasso_ivw_DP)
  theta_corrected_uv_egger_ivw_se=sd(theta_corrected_egger_ivw_DP)
  theta_corrected_uv_median_ivw_se=sd(theta_corrected_median_ivw_DP)
  
  theta_corrected_cml_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_dho_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_dho=numeric()
    Y_DP_dho=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_dho[j],Y_dho[j]),Sigma = SIG_dho[[j]])
      X_DP_dho[j]=DP_IV[1]
      Y_DP_dho[j]=DP_IV[2]
    }
    obj_corrected_dho=mr_mvinput(bx=as.matrix(X_DP_dho),bxse = as.matrix(seX_dho),by=Y_DP_dho,byse = seY_dho)
    obj_corrected_dho1=mr_input(bx=X_DP_dho,bxse = seX_dho,by=Y_DP_dho,byse = seY_dho)
    theta_corrected_ivw_dho_DP[i,]=mr_ivw(object = obj_corrected_dho1)$Estimate[1]
    theta_corrected_median_dho_DP[i,]=mr_median(object = obj_corrected_dho1,iterations = 200)$Estimate[1]
    theta_corrected_egger_dho_DP[i,]=mr_mvegger(object = obj_corrected_dho,orientate = 1)$Estimate[1]
    theta_corrected_lasso_dho_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_dho),by=Y_DP_dho,seby = seY_dho)$th_post[1]
    theta_corrected_cml_dho_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_dho),b_out =as.matrix(Y_DP_dho),se_bx = as.matrix(seX_dho),Sig_inv_l = SIGinv_dho,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  
  
   
  theta_corrected_uv_cml_dho_se=sd(theta_corrected_cml_dho_DP)
  theta_corrected_uv_ivw_dho_se=sd(theta_corrected_ivw_dho_DP)
  theta_corrected_uv_lasso_dho_se=sd(theta_corrected_lasso_dho_DP)
  theta_corrected_uv_egger_dho_se=sd(theta_corrected_egger_dho_DP)
  theta_corrected_uv_median_dho_se=sd(theta_corrected_median_dho_DP)
  
  theta_corrected_cml_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_sh_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_sh=numeric()
    Y_DP_sh=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_sh[j],Y_sh[j]),Sigma = SIG_sh[[j]])
      X_DP_sh[j]=DP_IV[1]
      Y_DP_sh[j]=DP_IV[2]
    }
    obj_corrected_sh=mr_mvinput(bx=as.matrix(X_DP_sh),bxse = as.matrix(seX_sh),by=Y_DP_sh,byse = seY_sh)
    obj_corrected_sh1=mr_input(bx=X_DP_sh,bxse = seX_sh,by=Y_DP_sh,byse = seY_sh)
    theta_corrected_ivw_sh_DP[i,]=mr_ivw(object = obj_corrected_sh1)$Estimate[1]
    theta_corrected_median_sh_DP[i,]=mr_median(object = obj_corrected_sh1,iterations = 200)$Estimate[1]
    theta_corrected_egger_sh_DP[i,]=mr_mvegger(object = obj_corrected_sh,orientate = 1)$Estimate[1]
    theta_corrected_lasso_sh_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_sh),by=Y_DP_sh,seby = seY_sh)$th_post[1]
    theta_corrected_cml_sh_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_sh),b_out =as.matrix(Y_DP_sh),se_bx = as.matrix(seX_sh),Sig_inv_l = SIGinv_sh,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  
  
   
  theta_corrected_uv_cml_sh_se=sd(theta_corrected_cml_sh_DP)
  theta_corrected_uv_ivw_sh_se=sd(theta_corrected_ivw_sh_DP)
  theta_corrected_uv_lasso_sh_se=sd(theta_corrected_lasso_sh_DP)
  theta_corrected_uv_egger_sh_se=sd(theta_corrected_egger_sh_DP)
  theta_corrected_uv_median_sh_se=sd(theta_corrected_median_sh_DP)
  
  num=1:Nsnps
  nullz_mv_cml=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_cml$p>0.1&gwas_Y_unadjusted$p>0.1
  nullz_mv_ivw=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_ivw$p>0.1&gwas_Y_unadjusted$p>0.1
  nullz_mv_dho=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_dho$p>0.1&gwas_Y_unadjusted$p>0.1
  nullz_mv_sh=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_sh$p>0.1&gwas_Y_unadjusted$p>0.1
  for(i in 1:d){
    nullz_mv_cml=nullz_mv_cml&gwas_H[[i]]$p>0.1
    nullz_mv_ivw=nullz_mv_ivw&gwas_H[[i]]$p>0.1
    nullz_mv_dho=nullz_mv_dho&gwas_H[[i]]$p>0.1
    nullz_mv_sh=nullz_mv_sh&gwas_H[[i]]$p>0.1
  }
  
  
  
  for(i in 1:d){
    if(i==1){
      null_z_mv_cml=cbind(gwas_X_corrected_cml$z[nullz_mv_cml],gwas_H[[i]]$z[nullz_mv_cml])
      null_z_mv_ivw=cbind(gwas_X_corrected_ivw$z[nullz_mv_ivw],gwas_H[[i]]$z[nullz_mv_ivw])
      null_z_mv_dho=cbind(gwas_X_corrected_dho$z[nullz_mv_dho],gwas_H[[i]]$z[nullz_mv_dho])
      null_z_mv_sh=cbind(gwas_X_corrected_sh$z[nullz_mv_sh],gwas_H[[i]]$z[nullz_mv_sh])
    }
    if(i>1){
      null_z_mv_cml=cbind(null_z_mv_cml,gwas_H[[i]]$z[nullz_mv_cml])
      null_z_mv_ivw=cbind(null_z_mv_ivw,gwas_H[[i]]$z[nullz_mv_ivw])
      null_z_mv_dho=cbind(null_z_mv_dho,gwas_H[[i]]$z[nullz_mv_dho])
      null_z_mv_sh=cbind(null_z_mv_sh,gwas_H[[i]]$z[nullz_mv_sh])
    }
  }
  
  null_z_mv_cml=cbind(null_z_mv_cml,gwas_Y_unadjusted$z[nullz_mv_cml])
  null_z_mv_ivw=cbind(null_z_mv_ivw,gwas_Y_unadjusted$z[nullz_mv_ivw])
  null_z_mv_dho=cbind(null_z_mv_dho,gwas_Y_unadjusted$z[nullz_mv_dho])
  null_z_mv_sh=cbind(null_z_mv_sh,gwas_Y_unadjusted$z[nullz_mv_sh])
  
  
  
  correlation_mv_cml=diag(1,nrow=d+2,ncol=d+2)
  correlation_mv_ivw=diag(1,nrow=d+2,ncol=d+2)
  correlation_mv_dho=diag(1,nrow=d+2,ncol=d+2)
  correlation_mv_sh=diag(1,nrow=d+2,ncol=d+2)
  for(i in 1:(d+2)){
    for(j in 1:(d+2)){
      correlation_mv_cml[i,j]=cor(null_z_mv_cml[,i],null_z_mv_cml[,j])
      correlation_mv_ivw[i,j]=cor(null_z_mv_ivw[,i],null_z_mv_ivw[,j])
      correlation_mv_dho[i,j]=cor(null_z_mv_dho[,i],null_z_mv_dho[,j])
      correlation_mv_sh[i,j]=cor(null_z_mv_sh[,i],null_z_mv_sh[,j])
    }
  }
  SIG_mv_cml=list()
  SIG_mv_ivw=list()
  SIG_mv_dho=list()
  SIG_mv_sh=list()
  
  SIGinv_mv_cml=list()
  SIGinv_mv_ivw=list()
  SIGinv_mv_dho=list()
  SIGinv_mv_sh=list()
  for(i in 1:N_non_null){
    a1=diag(c(gwas_X_corrected_cml$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_cml%*%diag(c(gwas_X_corrected_cml$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    a2=diag(c(gwas_X_corrected_ivw$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_ivw%*%diag(c(gwas_X_corrected_ivw$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    a3=diag(c(gwas_X_corrected_dho$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_dho%*%diag(c(gwas_X_corrected_dho$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    a4=diag(c(gwas_X_corrected_sh$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_sh%*%diag(c(gwas_X_corrected_sh$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    
    SIGinv_mv_cml[[i]]=solve(a1)
    SIGinv_mv_ivw[[i]]=solve(a2)
    SIGinv_mv_dho[[i]]=solve(a3)
    SIGinv_mv_sh[[i]]=solve(a4)
    
    SIG_mv_cml[[i]]=a1 
    SIG_mv_ivw[[i]]=a2 
    SIG_mv_dho[[i]]=a3 
    SIG_mv_sh[[i]]=a4 
  }
  
  
  X_mv_cml=cbind(gwas_X_corrected_cml$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  X_mv_ivw=cbind(gwas_X_corrected_ivw$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  Y_cml=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_ivw=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_mv_cml=cbind(gwas_X_corrected_cml$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seX_mv_ivw=cbind(gwas_X_corrected_ivw$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seY_cml=gwas_Y_unadjusted$se[1:N_non_null]
  seY_ivw=gwas_Y_unadjusted$se[1:N_non_null]
  
  X_mv_dho=cbind(gwas_X_corrected_dho$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  X_mv_sh=cbind(gwas_X_corrected_sh$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  Y_dho=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_sh=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_mv_dho=cbind(gwas_X_corrected_dho$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seX_mv_sh=cbind(gwas_X_corrected_sh$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seY_dho=gwas_Y_unadjusted$se[1:N_non_null]
  seY_sh=gwas_Y_unadjusted$se[1:N_non_null]
  
  obj_mv_corrected_cml=mr_mvinput(bx=X_mv_cml,bxse =seX_mv_cml,by=Y_cml,byse = seY_cml )
  obj_mv_corrected_ivw=mr_mvinput(bx=X_mv_ivw,bxse =seX_mv_ivw,by=Y_ivw,byse = seY_ivw )
  obj_mv_corrected_dho=mr_mvinput(bx=X_mv_dho,bxse =seX_mv_dho,by=Y_dho,byse = seY_dho )
  obj_mv_corrected_sh=mr_mvinput(bx=X_mv_sh,bxse =seX_mv_sh,by=Y_dho,byse = seY_sh )
  
  est=MVmr_cML(b_exp =X_mv_cml,b_out = matrix(Y_cml,ncol = 1),se_bx =seX_mv_cml,
               Sig_inv_l = SIGinv_mv_cml,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_cml,b_out = matrix(Y_cml,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_cml)
  theta_corrected_mv_cml_cml=theta_all[1]
  theta_corrected_mv_cml_cml_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_cml )
  theta_corrected_mv_ivw_cml=est$Estimate[1]
  theta_corrected_mv_ivw_cml_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_cml,seby = seY_cml,by=Y_cml)
  theta_corrected_mv_lasso_cml=est$th_post[1]
  theta_corrected_mv_lasso_cml_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_cml,orientate = 1)
  theta_corrected_mv_egger_cml=est$Estimate[1]
  theta_corrected_mv_egger_cml_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_cml,iterations = 100)
  theta_corrected_mv_median_cml=est$Estimate[1]
  theta_corrected_mv_median_cml_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =X_mv_ivw,b_out = matrix(Y_ivw,ncol = 1),se_bx =seX_mv_ivw,
               Sig_inv_l = SIGinv_mv_ivw,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_ivw,b_out = matrix(Y_ivw,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_ivw)
  theta_corrected_mv_cml_ivw=theta_all[1]
  theta_corrected_mv_cml_ivw_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_ivw )
  theta_corrected_mv_ivw_ivw=est$Estimate[1]
  theta_corrected_mv_ivw_ivw_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_ivw,seby = seY_ivw,by=Y_ivw)
  theta_corrected_mv_lasso_ivw=est$th_post[1]
  theta_corrected_mv_lasso_ivw_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_ivw,orientate = 1)
  theta_corrected_mv_egger_ivw=est$Estimate[1]
  theta_corrected_mv_egger_ivw_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_ivw,iterations = 100)
  theta_corrected_mv_median_ivw=est$Estimate[1]
  theta_corrected_mv_median_ivw_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =X_mv_dho,b_out = matrix(Y_dho,ncol = 1),se_bx =seX_mv_dho,
               Sig_inv_l = SIGinv_mv_dho,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_dho,b_out = matrix(Y_dho,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_dho)
  theta_corrected_mv_cml_dho=theta_all[1]
  theta_corrected_mv_cml_dho_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_dho )
  theta_corrected_mv_ivw_dho=est$Estimate[1]
  theta_corrected_mv_ivw_dho_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_dho,seby = seY_dho,by=Y_dho)
  theta_corrected_mv_lasso_dho=est$th_post[1]
  theta_corrected_mv_lasso_dho_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_dho,orientate = 1)
  theta_corrected_mv_egger_dho=est$Estimate[1]
  theta_corrected_mv_egger_dho_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_dho,iterations = 100)
  theta_corrected_mv_median_dho=est$Estimate[1]
  theta_corrected_mv_median_dho_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =X_mv_sh,b_out = matrix(Y_sh,ncol = 1),se_bx =seX_mv_sh,
               Sig_inv_l = SIGinv_mv_sh,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_sh,b_out = matrix(Y_sh,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_sh)
  theta_corrected_mv_cml_sh=theta_all[1]
  theta_corrected_mv_cml_sh_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_sh )
  theta_corrected_mv_ivw_sh=est$Estimate[1]
  theta_corrected_mv_ivw_sh_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_sh,seby = seY_sh,by=Y_sh)
  theta_corrected_mv_lasso_sh=est$th_post[1]
  theta_corrected_mv_lasso_sh_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_sh,orientate = 1)
  theta_corrected_mv_egger_sh=est$Estimate[1]
  theta_corrected_mv_egger_sh_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_sh,iterations = 100)
  theta_corrected_mv_median_sh=est$Estimate[1]
  theta_corrected_mv_median_sh_se_raw=est$StdError[1]
  
  
  theta_corrected_mv_cml_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_cml_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_cml=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_cml=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_cml[j,],Y_cml[j]),Sigma = SIG_mv_cml[[j]])
      X_DP_cml[j,]=DP_IV[1:(d+1)]
      Y_DP_cml[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_cml[j,1]<0){
        X_DP_cml[j,]=-X_DP_cml[j,]
        Y_DP_cml[j]=-Y_DP_cml[j]
      }
    }
    obj_corrected_cml=mr_mvinput(bx=X_DP_cml,bxse = seX_mv_cml,by=Y_DP_cml,byse = seY_cml)
    
    theta_corrected_mv_ivw_cml_DP[i,]=mr_mvivw(object = obj_corrected_cml)$Estimate[1]
    theta_corrected_mv_median_cml_DP[i,]=mr_mvmedian(object = obj_corrected_cml,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_cml_DP[i,]=mr_mvegger(object = obj_corrected_cml,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_cml_DP[i,]=mvmr_lasso(bx=X_DP_cml,by=Y_DP_cml,seby = seY_cml)$th_post[1]
    theta_corrected_mv_cml_cml_DP[i,]=MVmr_cML(b_exp=X_DP_cml,b_out =as.matrix(Y_DP_cml),se_bx =seX_mv_cml,Sig_inv_l = SIGinv_mv_cml,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
   
  theta_corrected_mv_cml_cml_se=sd(theta_corrected_mv_cml_cml_DP)
  theta_corrected_mv_ivw_cml_se=sd(theta_corrected_mv_ivw_cml_DP)
  theta_corrected_mv_lasso_cml_se=sd(theta_corrected_mv_lasso_cml_DP)
  theta_corrected_mv_egger_cml_se=sd(theta_corrected_mv_egger_cml_DP)
  theta_corrected_mv_median_cml_se=sd(theta_corrected_mv_median_cml_DP)
  
  theta_corrected_mv_cml_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_ivw=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_ivw=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_ivw[j,],Y_ivw[j]),Sigma = SIG_mv_ivw[[j]])
      X_DP_ivw[j,]=DP_IV[1:(d+1)]
      Y_DP_ivw[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_ivw[j,1]<0){
        X_DP_ivw[j,]=-X_DP_ivw[j,]
        Y_DP_ivw[j]=-Y_DP_ivw[j]
      }
    }
    obj_corrected_ivw=mr_mvinput(bx=X_DP_ivw,bxse = seX_mv_ivw,by=Y_DP_ivw,byse = seY_ivw)
    
    theta_corrected_mv_ivw_ivw_DP[i,]=mr_mvivw(object = obj_corrected_ivw)$Estimate[1]
    theta_corrected_mv_median_ivw_DP[i,]=mr_mvmedian(object = obj_corrected_ivw,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_ivw_DP[i,]=mr_mvegger(object = obj_corrected_ivw,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_ivw_DP[i,]=mvmr_lasso(bx=X_DP_ivw,by=Y_DP_ivw,seby = seY_ivw)$th_post[1]
    theta_corrected_mv_cml_ivw_DP[i,]=MVmr_cML(b_exp=X_DP_ivw,b_out =as.matrix(Y_DP_ivw),se_bx =seX_mv_ivw,Sig_inv_l = SIGinv_mv_ivw,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  theta_corrected_mv_cml_ivw_se=sd(theta_corrected_mv_cml_ivw_DP)
  theta_corrected_mv_ivw_ivw_se=sd(theta_corrected_mv_ivw_ivw_DP)
  theta_corrected_mv_lasso_ivw_se=sd(theta_corrected_mv_lasso_ivw_DP)
  theta_corrected_mv_egger_ivw_se=sd(theta_corrected_mv_egger_ivw_DP)
  theta_corrected_mv_median_ivw_se=sd(theta_corrected_mv_median_ivw_DP)
  
  theta_corrected_mv_cml_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_dho_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_dho_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_dho=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_dho=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_dho[j,],Y_dho[j]),Sigma = SIG_mv_dho[[j]])
      X_DP_dho[j,]=DP_IV[1:(d+1)]
      Y_DP_dho[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_dho[j,1]<0){
        X_DP_dho[j, ]=-X_DP_dho[j, ]
        Y_DP_dho[j]=-Y_DP_dho[j]
      }
    }
    obj_corrected_dho=mr_mvinput(bx=X_DP_dho,bxse = seX_mv_dho,by=Y_DP_dho,byse = seY_dho)
    
    theta_corrected_mv_ivw_dho_DP[i,]=mr_mvivw(object = obj_corrected_dho)$Estimate[1]
    theta_corrected_mv_median_dho_DP[i,]=mr_mvmedian(object = obj_corrected_dho,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_dho_DP[i,]=mr_mvegger(object = obj_corrected_dho,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_dho_DP[i,]=mvmr_lasso(bx=X_DP_dho,by=Y_DP_dho,seby = seY_ivw)$th_post[1]
    theta_corrected_mv_cml_dho_DP[i,]=MVmr_cML(b_exp=X_DP_dho,b_out =as.matrix(Y_DP_dho),se_bx =seX_mv_dho,Sig_inv_l = SIGinv_mv_dho,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
   
  theta_corrected_mv_cml_dho_se=sd(theta_corrected_mv_cml_dho_DP)
  theta_corrected_mv_ivw_dho_se=sd(theta_corrected_mv_ivw_dho_DP)
  theta_corrected_mv_lasso_dho_se=sd(theta_corrected_mv_lasso_dho_DP)
  theta_corrected_mv_egger_dho_se=sd(theta_corrected_mv_egger_dho_DP)
  theta_corrected_mv_median_dho_se=sd(theta_corrected_mv_median_dho_DP)
  
  theta_corrected_mv_cml_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_sh_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_sh_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_sh=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_sh=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_sh[j,],Y_sh[j]),Sigma = SIG_mv_sh[[j]])
      X_DP_sh[j,]=DP_IV[1:(d+1)]
      Y_DP_sh[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_sh[j,1]<0){
        X_DP_sh[j, ]=-X_DP_sh[j, ]
        Y_DP_sh[j]=-Y_DP_sh[j]
      }
    }
    obj_corrected_sh=mr_mvinput(bx=X_DP_sh,bxse = seX_mv_sh,by=Y_DP_sh,byse = seY_sh)
    
    theta_corrected_mv_ivw_sh_DP[i,]=mr_mvivw(object = obj_corrected_sh)$Estimate[1]
    theta_corrected_mv_median_sh_DP[i,]=mr_mvmedian(object = obj_corrected_sh,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_sh_DP[i,]=mr_mvegger(object = obj_corrected_sh,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_sh_DP[i,]=mvmr_lasso(bx=X_DP_sh,by=Y_DP_sh,seby = seY_ivw)$th_post[1]
    theta_corrected_mv_cml_sh_DP[i,]=MVmr_cML(b_exp=X_DP_sh,b_out =as.matrix(Y_DP_sh),se_bx =seX_mv_sh,Sig_inv_l = SIGinv_mv_sh,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start = 1,min_theta_range =-0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  theta_corrected_mv_cml_sh_se=sd(theta_corrected_mv_cml_sh_DP)
  theta_corrected_mv_ivw_sh_se=sd(theta_corrected_mv_ivw_sh_DP)
  theta_corrected_mv_lasso_sh_se=sd(theta_corrected_mv_lasso_sh_DP)
  theta_corrected_mv_egger_sh_se=sd(theta_corrected_mv_egger_sh_DP)
  theta_corrected_mv_median_sh_se=sd(theta_corrected_mv_median_sh_DP)
  
}


if(d>1){
  library(SlopeHunter)
  library(indexevent)
  
  H_exp=matrix(0,ncol=d,nrow=N_H_only)
  H_exp_se=matrix(0,ncol=d,nrow=N_H_only)
  for(i in 1:d){
    H_exp[,i]=gwas_H[[i]]$beta[(N_non_null+1):(N_non_null+N_H_only)]
    H_exp_se[,i]=gwas_H[[i]]$se[(N_non_null+1):(N_non_null+N_H_only)]
  }
  X_out=gwas_X_adjusted$beta[(N_non_null+1):(N_non_null+N_H_only)]
  X_out_se=gwas_X_adjusted$se[(N_non_null+1):(N_non_null+N_H_only)]
  
  
  meanHX=list()
  covHX=list()
  
  cov_invHX=list()
  
  
  
  correlationHX=diag(x=1,nrow=d+1,ncol=d+1)
  
  null_zHX=cbind(null_z_X_adj[,2:(2+d-1)],null_z_X_adj[,1])
  
  
  for(i in 1:(d+1)){
    for(j in 1:(d+1)){
      correlationHX[i,j]=cor(null_zHX[,i],null_zHX[,j])
      
    }
  }
  for(i in 1:N_H_only){
    meanHX[[i]]=c(H_exp[i,],X_out[i])
    covHX[[i]]=diag(c(H_exp_se[i,],X_out_se[i]),nrow=d+1,ncol=d+1)%*%correlationHX%*%diag(c(H_exp_se[i,],X_out_se[i]),nrow=d+1,ncol=d+1)
    cov_invHX[[i]]=l=solve(covHX[[i]])
    
  }
  b_DP_X_cml=matrix(0,nrow=DP,ncol = d)
  b_DP_X_ivw=matrix(0,nrow=DP,ncol = d)
  
  for(i in 1:DP){
    H_exp_DP=matrix(0,nrow=N_H_only,ncol = d)
    X_out_DP=numeric()
    
    for(j in 1:N_H_only){
      gwas_DP=rmvnorm(n=1,mean=meanHX[[j]],sigma = covHX[[j]])
      H_exp_DP[j,]=gwas_DP[,1:(ncol(gwas_DP)-1)]
      X_out_DP[j]=gwas_DP[,ncol(gwas_DP)]
      
    }
    object_DP_X=mr_mvinput(bx=H_exp_DP,bxse=H_exp_se,by=X_out_DP,byse = X_out_se)
    
    hunterX=data.frame(snp=1:N_H_only,xbeta=H_exp_DP,xse=H_exp_se,ybeta=X_out_DP,yse = X_out_se)
    
    b_DP_X_cml[i,]=MVmr_cML(b_exp =H_exp_DP,b_out = as.matrix(X_out_DP),se_bx = H_exp_se,Sig_inv_l = cov_invHX,n=sample_size,random_start =1 ,min_theta_range = -0.4,max_theta_range = -0.4,maxit = 200,thres = 0.03 )$BIC_theta
    b_DP_X_ivw[i,]=mr_mvivw(object =object_DP_X )$Estimate
    print(i)
  }
  b_X_cml=colMeans(b_DP_X_cml) 
  b_X_ivw=colMeans(b_DP_X_ivw) 
  
    
  
  covb_X_cml=cov(b_DP_X_cml) 
  covb_X_ivw=cov(b_DP_X_ivw) 
  
  namesX_b_cml=numeric()
  namesX_b_ivw=numeric()
  namesX_b_cml_se=numeric()
  namesX_b_ivw_se=numeric()
   
  for(i in 1:d){
    namesX_b_cml[i]=paste("b_X_cml",i,sep="")
    namesX_b_ivw[i]=paste("b_X_ivw",i,sep="")
    namesX_b_cml_se[i]=paste("seb_X_cml",i,sep="")
    namesX_b_ivw_se[i]=paste("seb_X_ivw",i,sep="")
     
  }
  namesb=c(namesX_b_cml,namesX_b_ivw, 
           namesX_b_cml_se,namesX_b_ivw_se )
  btotal=matrix(c(b_X_cml=b_X_cml,b_X_ivw=b_X_ivw, 
                  se_b_X_cml=diag(covb_X_cml),se_b_X_ivw=diag(covb_X_ivw)),nrow=1,ncol=4*d )
  colnames(btotal)=namesb   
  
    
  
  gwas_X_corrected_cml=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  gwas_X_corrected_ivw=data.frame(beta=rep(0,times=Nsnps),se=rep(0,times=Nsnps),z=rep(0,times=Nsnps),p=rep(0,times=Nsnps))
  
  H_se=matrix(0,ncol=d,nrow=Nsnps)
  H_beta=matrix(0,ncol=d,nrow=Nsnps)
  for(i in 1:d){
    H_se[,i]=gwas_H[[i]]$se
    H_beta[,i]=gwas_H[[i]]$beta
  }
  for(i in 1:Nsnps){
    
    gwas_X_corrected_cml$beta[i]=gwas_X_adjusted$beta[i]-b_X_cml%*%H_beta[i,]
    gwas_X_corrected_cml$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_cml)%*%H_se[i,]+t(b_X_cml)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_cml+t(H_beta[i,])%*%covb_X_cml%*%H_beta[i,])
    gwas_X_corrected_cml$z[i]=gwas_X_corrected_cml$beta[i]/gwas_X_corrected_cml$se[i]
    gwas_X_corrected_cml$p[i]=2*(1-pnorm(abs(gwas_X_corrected_cml$z[i])))
    gwas_X_corrected_ivw$beta[i]=gwas_X_adjusted$beta[i]-b_X_ivw%*%H_beta[i,]
    gwas_X_corrected_ivw$se[i]=sqrt(gwas_X_adjusted$se[i]^2+t(H_se[i,])%*%(diag(1,nrow=d,ncol=d)*covb_X_ivw)%*%H_se[i,]+t(b_X_ivw)%*%diag(H_se[i,],nrow=d,ncol=d)%*%b_X_ivw+t(H_beta[i,])%*%covb_X_ivw%*%H_beta[i,])
    gwas_X_corrected_ivw$z[i]=gwas_X_corrected_ivw$beta[i]/gwas_X_corrected_ivw$se[i]
    gwas_X_corrected_ivw$p[i]=2*(1-pnorm(abs(gwas_X_corrected_ivw$z[i])))
    
  }
  null_z_cml=cbind(gwas_X_corrected_cml$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_cml$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                   gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_cml$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  null_z_ivw=cbind(gwas_X_corrected_ivw$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_ivw$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1],
                   gwas_Y_unadjusted$z[(N_non_null+N_H_only+1):Nsnps][gwas_X_corrected_ivw$p[(N_non_null+N_H_only+1):Nsnps]>0.1&gwas_Y_unadjusted$p[(N_non_null+N_H_only+1):Nsnps]>0.1])
  
  correlation_cml=diag(1,nrow=2,ncol=2)
  correlation_ivw=diag(1,nrow=2,ncol=2)
  
  for(i in 1:2){
    for(j in 1:2){
      correlation_cml[i,j]=cor(null_z_cml[,i],null_z_cml[,j])
      correlation_ivw[i,j]=cor(null_z_ivw[,i],null_z_ivw[,j])
      
    }
  }
  SIG_cml=list()
  SIG_ivw=list()
  
  
  SIGinv_cml=list()
  SIGinv_ivw=list()
  
  for(i in 1:N_non_null){
    a1=diag(c(gwas_X_corrected_cml$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_cml%*%diag(c(gwas_X_corrected_cml$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    a2=diag(c(gwas_X_corrected_ivw$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)%*%correlation_ivw%*%diag(c(gwas_X_corrected_ivw$se[i],gwas_Y_unadjusted$se[i]),nrow=2,ncol=2)
    
    SIGinv_cml[[i]]=solve(a1)
    SIGinv_ivw[[i]]=solve(a2)
    
    SIG_cml[[i]]=a1 
    SIG_ivw[[i]]=a2 
  }
  
  
  X_cml=gwas_X_corrected_cml$beta[1:N_non_null]
  X_ivw=gwas_X_corrected_ivw$beta[1:N_non_null]
  Y_cml=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_ivw=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_cml=gwas_X_corrected_cml$se[1:N_non_null]
  seX_ivw=gwas_X_corrected_ivw$se[1:N_non_null]
  seY_cml=gwas_Y_unadjusted$se[1:N_non_null]
  seY_ivw=gwas_Y_unadjusted$se[1:N_non_null]
  
  
  obj_uv_corrected_cml=mr_mvinput(bx=matrix(X_cml),bxse =matrix(seX_cml),by=Y_cml,byse = seY_cml )
  obj_uv_corrected_ivw=mr_mvinput(bx=matrix(X_ivw),bxse =matrix(seX_ivw),by=Y_ivw,byse = seY_ivw )
  
  obj_uv_corrected_cml1=mr_input(bx=X_cml,bxse =seX_cml,by=Y_cml,byse = seY_cml )
  obj_uv_corrected_ivw1=mr_input(bx=X_ivw,bxse =seX_ivw,by=Y_ivw,byse = seY_ivw )
  
  est=MVmr_cML(b_exp =matrix(X_cml,ncol=1),b_out = matrix(Y_cml,ncol = 1),se_bx = matrix(seX_cml),
               Sig_inv_l = SIGinv_cml,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_cml,ncol=1),b_out = matrix(Y_cml,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_cml)
  theta_corrected_uv_cml_cml=theta_all[1]
  theta_corrected_uv_cml_cml_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_cml1 )
  theta_corrected_uv_ivw_cml=est$Estimate[1]
  theta_corrected_uv_ivw_cml_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_cml),seby = seY_cml,by=Y_cml)
  theta_corrected_uv_lasso_cml=est$th_post[1]
  theta_corrected_uv_lasso_cml_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_cml,orientate = 1)
  theta_corrected_uv_egger_cml=est$Estimate[1]
  theta_corrected_uv_egger_cml_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_cml1,iterations = 100)
  theta_corrected_uv_median_cml=est$Estimate[1]
  theta_corrected_uv_median_cml_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =matrix(X_ivw,ncol=1),b_out = matrix(Y_ivw,ncol = 1),se_bx = matrix(seX_ivw),
               Sig_inv_l = SIGinv_ivw,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = matrix(X_ivw,ncol=1),b_out = matrix(Y_ivw,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_ivw)
  theta_corrected_uv_cml_ivw=theta_all[1]
  theta_corrected_uv_cml_ivw_se_raw=se_est[1]
  est=mr_ivw(object =obj_uv_corrected_ivw1 )
  theta_corrected_uv_ivw_ivw=est$Estimate[1]
  theta_corrected_uv_ivw_ivw_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=matrix(X_ivw),seby = seY_ivw,by=Y_ivw)
  theta_corrected_uv_lasso_ivw=est$th_post[1]
  theta_corrected_uv_lasso_ivw_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_uv_corrected_ivw,orientate = 1)
  theta_corrected_uv_egger_ivw=est$Estimate[1]
  theta_corrected_uv_egger_ivw_se_raw=est$StdError.Est[1]
  est=mr_median(object = obj_uv_corrected_ivw1,iterations = 100)
  theta_corrected_uv_median_ivw=est$Estimate[1]
  theta_corrected_uv_median_ivw_se_raw=est$StdError[1]
  
  
  theta_corrected_cml_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_cml_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_cml=numeric()
    Y_DP_cml=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_cml[j],Y_cml[j]),Sigma = SIG_cml[[j]])
      X_DP_cml[j]=DP_IV[1]
      Y_DP_cml[j]=DP_IV[2]
    }
    obj_corrected_cml=mr_mvinput(bx=as.matrix(X_DP_cml),bxse = as.matrix(seX_cml),by=Y_DP_cml,byse = seY_cml)
    obj_corrected_cml1=mr_input(bx=X_DP_cml,bxse = seX_cml,by=Y_DP_cml,byse = seY_cml)
    theta_corrected_ivw_cml_DP[i,]=mr_ivw(object = obj_corrected_cml1)$Estimate[1]
    theta_corrected_median_cml_DP[i,]=mr_median(object = obj_corrected_cml1,iterations = 200)$Estimate[1]
    theta_corrected_egger_cml_DP[i,]=mr_mvegger(object = obj_corrected_cml,orientate = 1)$Estimate[1]
    theta_corrected_lasso_cml_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_cml),by=Y_DP_cml,seby = seY_cml)$th_post[1]
    theta_corrected_cml_cml_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_cml),b_out =as.matrix(Y_DP_cml),se_bx = as.matrix(seX_cml),Sig_inv_l = SIGinv_cml,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  theta_corrected_uv_cml_cml_se=sd(theta_corrected_cml_cml_DP)
  theta_corrected_uv_ivw_cml_se=sd(theta_corrected_ivw_cml_DP)
  theta_corrected_uv_lasso_cml_se=sd(theta_corrected_lasso_cml_DP)
  theta_corrected_uv_egger_cml_se=sd(theta_corrected_egger_cml_DP)
  theta_corrected_uv_median_cml_se=sd(theta_corrected_median_cml_DP)
  
  theta_corrected_cml_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_ivw_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_egger_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_lasso_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_median_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_ivw=numeric()
    Y_DP_ivw=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_ivw[j],Y_ivw[j]),Sigma = SIG_ivw[[j]])
      X_DP_ivw[j]=DP_IV[1]
      Y_DP_ivw[j]=DP_IV[2]
    }
    obj_corrected_ivw=mr_mvinput(bx=as.matrix(X_DP_ivw),bxse = as.matrix(seX_ivw),by=Y_DP_ivw,byse = seY_ivw)
    obj_corrected_ivw1=mr_input(bx=X_DP_ivw,bxse = seX_ivw,by=Y_DP_ivw,byse = seY_ivw)
    theta_corrected_ivw_ivw_DP[i,]=mr_ivw(object = obj_corrected_ivw1)$Estimate[1]
    theta_corrected_median_ivw_DP[i,]=mr_median(object = obj_corrected_ivw1,iterations = 200)$Estimate[1]
    theta_corrected_egger_ivw_DP[i,]=mr_mvegger(object = obj_corrected_ivw,orientate = 1)$Estimate[1]
    theta_corrected_lasso_ivw_DP[i,]=mvmr_lasso(bx=as.matrix(X_DP_ivw),by=Y_DP_ivw,seby = seY_ivw)$th_post[1]
    theta_corrected_cml_ivw_DP[i,]=MVmr_cML(b_exp=as.matrix(X_DP_ivw),b_out =as.matrix(Y_DP_ivw),se_bx = as.matrix(seX_ivw),Sig_inv_l = SIGinv_ivw,n=sample_size,K_vec = c(0:(N_non_null-2)),random_start =1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  
  
  
  theta_corrected_uv_cml_ivw_se=sd(theta_corrected_cml_ivw_DP)
  theta_corrected_uv_ivw_ivw_se=sd(theta_corrected_ivw_ivw_DP)
  theta_corrected_uv_lasso_ivw_se=sd(theta_corrected_lasso_ivw_DP)
  theta_corrected_uv_egger_ivw_se=sd(theta_corrected_egger_ivw_DP)
  theta_corrected_uv_median_ivw_se=sd(theta_corrected_median_ivw_DP)
  
  
  
  num=1:Nsnps
  nullz_mv_cml=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_cml$p>0.1&gwas_Y_unadjusted$p>0.1
  nullz_mv_ivw=(num>=N_non_null+N_H_only+1)&gwas_X_corrected_ivw$p>0.1&gwas_Y_unadjusted$p>0.1
  
  for(i in 1:d){
    nullz_mv_cml=nullz_mv_cml&gwas_H[[i]]$p>0.1
    nullz_mv_ivw=nullz_mv_ivw&gwas_H[[i]]$p>0.1
    
  }
  
  
  
  for(i in 1:d){
    if(i==1){
      null_z_mv_cml=cbind(gwas_X_corrected_cml$z[nullz_mv_cml],gwas_H[[i]]$z[nullz_mv_cml])
      null_z_mv_ivw=cbind(gwas_X_corrected_ivw$z[nullz_mv_ivw],gwas_H[[i]]$z[nullz_mv_ivw])
      
    }
    if(i>1){
      null_z_mv_cml=cbind(null_z_mv_cml,gwas_H[[i]]$z[nullz_mv_cml])
      null_z_mv_ivw=cbind(null_z_mv_ivw,gwas_H[[i]]$z[nullz_mv_ivw])
      
    }
  }
  
  null_z_mv_cml=cbind(null_z_mv_cml,gwas_Y_unadjusted$z[nullz_mv_cml])
  null_z_mv_ivw=cbind(null_z_mv_ivw,gwas_Y_unadjusted$z[nullz_mv_ivw])
  
  
  
  correlation_mv_cml=diag(1,nrow=d+2,ncol=d+2)
  correlation_mv_ivw=diag(1,nrow=d+2,ncol=d+2)
  
  for(i in 1:(d+2)){
    for(j in 1:(d+2)){
      correlation_mv_cml[i,j]=cor(null_z_mv_cml[,i],null_z_mv_cml[,j])
      correlation_mv_ivw[i,j]=cor(null_z_mv_ivw[,i],null_z_mv_ivw[,j])
      
    }
  }
  SIG_mv_cml=list()
  SIG_mv_ivw=list()
  
  
  SIGinv_mv_cml=list()
  SIGinv_mv_ivw=list()
  
  for(i in 1:N_non_null){
    a1=diag(c(gwas_X_corrected_cml$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_cml%*%diag(c(gwas_X_corrected_cml$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    a2=diag(c(gwas_X_corrected_ivw$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)%*%correlation_mv_ivw%*%diag(c(gwas_X_corrected_ivw$se[i],se_exp_adj[i,][2:(1+d)],gwas_Y_unadjusted$se[i]),nrow=d+2,ncol=d+2)
    
    SIGinv_mv_cml[[i]]=solve(a1)
    SIGinv_mv_ivw[[i]]=solve(a2)
    
    
    SIG_mv_cml[[i]]=a1 
    SIG_mv_ivw[[i]]=a2 
    
  }
  
  
  X_mv_cml=cbind(gwas_X_corrected_cml$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  X_mv_ivw=cbind(gwas_X_corrected_ivw$beta[1:N_non_null],beta_exp_adj[,2:(d+1)])
  Y_cml=gwas_Y_unadjusted$beta[1:N_non_null]
  Y_ivw=gwas_Y_unadjusted$beta[1:N_non_null]
  seX_mv_cml=cbind(gwas_X_corrected_cml$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seX_mv_ivw=cbind(gwas_X_corrected_ivw$se[1:N_non_null],se_exp_adj[,2:(d+1)])
  seY_cml=gwas_Y_unadjusted$se[1:N_non_null]
  seY_ivw=gwas_Y_unadjusted$se[1:N_non_null]
  
  obj_mv_corrected_cml=mr_mvinput(bx=X_mv_cml,bxse =seX_mv_cml,by=Y_cml,byse = seY_cml )
  obj_mv_corrected_ivw=mr_mvinput(bx=X_mv_ivw,bxse =seX_mv_ivw,by=Y_ivw,byse = seY_ivw )
  
  
  est=MVmr_cML(b_exp =X_mv_cml,b_out = matrix(Y_cml,ncol = 1),se_bx =seX_mv_cml,
               Sig_inv_l = SIGinv_mv_cml,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_cml,b_out = matrix(Y_cml,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_cml)
  theta_corrected_mv_cml_cml=theta_all[1]
  theta_corrected_mv_cml_cml_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_cml )
  theta_corrected_mv_ivw_cml=est$Estimate[1]
  theta_corrected_mv_ivw_cml_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_cml,seby = seY_cml,by=Y_cml)
  theta_corrected_mv_lasso_cml=est$th_post[1]
  theta_corrected_mv_lasso_cml_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_cml,orientate = 1)
  theta_corrected_mv_egger_cml=est$Estimate[1]
  theta_corrected_mv_egger_cml_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_cml,iterations = 100)
  theta_corrected_mv_median_cml=est$Estimate[1]
  theta_corrected_mv_median_cml_se_raw=est$StdError[1]
  est=MVmr_cML(b_exp =X_mv_ivw,b_out = matrix(Y_ivw,ncol = 1),se_bx =seX_mv_ivw,
               Sig_inv_l = SIGinv_mv_ivw,K_vec = c(1:(N_non_null-1)),n=sample_size,min_theta_range = -0.5,
               max_theta_range = 0.5,thres = 0.03,random_start = 4)
  valid=setdiff(c(1:N_non_null),est$BIC_invalid)
  theta_all=est$BIC_theta
  se_est=MVcML_SdTheta(b_exp = X_mv_ivw,b_out = matrix(Y_ivw,ncol = 1),theta = theta_all,zero_ind = valid,r_vec =NULL,Sig_inv_l = SIGinv_mv_ivw)
  theta_corrected_mv_cml_ivw=theta_all[1]
  theta_corrected_mv_cml_ivw_se_raw=se_est[1]
  est=mr_mvivw(object =obj_mv_corrected_ivw )
  theta_corrected_mv_ivw_ivw=est$Estimate[1]
  theta_corrected_mv_ivw_ivw_se_raw=est$StdError[1]
  est=mvmr_lasso(bx=X_mv_ivw,seby = seY_ivw,by=Y_ivw)
  theta_corrected_mv_lasso_ivw=est$th_post[1]
  theta_corrected_mv_lasso_ivw_se_raw=est$se_post[1]
  est=mr_mvegger(object = obj_mv_corrected_ivw,orientate = 1)
  theta_corrected_mv_egger_ivw=est$Estimate[1]
  theta_corrected_mv_egger_ivw_se_raw=est$StdError.Est[1]
  est=mr_mvmedian(object = obj_mv_corrected_ivw,iterations = 100)
  theta_corrected_mv_median_ivw=est$Estimate[1]
  theta_corrected_mv_median_ivw_se_raw=est$StdError[1]
  
  
  
  theta_corrected_mv_cml_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_cml_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_cml_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_cml=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_cml=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_cml[j,],Y_cml[j]),Sigma = SIG_mv_cml[[j]])
      X_DP_cml[j,]=DP_IV[1:(d+1)]
      Y_DP_cml[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_cml[j,1]<0){
        X_DP_cml[j,]=-X_DP_cml[j,]
        Y_DP_cml[j]=-Y_DP_cml[j]
      }
    }
    obj_corrected_cml=mr_mvinput(bx=X_DP_cml,bxse = seX_mv_cml,by=Y_DP_cml,byse = seY_cml)
    
    theta_corrected_mv_ivw_cml_DP[i,]=mr_mvivw(object = obj_corrected_cml)$Estimate[1]
    theta_corrected_mv_median_cml_DP[i,]=mr_mvmedian(object = obj_corrected_cml,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_cml_DP[i,]=mr_mvegger(object = obj_corrected_cml,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_cml_DP[i,]=mvmr_lasso(bx=X_DP_cml,by=Y_DP_cml,seby = seY_cml)$th_post[1]
    theta_corrected_mv_cml_cml_DP[i,]=MVmr_cML(b_exp=X_DP_cml,b_out =as.matrix(Y_DP_cml),se_bx =seX_mv_cml,Sig_inv_l = SIGinv_mv_cml,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  theta_corrected_mv_cml_cml_se=sd(theta_corrected_mv_cml_cml_DP)
  theta_corrected_mv_ivw_cml_se=sd(theta_corrected_mv_ivw_cml_DP)
  theta_corrected_mv_lasso_cml_se=sd(theta_corrected_mv_lasso_cml_DP)
  theta_corrected_mv_egger_cml_se=sd(theta_corrected_mv_egger_cml_DP)
  theta_corrected_mv_median_cml_se=sd(theta_corrected_mv_median_cml_DP)
  
  theta_corrected_mv_cml_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_ivw_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_egger_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_lasso_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  theta_corrected_mv_median_ivw_DP=matrix(0,nrow=DP,ncol = 1)
  for(i in 1:DP){
    X_DP_ivw=matrix(0,n=N_non_null,ncol = d+1)
    Y_DP_ivw=numeric()
    for(j in 1:N_non_null){
      DP_IV=mvrnorm(n=1,mu=c(X_mv_ivw[j,],Y_ivw[j]),Sigma = SIG_mv_ivw[[j]])
      X_DP_ivw[j,]=DP_IV[1:(d+1)]
      Y_DP_ivw[j]=DP_IV[length(DP_IV)]
    }
    for(j in 1:N_non_null){
      if(X_DP_ivw[j,1]<0){
        X_DP_ivw[j,]=-X_DP_ivw[j,]
        Y_DP_ivw[j]=-Y_DP_ivw[j]
      }
    }
    obj_corrected_ivw=mr_mvinput(bx=X_DP_ivw,bxse = seX_mv_ivw,by=Y_DP_ivw,byse = seY_ivw)
    
    theta_corrected_mv_ivw_ivw_DP[i,]=mr_mvivw(object = obj_corrected_ivw)$Estimate[1]
    theta_corrected_mv_median_ivw_DP[i,]=mr_mvmedian(object = obj_corrected_ivw,iterations = 200)$Estimate[1]
    theta_corrected_mv_egger_ivw_DP[i,]=mr_mvegger(object = obj_corrected_ivw,orientate = 1)$Estimate[1]
    theta_corrected_mv_lasso_ivw_DP[i,]=mvmr_lasso(bx=X_DP_ivw,by=Y_DP_ivw,seby = seY_ivw)$th_post[1]
    theta_corrected_mv_cml_ivw_DP[i,]=MVmr_cML(b_exp=X_DP_ivw,b_out =as.matrix(Y_DP_ivw),se_bx =seX_mv_ivw,Sig_inv_l = SIGinv_mv_ivw,n=sample_size,K_vec = c(0:(N_non_null-d-1)),random_start =1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03 )$BIC_theta[1]
  }
  
  theta_corrected_mv_cml_ivw_se=sd(theta_corrected_mv_cml_ivw_DP)
  theta_corrected_mv_ivw_ivw_se=sd(theta_corrected_mv_ivw_ivw_DP)
  theta_corrected_mv_lasso_ivw_se=sd(theta_corrected_mv_lasso_ivw_DP)
  theta_corrected_mv_egger_ivw_se=sd(theta_corrected_mv_egger_ivw_DP)
  theta_corrected_mv_median_ivw_se=sd(theta_corrected_mv_median_ivw_DP)
  
  
}






object_X_adj=mr_mvinput(bx=beta_exp_adj,bxse=se_exp_adj,by=as.vector(beta_Y_no_adj),byse = as.vector(se_Y_no_adj))
object_uv_X_adj=mr_mvinput(bx=beta_X_adj,bxse=se_X_adj,by=as.vector(beta_Y_no_adj),byse = as.vector(se_Y_no_adj))
object_uv_X_adj1=mr_input(bx=as.vector(beta_X_adj),bxse=as.vector(se_X_adj),by=as.vector(beta_Y_no_adj),byse = as.vector(se_Y_no_adj))




est=MVmr_cML(b_exp = beta_exp_adj,b_out = matrix(beta_Y_no_adj,nco=1),se_bx = se_exp_adj,Sig_inv_l = SIG_inv_X_adjust,n=sample_size,K_vec = (0:(N_non_null-d-1)),random_start = 4,min_theta_range = -0.5,max_theta_range = 0.5,maxit = 500,thres = 0.03)
theta_all=est$BIC_theta
valid=setdiff(c(1:N_non_null),est$BIC_invalid)
se_est=MVcML_SdTheta(b_exp =beta_exp_adj,b_out =matrix(beta_Y_no_adj,nco=1),Sig_inv_l = SIG_inv_X_adjust,theta = theta_all,zero_ind = valid,r_vec =NULL   )
theta_mv_cml_X_adj=theta_all[1]
se_theta_mv_cml_X_adj_raw=se_est[1]





est=MVmr_cML(b_exp = beta_X_adj,b_out = matrix(beta_Y_no_adj,nco=1),se_bx = se_X_adj,Sig_inv_l = SIG_inv_uv_X_adjust,n=sample_size,K_vec = (0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,maxit = 500,thres = 0.03)
theta_all=est$BIC_theta
valid=setdiff(c(1:N_non_null),est$BIC_invalid)
se_est=MVcML_SdTheta(b_exp =beta_X_adj,b_out =matrix(beta_Y_no_adj,nco=1),Sig_inv_l = SIG_inv_uv_X_adjust,theta = theta_all,zero_ind = valid,r_vec =NULL   )
theta_uv_cml_X_adj=theta_all[1]
se_theta_uv_cml_X_adj_raw=se_est[1]





est=mr_mvivw(object = object_X_adj)
theta_mv_ivw_X_adj=est$Estimate[1]
se_theta_mv_ivw_X_adj_raw=est$StdError[1]


est=mr_ivw(object = object_uv_X_adj1)
theta_uv_ivw_X_adj=est$Estimate[1]
se_theta_uv_ivw_X_adj_raw=est$StdError[1]




est=mvmr_lasso(bx=beta_exp_adj,by=matrix(beta_Y_no_adj,ncol=1),seby = as.vector(se_Y_no_adj))
theta_mv_lasso_X_adj=est$th_post[1]
se_theta_mv_lasso_X_adj_raw=est$se_post[1]




est=mvmr_lasso(bx=beta_X_adj,by=matrix(beta_Y_no_adj,ncol=1),seby = as.vector(se_Y_no_adj))
theta_uv_lasso_X_adj=est$th_post[1]
se_theta_uv_lasso_X_adj_raw=est$se_post[1]






est=mr_mvegger(object = object_X_adj,orientate = 1)
theta_mv_egger_X_adj=est$Estimate[1]
se_theta_mv_egger_X_adj_raw=est$StdError.Est[1]




est=mr_mvegger(object = object_uv_X_adj,orientate = 1)
theta_uv_egger_X_adj=est$Estimate[1]
se_theta_uv_egger_X_adj_raw=est$StdError.Est[1]





est=mr_mvmedian(object = object_X_adj,iterations = 100)
theta_mv_median_X_adj=est$Estimate[1]
se_theta_mv_median_X_adj_raw=est$StdError[1]



est=mr_median(object = object_uv_X_adj1,iterations = 100)
theta_uv_median_X_adj=est$Estimate[1]
se_theta_uv_median_X_adj_raw=est$StdError[1]






est_mv_ivw_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_mv_egger_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_mv_cml_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_mv_lasso_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_mv_median_X_adj_DP=matrix(0,nrow=DP,ncol=1)

est_uv_ivw_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_uv_egger_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_uv_cml_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_uv_lasso_X_adj_DP=matrix(0,nrow=DP,ncol=1)
est_uv_median_X_adj_DP=matrix(0,nrow=DP,ncol=1)


for(i in 1:DP){
   
  data_uv_X_adj=matrix(0,nrow=N_non_null,ncol=2)
  
   
  data_mv_X_adj=matrix(0,nrow=N_non_null,ncol=d+2)
  
  for(j in 1:N_non_null){
    
    data_uv_X_adj[j,]=mvrnorm(n=1,mu=c(beta_X_adj[j,],beta_Y_no_adj[j,]),Sigma = SIG_uv_X_adjust[[j]])
   
    
     
    data_mv_X_adj[j,]=mvrnorm(n=1,mu=c(beta_exp_adj[j,],beta_Y_no_adj[j,]),Sigma = SIG_X_adjust[[j]])
     
  }
  
 
  object_X_adj=mr_mvinput(bx=data_mv_X_adj[,1:(ncol(data_mv_X_adj)-1)],bxse=se_exp_adj,by=data_mv_X_adj[,ncol(data_mv_X_adj)],byse = as.vector(se_Y_no_adj))
  object_uv_X_adj=mr_mvinput(bx=as.matrix(data_uv_X_adj[,1:(ncol(data_uv_X_adj)-1)],nrow=N_non_null,ncol=1),bxse=se_X_adj,by=data_uv_X_adj[,ncol(data_uv_X_adj)] ,byse = as.vector(se_Y_no_adj))
   
  object_uv_X_adj1=mr_input(bx= data_uv_X_adj[,1:(ncol(data_uv_X_adj)-1)] ,bxse=as.vector(se_X_adj),by=data_uv_X_adj[,ncol(data_uv_X_adj)] ,byse = as.vector(se_Y_no_adj))
     
  
  est_mv_lasso_X_adj=mvmr_lasso(bx=data_mv_X_adj[,1:(ncol(data_mv_X_adj)-1)],by=data_mv_X_adj[,ncol(data_mv_X_adj)],seby = as.vector(se_Y_no_adj))
  est_uv_lasso_X_adj=mvmr_lasso(bx=as.matrix(data_uv_X_adj[,1:(ncol(data_uv_X_adj)-1)],nrow=N_non_null,ncol=1),by=data_uv_X_adj[,ncol(data_uv_X_adj)] ,seby = as.vector(se_Y_no_adj))
 
  est_mv_cml_X_adj=MVmr_cML(b_exp =data_mv_X_adj[,1:(ncol(data_mv_X_adj)-1)],b_out =as.matrix(data_mv_X_adj[,ncol(data_mv_X_adj)],nrow=N_non_null,ncol=1),se_bx =se_exp_adj,Sig_inv_l = SIG_inv_X_adjust,n=sample_size,K_vec =c(0:(N_non_null-d-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03  )
  est_uv_cml_X_adj=MVmr_cML(b_exp =as.matrix(data_uv_X_adj[,1:(ncol(data_uv_X_adj)-1)],nrow=N_non_null,ncol=1),b_out =as.matrix(data_uv_X_adj[,ncol(data_uv_X_adj)],nrow=N_non_null,ncol=1),se_bx =se_exp_adj,Sig_inv_l = SIG_inv_uv_X_adjust,n=sample_size,K_vec =c(0:(N_non_null-2)),random_start = 1,min_theta_range = -0.5,max_theta_range = 0.5,thres = 0.03   )
    
  est_mv_ivw_X_adj=mr_mvivw(object = object_X_adj)
  est_uv_ivw_X_adj=mr_mvivw(object = object_uv_X_adj)
    
  
  
  
  est_mv_egger_X_adj=mr_mvegger(object = object_X_adj,orientate = 1)
  est_uv_egger_X_adj=mr_mvegger(object = object_uv_X_adj,orientate = 1)
    
  
  est_mv_median_X_adj=mr_mvmedian(object = object_X_adj,iterations = 200 )
  est_uv_median_X_adj=mr_median(object = object_uv_X_adj1,iterations = 200 )
     
  
 
  est_mv_ivw_X_adj_DP[i,]=est_mv_ivw_X_adj$Estimate[1]
  est_mv_egger_X_adj_DP[i,]=est_mv_egger_X_adj$Estimate[1]
  est_mv_cml_X_adj_DP[i,]=est_mv_cml_X_adj$BIC_theta[1]
  est_mv_lasso_X_adj_DP[i,]=est_mv_lasso_X_adj$th_post[1]
  est_mv_median_X_adj_DP[i,]=est_mv_median_X_adj$Estimate[1]
  est_uv_ivw_X_adj_DP[i,]=est_uv_ivw_X_adj$Estimate[1]
  est_uv_egger_X_adj_DP[i,]=est_uv_egger_X_adj$Estimate[1]
  est_uv_cml_X_adj_DP[i,]=est_uv_cml_X_adj$BIC_theta[1]
  est_uv_lasso_X_adj_DP[i,]=est_uv_lasso_X_adj$th_post[1]
  est_uv_median_X_adj_DP[i,]=est_uv_median_X_adj$Estimate[1]
  
  print(i)
  
  
}

 
 
se_theta_mv_cml_X_adj=sd(est_mv_cml_X_adj_DP)
 
se_theta_mv_ivw_X_adj=sd(est_mv_ivw_X_adj_DP)
 
se_theta_mv_lasso_X_adj=sd(est_mv_lasso_X_adj_DP)
 
se_theta_mv_median_X_adj=sd(est_mv_median_X_adj_DP)
 
se_theta_mv_egger_X_adj=sd(est_mv_egger_X_adj_DP)
 
 
se_theta_uv_cml_X_adj=sd(est_uv_cml_X_adj_DP)
 
se_theta_uv_ivw_X_adj=sd(est_uv_ivw_X_adj_DP)
 
se_theta_uv_lasso_X_adj=sd(est_uv_lasso_X_adj_DP)
 
se_theta_uv_median_X_adj=sd(est_uv_median_X_adj_DP)
 
se_theta_uv_egger_X_adj=sd(est_uv_egger_X_adj_DP)
 


  

 

if(d==1){
  result=data.frame( 
                    theta_mv_cml_X_adj,se_theta_mv_cml_X_adj,theta_uv_cml_X_adj,se_theta_uv_cml_X_adj,
                  
                    theta_mv_egger_X_adj,se_theta_mv_egger_X_adj,theta_uv_egger_X_adj,se_theta_uv_egger_X_adj,
                    theta_mv_lasso_X_adj,se_theta_mv_lasso_X_adj,theta_uv_lasso_X_adj,se_theta_uv_lasso_X_adj,
                    theta_mv_median_X_adj,se_theta_mv_median_X_adj,theta_uv_median_X_adj,se_theta_uv_median_X_adj,
                  
                    theta_mv_ivw_X_adj,se_theta_mv_ivw_X_adj,theta_uv_ivw_X_adj,se_theta_uv_ivw_X_adj,
                    
                    theta_corrected_uv_cml_ivw=theta_corrected_uv_cml_ivw,theta_corrected_uv_cml_ivw_se=theta_corrected_uv_cml_ivw_se,
                    theta_corrected_uv_ivw_ivw=theta_corrected_uv_ivw_ivw,theta_corrected_uv_ivw_ivw_se=theta_corrected_uv_ivw_ivw_se,
                    theta_corrected_uv_lasso_ivw=theta_corrected_uv_lasso_ivw,theta_corrected_uv_lasso_ivw_se=theta_corrected_uv_lasso_ivw_se,
                    theta_corrected_uv_median_ivw=theta_corrected_uv_median_ivw,theta_corrected_uv_median_ivw_se=theta_corrected_uv_median_ivw_se,
                    theta_corrected_uv_egger_ivw=theta_corrected_uv_egger_ivw,theta_corrected_uv_egger_ivw_se=theta_corrected_uv_egger_ivw_se,
                    theta_corrected_uv_cml_dho=theta_corrected_uv_cml_dho,theta_corrected_uv_cml_dho_se=theta_corrected_uv_cml_dho_se,
                    theta_corrected_uv_ivw_dho=theta_corrected_uv_ivw_dho,theta_corrected_uv_ivw_dho_se=theta_corrected_uv_ivw_dho_se,
                    theta_corrected_uv_lasso_dho=theta_corrected_uv_lasso_dho,theta_corrected_uv_lasso_dho_se=theta_corrected_uv_lasso_dho_se,
                    theta_corrected_uv_median_dho=theta_corrected_uv_median_dho,theta_corrected_uv_median_dho_se=theta_corrected_uv_median_dho_se,
                    theta_corrected_uv_egger_dho=theta_corrected_uv_egger_dho,theta_corrected_uv_egger_dho_se=theta_corrected_uv_egger_dho_se,
                    theta_corrected_uv_cml_sh=theta_corrected_uv_cml_sh,theta_corrected_uv_cml_sh_se=theta_corrected_uv_cml_sh_se,
                    theta_corrected_uv_ivw_sh=theta_corrected_uv_ivw_sh,theta_corrected_uv_ivw_sh_se=theta_corrected_uv_ivw_sh_se,
                    theta_corrected_uv_lasso_sh=theta_corrected_uv_lasso_sh,theta_corrected_uv_lasso_sh_se=theta_corrected_uv_lasso_sh_se,
                    theta_corrected_uv_median_sh=theta_corrected_uv_median_sh,theta_corrected_uv_median_sh_se=theta_corrected_uv_median_sh_se,
                    theta_corrected_uv_egger_sh=theta_corrected_uv_egger_sh,theta_corrected_uv_egger_sh_se=theta_corrected_uv_egger_sh_se,
                    theta_corrected_uv_cml_cml=theta_corrected_uv_cml_cml,theta_corrected_uv_cml_cml_se=theta_corrected_uv_cml_cml_se,
                    theta_corrected_uv_ivw_cml=theta_corrected_uv_ivw_cml,theta_corrected_uv_ivw_cml_se=theta_corrected_uv_ivw_cml_se,
                    theta_corrected_uv_lasso_cml=theta_corrected_uv_lasso_cml,theta_corrected_uv_lasso_cml_se=theta_corrected_uv_lasso_cml_se,
                    theta_corrected_uv_median_cml=theta_corrected_uv_median_cml,theta_corrected_uv_median_cml_se=theta_corrected_uv_median_cml_se,
                    theta_corrected_uv_egger_cml=theta_corrected_uv_egger_cml,theta_corrected_uv_egger_cml_se=theta_corrected_uv_egger_cml_se,
                    theta_corrected_mv_cml_ivw=theta_corrected_mv_cml_ivw,theta_corrected_mv_cml_ivw_se=theta_corrected_mv_cml_ivw_se,
                    theta_corrected_mv_ivw_ivw=theta_corrected_mv_ivw_ivw,theta_corrected_mv_ivw_ivw_se=theta_corrected_mv_ivw_ivw_se,
                    theta_corrected_mv_lasso_ivw=theta_corrected_mv_lasso_ivw,theta_corrected_mv_lasso_ivw_se=theta_corrected_mv_lasso_ivw_se,
                    theta_corrected_mv_median_ivw=theta_corrected_mv_median_ivw,theta_corrected_mv_median_ivw_se=theta_corrected_mv_median_ivw_se,
                    theta_corrected_mv_egger_ivw=theta_corrected_mv_egger_ivw,theta_corrected_mv_egger_ivw_se=theta_corrected_mv_egger_ivw_se,
                    theta_corrected_mv_cml_dho=theta_corrected_mv_cml_dho,theta_corrected_mv_cml_dho_se=theta_corrected_mv_cml_dho_se,
                    theta_corrected_mv_ivw_dho=theta_corrected_mv_ivw_dho,theta_corrected_mv_ivw_dho_se=theta_corrected_mv_ivw_dho_se,
                    theta_corrected_mv_lasso_dho=theta_corrected_mv_lasso_dho,theta_corrected_mv_lasso_dho_se=theta_corrected_mv_lasso_dho_se,
                    theta_corrected_mv_median_dho=theta_corrected_mv_median_dho,theta_corrected_mv_median_dho_se=theta_corrected_mv_median_dho_se,
                    theta_corrected_mv_egger_dho=theta_corrected_mv_egger_dho,theta_corrected_mv_egger_dho_se=theta_corrected_mv_egger_dho_se,
                    theta_corrected_mv_cml_sh=theta_corrected_mv_cml_sh,theta_corrected_mv_cml_sh_se=theta_corrected_mv_cml_sh_se,
                    theta_corrected_mv_ivw_sh=theta_corrected_mv_ivw_sh,theta_corrected_mv_ivw_sh_se=theta_corrected_mv_ivw_sh_se,
                    theta_corrected_mv_lasso_sh=theta_corrected_mv_lasso_sh,theta_corrected_mv_lasso_sh_se=theta_corrected_mv_lasso_sh_se,
                    theta_corrected_mv_median_sh=theta_corrected_mv_median_sh,theta_corrected_mv_median_sh_se=theta_corrected_mv_median_sh_se,
                    theta_corrected_mv_egger_sh=theta_corrected_mv_egger_sh,theta_corrected_mv_egger_sh_se=theta_corrected_mv_egger_sh_se,
                    theta_corrected_mv_cml_cml=theta_corrected_mv_cml_cml,theta_corrected_mv_cml_cml_se=theta_corrected_mv_cml_cml_se,
                    theta_corrected_mv_ivw_cml=theta_corrected_mv_ivw_cml,theta_corrected_mv_ivw_cml_se=theta_corrected_mv_ivw_cml_se,
                    theta_corrected_mv_lasso_cml=theta_corrected_mv_lasso_cml,theta_corrected_mv_lasso_cml_se=theta_corrected_mv_lasso_cml_se,
                    theta_corrected_mv_median_cml=theta_corrected_mv_median_cml,theta_corrected_mv_median_cml_se=theta_corrected_mv_median_cml_se,
                    theta_corrected_mv_egger_cml=theta_corrected_mv_egger_cml,theta_corrected_mv_egger_cml_se=theta_corrected_mv_egger_cml_se,
                    se_theta_mv_cml_X_adj_raw,
                    se_theta_mv_egger_X_adj_raw,
                    se_theta_mv_lasso_X_adj_raw,
                    se_theta_mv_ivw_X_adj_raw,
                    se_theta_mv_median_X_adj_raw,
                    se_theta_uv_cml_X_adj_raw,
                    se_theta_uv_egger_X_adj_raw,
                    se_theta_uv_lasso_X_adj_raw,
                    se_theta_uv_ivw_X_adj_raw,
                    se_theta_uv_median_X_adj_raw,
                    theta_corrected_uv_cml_cml_se_raw,
                    theta_corrected_uv_ivw_cml_se_raw,
                    theta_corrected_uv_egger_cml_se_raw,
                    theta_corrected_uv_lasso_cml_se_raw,
                    theta_corrected_uv_median_cml_se_raw,
                    theta_corrected_uv_cml_ivw_se_raw,
                    theta_corrected_uv_ivw_ivw_se_raw,
                    theta_corrected_uv_egger_ivw_se_raw,
                    theta_corrected_uv_lasso_ivw_se_raw,
                    theta_corrected_uv_median_ivw_se_raw,
                    theta_corrected_uv_cml_dho_se_raw,
                    theta_corrected_uv_ivw_dho_se_raw,
                    theta_corrected_uv_egger_dho_se_raw,
                    theta_corrected_uv_lasso_dho_se_raw,
                    theta_corrected_uv_median_dho_se_raw,
                    theta_corrected_uv_cml_sh_se_raw,
                    theta_corrected_uv_ivw_sh_se_raw,
                    theta_corrected_uv_egger_sh_se_raw,
                    theta_corrected_uv_lasso_sh_se_raw,
                    theta_corrected_uv_median_sh_se_raw,
                    theta_corrected_mv_cml_cml_se_raw,
                    theta_corrected_mv_ivw_cml_se_raw,
                    theta_corrected_mv_egger_cml_se_raw,
                    theta_corrected_mv_lasso_cml_se_raw,
                    theta_corrected_mv_median_cml_se_raw,
                    theta_corrected_mv_cml_ivw_se_raw,
                    theta_corrected_mv_ivw_ivw_se_raw,
                    theta_corrected_mv_egger_ivw_se_raw,
                    theta_corrected_mv_lasso_ivw_se_raw,
                    theta_corrected_mv_median_ivw_se_raw,
                    theta_corrected_mv_cml_dho_se_raw,
                    theta_corrected_mv_ivw_dho_se_raw,
                    theta_corrected_mv_egger_dho_se_raw,
                    theta_corrected_mv_lasso_dho_se_raw,
                    theta_corrected_mv_median_dho_se_raw,
                    theta_corrected_mv_cml_sh_se_raw,
                    theta_corrected_mv_ivw_sh_se_raw,
                    theta_corrected_mv_egger_sh_se_raw,
                    theta_corrected_mv_lasso_sh_se_raw,
                    theta_corrected_mv_median_sh_se_raw)
                     
  
  
}
if(d>1){
  result=data.frame(theta_mv_cml_X_adj,se_theta_mv_cml_X_adj,theta_uv_cml_X_adj,se_theta_uv_cml_X_adj,
                    
                    theta_mv_egger_X_adj,se_theta_mv_egger_X_adj,theta_uv_egger_X_adj,se_theta_uv_egger_X_adj,
                    theta_mv_lasso_X_adj,se_theta_mv_lasso_X_adj,theta_uv_lasso_X_adj,se_theta_uv_lasso_X_adj,
                    theta_mv_median_X_adj,se_theta_mv_median_X_adj,theta_uv_median_X_adj,se_theta_uv_median_X_adj,
                    
                    theta_mv_ivw_X_adj,se_theta_mv_ivw_X_adj,theta_uv_ivw_X_adj,se_theta_uv_ivw_X_adj,
                    
                    theta_corrected_uv_cml_ivw=theta_corrected_uv_cml_ivw,theta_corrected_uv_cml_ivw_se=theta_corrected_uv_cml_ivw_se,
                    theta_corrected_uv_ivw_ivw=theta_corrected_uv_ivw_ivw,theta_corrected_uv_ivw_ivw_se=theta_corrected_uv_ivw_ivw_se,
                    theta_corrected_uv_lasso_ivw=theta_corrected_uv_lasso_ivw,theta_corrected_uv_lasso_ivw_se=theta_corrected_uv_lasso_ivw_se,
                    theta_corrected_uv_median_ivw=theta_corrected_uv_median_ivw,theta_corrected_uv_median_ivw_se=theta_corrected_uv_median_ivw_se,
                    theta_corrected_uv_egger_ivw=theta_corrected_uv_egger_ivw,theta_corrected_uv_egger_ivw_se=theta_corrected_uv_egger_ivw_se,
                    theta_corrected_uv_cml_cml=theta_corrected_uv_cml_cml,theta_corrected_uv_cml_cml_se=theta_corrected_uv_cml_cml_se,
                    theta_corrected_uv_ivw_cml=theta_corrected_uv_ivw_cml,theta_corrected_uv_ivw_cml_se=theta_corrected_uv_ivw_cml_se,
                    theta_corrected_uv_lasso_cml=theta_corrected_uv_lasso_cml,theta_corrected_uv_lasso_cml_se=theta_corrected_uv_lasso_cml_se,
                    theta_corrected_uv_median_cml=theta_corrected_uv_median_cml,theta_corrected_uv_median_cml_se=theta_corrected_uv_median_cml_se,
                    theta_corrected_uv_egger_cml=theta_corrected_uv_egger_cml,theta_corrected_uv_egger_cml_se=theta_corrected_uv_egger_cml_se,
                    theta_corrected_mv_cml_ivw=theta_corrected_mv_cml_ivw,theta_corrected_mv_cml_ivw_se=theta_corrected_mv_cml_ivw_se,
                    theta_corrected_mv_ivw_ivw=theta_corrected_mv_ivw_ivw,theta_corrected_mv_ivw_ivw_se=theta_corrected_mv_ivw_ivw_se,
                    theta_corrected_mv_lasso_ivw=theta_corrected_mv_lasso_ivw,theta_corrected_mv_lasso_ivw_se=theta_corrected_mv_lasso_ivw_se,
                    theta_corrected_mv_median_ivw=theta_corrected_mv_median_ivw,theta_corrected_mv_median_ivw_se=theta_corrected_mv_median_ivw_se,
                    theta_corrected_mv_egger_ivw=theta_corrected_mv_egger_ivw,theta_corrected_mv_egger_ivw_se=theta_corrected_mv_egger_ivw_se,
                    theta_corrected_mv_cml_cml=theta_corrected_mv_cml_cml,theta_corrected_mv_cml_cml_se=theta_corrected_mv_cml_cml_se,
                    theta_corrected_mv_ivw_cml=theta_corrected_mv_ivw_cml,theta_corrected_mv_ivw_cml_se=theta_corrected_mv_ivw_cml_se,
                    theta_corrected_mv_lasso_cml=theta_corrected_mv_lasso_cml,theta_corrected_mv_lasso_cml_se=theta_corrected_mv_lasso_cml_se,
                    theta_corrected_mv_median_cml=theta_corrected_mv_median_cml,theta_corrected_mv_median_cml_se=theta_corrected_mv_median_cml_se,
                    theta_corrected_mv_egger_cml=theta_corrected_mv_egger_cml,theta_corrected_mv_egger_cml_se=theta_corrected_mv_egger_cml_se,
                    se_theta_mv_cml_X_adj_raw,
                    se_theta_mv_egger_X_adj_raw,
                    se_theta_mv_lasso_X_adj_raw,
                    se_theta_mv_ivw_X_adj_raw,
                    se_theta_mv_median_X_adj_raw,
                    se_theta_uv_cml_X_adj_raw,
                    se_theta_uv_egger_X_adj_raw,
                    se_theta_uv_lasso_X_adj_raw,
                    se_theta_uv_ivw_X_adj_raw,
                    se_theta_uv_median_X_adj_raw,
                    theta_corrected_uv_cml_cml_se_raw,
                    theta_corrected_uv_ivw_cml_se_raw,
                    theta_corrected_uv_egger_cml_se_raw,
                    theta_corrected_uv_lasso_cml_se_raw,
                    theta_corrected_uv_median_cml_se_raw,
                    theta_corrected_uv_cml_ivw_se_raw,
                    theta_corrected_uv_ivw_ivw_se_raw,
                    theta_corrected_uv_egger_ivw_se_raw,
                    theta_corrected_uv_lasso_ivw_se_raw,
                    theta_corrected_uv_median_ivw_se_raw,
                    theta_corrected_mv_cml_cml_se_raw,
                    theta_corrected_mv_ivw_cml_se_raw,
                    theta_corrected_mv_egger_cml_se_raw,
                    theta_corrected_mv_lasso_cml_se_raw,
                    theta_corrected_mv_median_cml_se_raw,
                    theta_corrected_mv_cml_ivw_se_raw,
                    theta_corrected_mv_ivw_ivw_se_raw,
                    theta_corrected_mv_egger_ivw_se_raw,
                    theta_corrected_mv_lasso_ivw_se_raw,
                    theta_corrected_mv_median_ivw_se_raw)
                 
}     




result=cbind(result, btotal)




    }
    
    output=list(total_estimate=total_estimate)
    save(output,file=output_file[i1,i2])
  }
 
  
  
}

