setwd("D:/paper2/real_data/BMI_BP")

d=c(1,2,5)
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

DP=100
library(MendelianRandomization)
cov_name=list()
exp_name=numeric()

for(i in 1:length(d)){
  cov_name[[i]]=numeric()
  exp_name[i]=paste("GWAS_",d[i],"MetPCNoGenadj.BMI.glm.linear",sep="")
}
for(i in 1:length(d)){
  for(j in 1:d[i]){
    cov_name[[i]][j]=paste("GWAS.metPCnoGen",j,".glm.linear",sep="")
  }
}
exp=list()
cov=list()
for(i in 1:length(d)){
  cov[[i]]=list()
}
names1=c("chr","pos","snp","ref","alt","a1","test","n","beta","se","z","p")
for(i in 1:length(d)){
  exp[[i]]=read.table(file=exp_name[i],sep="\t",header = FALSE)
  colnames(exp[[i]])=names1
}
for(i in 1:length(d)){
  for(j in 1:d[i]){
    cov[[i]][[j]]=read.table(file=cov_name[[i]][j],sep="\t",header = FALSE)
    colnames(cov[[i]][[j]])=names1
  }
}
SBP=read.table(file="GWAS.SBP.glm.linear",header=FALSE,sep="\t")

colnames(SBP)=names1
snp=SBP$snp
for(i in 1:length(d)){
  snp=intersect(snp,exp[[i]]$snp)
  for(j in 1:d[[i]]){
    snp=intersect(snp,cov[[i]][[j]]$snp)
  }
}
SBP=SBP[SBP$snp%in%snp,]
for(i in 1:length(d)){
  exp[[i]]=exp[[i]][exp[[i]]$snp%in%snp,]
  for(j in 1:d[[i]]){
    cov[[i]][[j]]=cov[[i]][[j]][cov[[i]][[j]]$snp%in%snp,]
  }
}
indsnp=read.table(file="110K_QCed0.001.bim",header=FALSE,sep="\t")[,2]
IV=list()
IVuv=list()

for(i in 1:length(d)){
  IVuv[[i]]=exp[[i]]$snp[exp[[i]]$snp%in%indsnp&exp[[i]]$p<5e-4]
  IV[[i]]=c(exp[[i]]$snp[exp[[i]]$snp%in%indsnp&exp[[i]]$p<5e-4])
  for(j in 1:d[[i]]){
    IV[[i]]=c(IV[[i]],cov[[i]][[j]]$snp[cov[[i]][[j]]$snp%in%indsnp&cov[[i]][[j]]$p<5e-4])
  }
  IVuv[[i]]=unique(IVuv[[i]])
  IV[[i]]=unique(IV[[i]])
  
}
beta_exp=list()
beta_exp_uv=list()
se_exp=list()
se_exp_uv=list()
beta_out=list()
beta_out_uv=list()
se_out=list()
se_out_uv=list()
SIG=list()
SIGexp=list()
SIGinv=list()
SIG_uv=list()
SIGinv_uv=list()
for(i in 1:length(d)){
  beta_exp[[i]]=matrix(0,nrow=length(IV[[i]]),ncol =d[i]+1)
  se_exp[[i]]=matrix(0,nrow=length(IV[[i]]),ncol =d[i]+1)
  SIG[[i]]=list()
  SIGexp[[i]]=list()
  SIGinv[[i]]=list()
  SIG_uv[[i]]=list()
  SIGinv_uv[[i]]=list()
}
for(i in 1:length(d)){
  beta_exp[[i]][,1]=exp[[i]]$beta[exp[[i]]$snp%in%IV[[i]]]
  for(j in 2:(d[i]+1)){
    beta_exp[[i]][,j]=cov[[i]][[j-1]]$beta[cov[[i]][[j-1]]$snp%in%IV[[i]]]
  }
  se_exp[[i]][,1]=exp[[i]]$se[exp[[i]]$snp%in%IV[[i]]]
  for(j in 2:(d[i]+1)){
    se_exp[[i]][,j]=cov[[i]][[j-1]]$se[cov[[i]][[j-1]]$snp%in%IV[[i]]]
  }
  beta_exp_uv[[i]]=exp[[i]]$beta[exp[[i]]$snp%in%IVuv[[i]]]
  se_exp_uv[[i]]=exp[[i]]$se[exp[[i]]$snp%in%IVuv[[i]]]
  beta_out_uv[[i]]=SBP$beta[SBP$snp%in%IVuv[[i]]]
  se_out_uv[[i]]=SBP$se[SBP$snp%in%IVuv[[i]]]
  beta_out[[i]]=SBP$beta[SBP$snp%in%IV[[i]]]
  se_out[[i]]=SBP$se[SBP$snp%in%IV[[i]]]
}
correlation=list()
correlation_uv=list()
null_z=(SBP$p>0.1)&(SBP$snp%in%indsnp)
for(i in 1:length(d)){
  null_z=null_z&(exp[[i]]$p>0.1)&(exp[[i]]$snp%in%indsnp)
  for(j in 1:d[[i]]){
    null_z=null_z&(cov[[i]][[j]]$p>0.1)&(cov[[i]][[j]]$snp%in%indsnp)
  }
}
nullz=list()
nullz_uv=list()
for(i in 1:length(d)){
  nullz[[i]]=cbind(exp[[i]]$z[null_z])
  for(j in 1:d[i]){
    nullz[[i]]=cbind(nullz[[i]],cov[[i]][[j]]$z[null_z])
  }
  nullz[[i]]=cbind(nullz[[i]],SBP$z[null_z])
  nullz_uv[[i]]=cbind(exp[[i]]$z[null_z],SBP$z[null_z])
}
for(i in 1:length(d)){
  correlation[[i]]=cor(nullz[[i]])
  correlation_uv[[i]]=cor(nullz_uv[[i]])
}
for(i in 1:length(d)){
  for(j in 1:length(IV[[i]])){
    SIG[[i]][[j]]=diag(c(se_exp[[i]][j,],se_out[[i]][j]),nrow=d[i]+2,ncol = d[i]+2)%*%correlation[[i]]%*%diag(c(se_exp[[i]][j,],se_out[[i]][j]),nrow=d[i]+2,ncol = d[i]+2)
    SIGinv[[i]][[j]]=solve(SIG[[i]][[j]])
    SIGexp[[i]][[j]]=diag(c(se_exp[[i]][j,] ),nrow=d[i]+1,ncol = d[i]+1)%*%correlation[[i]][1:(d[i]+1),1:(d[i]+1)]%*%diag(c(se_exp[[i]][j,] ),nrow=d[i]+1,ncol = d[i]+1)
  }
  for(j in 1:length(IVuv[[i]])){
    SIG_uv[[i]][[j]]=diag(c(se_exp_uv[[i]][j],se_out_uv[[i]][j]),nrow=2,ncol = 2)%*%correlation_uv[[i]]%*%diag(c(se_exp_uv[[i]][j],se_out_uv[[i]][j]),nrow=2,ncol =2)
    SIGinv_uv[[i]][[j]]=solve(SIG_uv[[i]][[j]])
  }
  
}
theta_cml=numeric()
theta_cmluv=numeric()
theta_egger=numeric()
theta_eggeruv=numeric()
theta_lasso=numeric()
theta_lassouv=numeric()
theta_median=numeric()
theta_medianuv=numeric()
theta_ivw=numeric()
theta_ivwuv=numeric()
se_cml=numeric()
se_cmluv=numeric()
se_egger=numeric()
se_eggeruv=numeric()
se_lasso=numeric()
se_lassouv=numeric()
se_median=numeric()
se_medianuv=numeric()
se_ivw=numeric()
se_ivwuv=numeric()


theta_cmlDP=numeric()
theta_cmluvDP=numeric()
theta_eggerDP=numeric()
theta_eggeruvDP=numeric()
theta_lassoDP=numeric()
theta_lassouvDP=numeric()
theta_medianDP=numeric()
theta_medianuvDP=numeric()
theta_ivwDP=numeric()
theta_ivwuvDP=numeric()

for(i1 in 1:length(d)){
  dDP=d[i1]
  
  niv=length(IV[[i1]])
  niv_uv=length(IVuv[[i1]])
  SIGDP=SIG[[i1]]
  SIGinvDP=SIGinv[[i1]]
  SIGDP_uv=SIG_uv[[i1]]
  SIGinvDP_uv=SIGinv_uv[[i1]]
  se_expDP=se_exp[[i1]]
  se_expDP_uv=se_exp_uv[[i1]]
  se_outDP=se_out[[i1]]
  se_outDP_uv=se_out_uv[[i1]]
  beta_exp_d=beta_exp[[i1]]
  beta_exp_d_uv=beta_exp_uv[[i1]]
  beta_out_d=beta_out[[i1]]
  beta_out_d_uv=beta_out_uv[[i1]]
  for(i2 in 1:DP){
    beta_expDP=matrix(0,nrow = niv,ncol = dDP+1)
    beta_outDP=matrix(0,nrow = niv,ncol =1)
    beta_expDP_uv=matrix(0,nrow = niv_uv,ncol =1)
    beta_outDP_uv=matrix(0,nrow = niv_uv,ncol =1)
    for(i in 1:niv){
      beta=mvrnorm(n=1,mu=c(beta_exp_d[i,],beta_out_d[i]),Sigma = SIGDP[[i]])
      beta_expDP[i,]=beta[1:(length(beta)-1)]
      beta_outDP[i,]=beta[length(beta)]
      if(beta_expDP[i,1]<0){
        beta_expDP[i,]=-beta_expDP[i,]
        beta_outDP[i,]=-beta_outDP[i,]
      }
    }
    for(i in 1:niv_uv){
      beta_uv=mvrnorm(n=1,mu=c(beta_exp_d_uv[i],beta_out_d_uv[i]),Sigma = SIGDP_uv[[i]])
      beta_expDP_uv[i,]=beta_uv[1:(length(beta_uv)-1)]
      beta_outDP_uv[i,]=beta_uv[length(beta_uv)]
      if(beta_expDP_uv[i,1]<0){
        beta_expDP_uv[i,]=-beta_expDP_uv[i,]
        beta_outDP_uv[i,]=-beta_outDP_uv[i,]
      }
    }
    obj=mr_mvinput(bx=beta_expDP,bxse = se_expDP,by=beta_outDP[,1],byse = se_outDP)
    obj_uv=mr_input(bx=beta_expDP_uv[,1],bxse = se_expDP_uv,by=beta_outDP_uv[,1],byse = se_outDP_uv)
    theta_cmlDP[i2]=MVmr_cML(b_exp =beta_expDP,b_out=beta_outDP,se_bx = se_expDP,Sig_inv_l = SIGinvDP,
                             n=100000,K_vec = c(0:(niv-dDP-2)),random_start = 1,min_theta_range = 0.1,
                             max_theta_range = 0.1,maxit = 500,thres = 0.01)$BIC_theta[1]
    theta_cmluvDP[i2]=MVmr_cML(b_exp =beta_expDP_uv,b_out=matrix(beta_outDP_uv,ncol=1,nrow=length(beta_outDP_uv)),se_bx = matrix(se_expDP_uv,ncol = 1,nrow=length(se_expDP_uv)),Sig_inv_l = SIGinvDP_uv,
                               n=100000,K_vec = c(0:(niv_uv-2)),random_start = 1,min_theta_range = 0.1,
                               max_theta_range = 0.1,maxit = 500,thres = 0.01)$BIC_theta[1]
    theta_eggerDP[i2]=mr_mvegger(object = obj)$Estimate[1]
    theta_eggeruvDP[i2]=mr_egger(object = obj_uv)$Estimate[1]
    theta_lassoDP[i2]=mvmr_lasso(bx=beta_expDP,by=beta_outDP[,1],seby=se_outDP)$th_post[1]
    theta_lassouvDP[i2]=mvmr_lasso(bx=beta_expDP_uv,by=beta_outDP_uv[,1],seby=se_outDP_uv)$th_post[1]
    theta_ivwDP[i2]=mr_mvivw(object = obj)$Estimate[1]
    theta_ivwuvDP[i2]=mr_ivw(object = obj_uv)$Estimate[1]
    theta_medianDP[i2]=mr_mvmedian(object = obj,iterations = 200)$Estimate[1]
    theta_medianuvDP[i2]=mr_median(object = obj_uv,iterations = 200)$Estimate[1]
    print(i2)
    
  }
  theta_cml[i1]=mean(theta_cmlDP)
  theta_cmluv[i1]=mean(theta_cmluvDP)
  theta_egger[i1]=mean(theta_eggerDP)
  theta_eggeruv[i1]=mean(theta_eggeruvDP)
  theta_lasso[i1]=mean(theta_lassoDP)
  theta_lassouv[i1]=mean(theta_lassouvDP)
  theta_ivw[i1]=mean(theta_ivwDP)
  theta_ivwuv[i1]=mean(theta_ivwuvDP)
  theta_median[i1]=mean(theta_medianDP)
  theta_medianuv[i1]=mean(theta_medianuvDP)
  se_cml[i1]=sd(theta_cmlDP)
  se_cmluv[i1]=sd(theta_cmluvDP)
  se_egger[i1]=sd(theta_eggerDP)
  se_eggeruv[i1]=sd(theta_eggeruvDP)
  se_lasso[i1]=sd(theta_lassoDP)
  se_lassouv[i1]=sd(theta_lassouvDP)
  se_ivw[i1]=sd(theta_ivwDP)
  se_ivwuv[i1]=sd(theta_ivwuvDP)
  se_median[i1]=sd(theta_medianDP)
  se_medianuv[i1]=sd(theta_medianuvDP)
}

output=list(data.frame(theta_cml,theta_cmluv,theta_egger,theta_eggeruv,theta_lasso,theta_lassouv,theta_ivw,theta_ivwuv,
                  theta_median,theta_medianuv,se_cml,se_cmluv,se_egger,se_eggeruv,se_lasso,se_lassouv,se_ivw,se_ivwuv,
                  se_median,se_medianuv))

save(output,file="BMI_SBPnoGen.RData")

input=list()
for(i in 1:length(d)){
  input[[i]]=mr_mvinput(bx=beta_exp[[i]],by=beta_out[[i]],bxse=se_exp[[i]],byse = se_out[[i]])
}
for(i in 1:length(d)){
  print(mr_mvivw(object = input[[i]],nx=mean(exp[[i]]$n),correl.x =correlation[[i]][1:(d[i]+1),1:(d[i]+1)]))
}

 