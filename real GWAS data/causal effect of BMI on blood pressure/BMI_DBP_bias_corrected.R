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
library(MASS)
library(glmnet)
exp_name=numeric()

for(i in 1:length(d)){
  exp_name[i]=paste("GWAS_BMI",d[i],"MetPC_corrected.txt",sep="")
}

exp=list()


for(i in 1:length(d)){
  exp[[i]]=read.table(file=exp_name[i],sep="\t",header = TRUE)
  
}
names1=c("chr","pos","snp","ref","alt","a1","test","n","beta","se","z","p")
DBP=read.table(file="GWAS.DBP.glm.linear",header=FALSE,sep="\t")

colnames(DBP)=names1
snp=DBP$snp
for(i in 1:length(d)){
  snp=intersect(snp,exp[[i]]$snp)
}
DBP=DBP[DBP$snp%in%snp,]
for(i in 1:length(d)){
  exp[[i]]=exp[[i]][exp[[i]]$snp%in%snp,]
}
indsnp=read.table(file="110K_QCed0.001.bim",header=FALSE,sep="\t")[,2]

IVuv=list()

for(i in 1:length(d)){
  IVuv[[i]]=exp[[i]]$snp[exp[[i]]$snp%in%indsnp&exp[[i]]$p<5e-6]
  
  IVuv[[i]]=unique(IVuv[[i]])
  
  
}

beta_exp_uv=list()

se_exp_uv=list()

beta_out_uv=list()

se_out_uv=list()


SIG_uv=list()
SIGinv_uv=list()
for(i in 1:length(d)){
  
  
  
  SIG_uv[[i]]=list()
  SIGinv_uv[[i]]=list()
}
for(i in 1:length(d)){
  
  beta_exp_uv[[i]]=exp[[i]]$beta[exp[[i]]$snp%in%IVuv[[i]]]
  se_exp_uv[[i]]=exp[[i]]$se[exp[[i]]$snp%in%IVuv[[i]]]
  beta_out_uv[[i]]=DBP$beta[DBP$snp%in%IVuv[[i]]]
  se_out_uv[[i]]=DBP$se[DBP$snp%in%IVuv[[i]]]
  
}

correlation_uv=list()

null_z=(DBP$p>0.1)&(DBP$snp%in%indsnp)
for(i in 1:length(d)){
  null_z=null_z&(exp[[i]]$p>0.1)&(exp[[i]]$snp%in%indsnp)
  
}

nullz_uv=list()
for(i in 1:length(d)){
  
  nullz_uv[[i]]=cbind(exp[[i]]$z[null_z],DBP$z[null_z])
}
for(i in 1:length(d)){
  
  correlation_uv[[i]]=cor(nullz_uv[[i]])
}
for(i in 1:length(d)){
  
  for(j in 1:length(IVuv[[i]])){
    SIG_uv[[i]][[j]]=diag(c(se_exp_uv[[i]][j],se_out_uv[[i]][j]),nrow=2,ncol = 2)%*%correlation_uv[[i]]%*%diag(c(se_exp_uv[[i]][j],se_out_uv[[i]][j]),nrow=2,ncol =2)
    SIGinv_uv[[i]][[j]]=solve(SIG_uv[[i]][[j]])
  }
  
}

theta_cmluv=numeric()
theta_eggeruv=numeric()
theta_lassouv=numeric()
theta_medianuv=numeric()
theta_ivwuv=numeric()
se_cmluv=numeric()
se_eggeruv=numeric()
se_lassouv=numeric()
se_medianuv=numeric()
se_ivwuv=numeric()



theta_cmluvDP=numeric()
theta_eggeruvDP=numeric()
theta_lassouvDP=numeric()
theta_medianuvDP=numeric()
theta_ivwuvDP=numeric()

for(i1 in 1:length(d)){
  dDP=d[i1]
  
  
  niv_uv=length(IVuv[[i1]])
  
  SIGinvDP_uv=SIGinv_uv[[i1]]
  SIGDP_uv=SIG_uv[[i1]]
  
  se_expDP_uv=se_exp_uv[[i1]]
  se_outDP_uv=se_out_uv[[i1]]
  beta_exp_d_uv=beta_exp_uv[[i1]]
  beta_out_d_uv=beta_out_uv[[i1]]
  for(i2 in 1:DP){
    
    beta_expDP_uv=matrix(0,nrow = niv_uv,ncol =1)
    beta_outDP_uv=matrix(0,nrow = niv_uv,ncol =1)
    
    for(i in 1:niv_uv){
      beta_uv=mvrnorm(n=1,mu=c(beta_exp_d_uv[i],beta_out_d_uv[i]),Sigma = SIGDP_uv[[i]])
      beta_expDP_uv[i,]=beta_uv[1:(length(beta_uv)-1)]
      beta_outDP_uv[i,]=beta_uv[length(beta_uv)]
      if(beta_expDP_uv[i,1]<0){
        beta_expDP_uv[i,]=-beta_expDP_uv[i,]
        beta_outDP_uv[i,]=-beta_outDP_uv[i,]
      }
    }
    
    obj_uv=mr_input(bx=beta_expDP_uv[,1],bxse = se_expDP_uv,by=beta_outDP_uv[,1],byse = se_outDP_uv)
    theta_cmluvDP[i2]=MVmr_cML(b_exp =beta_expDP_uv,b_out=matrix(beta_outDP_uv,ncol=1,nrow=length(beta_outDP_uv)),se_bx = matrix(se_expDP_uv,ncol = 1,nrow=length(se_expDP_uv)),Sig_inv_l = SIGinvDP_uv,
                               n=100000,K_vec = c(0:(niv_uv-2)),random_start = 1,min_theta_range = 0.1,
                               max_theta_range = 0.1,maxit = 500,thres = 0.01)$BIC_theta[1]
    theta_eggeruvDP[i2]=mr_egger(object = obj_uv)$Estimate[1]
    theta_lassouvDP[i2]=mvmr_lasso(bx=beta_expDP_uv,by=beta_outDP_uv[,1],seby=se_outDP_uv)$th_post[1]
    theta_ivwuvDP[i2]=mr_ivw(object = obj_uv)$Estimate[1]
    theta_medianuvDP[i2]=mr_median(object = obj_uv,iterations = 200)$Estimate[1]
    print(i2)
    
  }
  
  theta_cmluv[i1]=mean(theta_cmluvDP)
  theta_eggeruv[i1]=mean(theta_eggeruvDP)
  theta_lassouv[i1]=mean(theta_lassouvDP)
  theta_ivwuv[i1]=mean(theta_ivwuvDP)
  theta_medianuv[i1]=mean(theta_medianuvDP)
  se_cmluv[i1]=sd(theta_cmluvDP)
  se_eggeruv[i1]=sd(theta_eggeruvDP)
  se_lassouv[i1]=sd(theta_lassouvDP)
  se_ivwuv[i1]=sd(theta_ivwuvDP)
  se_medianuv[i1]=sd(theta_medianuvDP)
}

output=list(data.frame(theta_cmluv,theta_eggeruv,theta_lassouv,theta_ivwuv,
                       theta_medianuv,se_cmluv,se_eggeruv,se_lassouv,se_ivwuv,
                       se_medianuv))
save(output,file="BMI_DBPcorrected.RData")