setwd("D:/paper2/real_data/BMI_BP")
d=c(1,2,5)
DP=100
library(MendelianRandomization)
library(MASS)
cov_name=list()
exp_name=numeric()

for(i in 1:length(d)){
  cov_name[[i]]=numeric()
  exp_name[i]=paste("GWAS_",d[i],"MetPCadj.BMI.glm.linear",sep="")
}
for(i in 1:length(d)){
  for(j in 1:d[i]){
    cov_name[[i]][j]=paste("GWAS.metPC",j,".glm.linear",sep="")
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
for(i in 1:length(d)){
  IV[[i]]=cov[[i]][[1]]$snp[cov[[i]][[1]]$snp%in%indsnp&cov[[i]][[1]]$p<5e-10]
  for(j in 1:d[i]){
    IV[[i]]=c(IV[[i]],cov[[i]][[j]]$snp[cov[[i]][[j]]$snp%in%indsnp&cov[[i]][[j]]$p<5e-10])
  }
  IV[[i]]=unique(IV[[i]])
}

beta_cov=list()
se_cov=list()
beta_exp=list()
se_exp=list()
for(i in 1:length(d)){
  beta_cov[[i]]=matrix(0,nrow=length(IV[[i]]),ncol=d[i])
  se_cov[[i]]=matrix(0,nrow=length(IV[[i]]),ncol=d[i])
  for(j in 1:d[i]){
    beta_cov[[i]][,j]=cov[[i]][[j]]$beta[cov[[i]][[j]]$snp%in%IV[[i]]]
    se_cov[[i]][,j]=cov[[i]][[j]]$se[cov[[i]][[j]]$snp%in%IV[[i]]]
  }
  beta_exp[[i]]=matrix(exp[[i]]$beta[exp[[i]]$snp%in%IV[[i]]],ncol=1)
  se_exp[[i]]=matrix(exp[[i]]$se[exp[[i]]$snp%in%IV[[i]]],ncol=1)
}

for(i in 1:length(d)){
  null_z=cov[[i]][[1]]$p>0.1&exp[[i]]$p>0.1&exp[[i]]$snp%in%indsnp
  for(j in 1:d[i]){
    null_z=cov[[i]][[j]]$p>0.1&cov[[i]][[j]]$snp%in%indsnp&null_z
  }
}
nullz=list()
for(i in 1:length(d)){
   
    nullz[[i]]=cbind(cov[[i]][[1]]$z[null_z])
 if(i>1){
   for(j in 2:d[i]){
    nullz[[i]]=cbind(nullz[[i]],cov[[i]][[j]]$z[null_z])
 }
  }
  
  nullz[[i]]=cbind(nullz[[i]],exp[[i]]$z[null_z])
}
correlation=list()
for(i in 1:length(d)){
  correlation[[i]]=cor(nullz[[i]])
}
SIG=list()
SIGinv=list()
for(i in 1:length(d)){
  SIG[[i]]=list()
  SIGinv[[i]]=list()
}
for(i in 1:length(d)){
  for(j in 1:length(IV[[i]])){
    SIG[[i]][[j]]=diag(c(se_cov[[i]][j,],se_exp[[i]][j]),nrow=d[i]+1,ncol=d[i]+1)%*%correlation[[i]]%*%diag(c(se_cov[[i]][j,],se_exp[[i]][j]),nrow=d[i]+1,ncol=d[i]+1)
    SIGinv[[i]][[j]]=solve(SIG[[i]][[j]])
    }
}


covbtotal=list()
btotal=list()
for(i1 in 1:length(d)){
  dDP=d[i1]
  b_DP=matrix(0,nrow = DP,ncol=dDP)
  beta_cov_d=beta_cov[[i1]]
  se_cov_d=se_cov[[i1]]
  beta_exp_d=beta_exp[[i1]]
  se_exp_d=se_exp[[i1]]
  SIGDP=SIG[[i1]]
  SIGinvDP=SIGinv[[i1]]
  niv=length(IV[[i1]])
  for(i2 in 1:DP){
    betaXDP=matrix(0,nrow=niv,ncol=dDP)
    betaYDP=matrix(0,nrow = niv,ncol=1)
    for(i in 1:niv){
      beta=mvrnorm(n=1,mu=c(beta_cov_d[i,],beta_exp_d[i]),Sigma = SIGDP[[i]])
      betaXDP[i,]=beta[1:(length(beta)-1)]
      betaYDP[i,]=beta[length(beta)]
    }
    obj=MVmr_cML(b_exp =betaXDP,b_out=betaYDP,se_bx = se_cov_d,
                 Sig_inv_l = SIGinvDP ,n=100000,K_vec = c(0:(niv-dDP-2)),random_start =  1,
                 min_theta_range =0,max_theta_range = 0,maxit = 500,thres = 0.01)
    b_DP[i2,]=obj$BIC_theta
  }
  btotal[[i1]]=colMeans(b_DP)
  covbtotal[[i1]]=cov(b_DP)
}
exp_corrected=list()
beta_covtotal=matrix(0,nrow=nrow(cov[[1]][[1]]),ncol = max(d))
se_covtotal=matrix(0,nrow=nrow(cov[[1]][[1]]),ncol = max(d))
for(i in 1:max(d)){
  beta_covtotal[,i]=cov[[length(d)]][[i]]$beta
  se_covtotal[,i]=cov[[length(d)]][[i]]$se
}
for(i in 1:length(d)){
  beta=numeric()
  se=numeric()
  se_old=exp[[i]]$se
  beta_old=exp[[i]]$beta
  beta_met=matrix(beta_covtotal[,1:d[i]],nrow=nrow(beta_covtotal))
  se_met=matrix(se_covtotal[,1:d[i]],nrow = nrow(se_covtotal))
  b=btotal[[i]]
  covb=covbtotal[[i]]
  for(j in 1:nrow(exp[[i]])){
    beta[j]=beta_old[j]-t(b)%*%beta_met[j,]
    se[j]=se_old[j]+sum(covb*diag(se_met[j,],nrow=d[i],ncol=d[i]))+sum(b^2*se_met[j,]^2)+t(beta_met[j,])%*%covb%*%beta_met[j,]
    print(j)
  }
  z=beta/se
  p=(1-pnorm(abs(z)))*2
  exp_corrected[[i]]=data.frame(
    chr=exp[[i]]$chr,
    pos=exp[[i]]$pos,
    snp=exp[[i]]$snp,
    ref=exp[[i]]$ref,
    alt=exp[[i]]$alt,
    a1=exp[[i]]$a1,
    n=exp[[i]]$n,
    test=exp[[i]]$test,
    beta=beta,
    se=se,
    z=z,
    p=p
  )
}
filenames=numeric()
for(i in 1:length(d)){
  filenames[i]=paste("GWAS_BMI",d[i],"MetPC_corrected.txt",sep="")
}
for(i in 1:length(d)){
  write.table(exp_corrected[[i]],sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE,file = filenames[i])
}