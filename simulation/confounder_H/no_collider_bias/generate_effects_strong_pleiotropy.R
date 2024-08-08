setwd("D:/paper2/simulation/strong_pleiotropy/")
N_non_null=40
N_H_only=40
N_null=3000
invalid1=0.3
invalid2=0.5
invalid_b=0.3
Nsnps=N_non_null+N_null+N_H_only

 
dtotal=c(1,2,4)
betaUXtotal=numeric()
for(i in 1:length(dtotal)){
  betaUXtotal[i]=1/(dtotal[i]+1)
}
betaUY=1
betaUHtotal=list()
for(i in 1:length(dtotal)){
  betaUHtotal[[i]]=rep(1/(dtotal[i]+1),times=dtotal[i])
}
betaGU_fixed_balanced=c(rep(0,times=N_non_null),rep(0,times=N_H_only),rep(0,times=N_null))
betaGU_fixed_directional=c(runif(n=N_non_null,min=0.2,max=0.8),rep(0,times=N_H_only),rep(0,times=N_null))
betaGH_total=rbind(matrix(runif(n=N_non_null*max(dtotal),min=0,max=0.22),nrow=N_non_null,ncol=max(dtotal)),
                   matrix(runif(n=N_H_only*max(dtotal),min=0,max=0.22),nrow=N_H_only,ncol=max(dtotal)),
                   matrix(0,nrow =N_null,ncol = max(dtotal) ))
                   
betaXH_total=rnorm(n=max(dtotal),mean=0.5,sd=0) 
 
betaGX=c(runif(n=N_non_null,min=0,max=0.22),runif(n=N_H_only*invalid_b,min=0,max=0.22),rep(0,times=(1-invalid_b)*N_H_only),rep(0,times=N_null))
betaGX_valid=c(runif(n=N_non_null,min=0,max=0.22),runif(n=N_H_only*invalid_b,min=0,max=0),rep(0,times=(1-invalid_b)*N_H_only),rep(0,times=N_null))
betaHX_total=rnorm(n=max(dtotal),mean=0.5,sd=0)
betaGY_fixed_balanced=c(rnorm(n=N_non_null,mean=0,sd=0.2),rep(0,times=N_H_only),rep(0,times=N_null))
betaGY_fixed_directional=c(rnorm(n=N_non_null,mean=0.5,sd=0.25),rep(0,times=N_H_only),rep(0,times=N_null))
betaHY_total=c(0.2,0.1,0.3,0.4)
theta=1


#When H is a confounder
 
# 30% invalid
 
 
# Directional and correlated pleiotropy
for(i in 1:length(dtotal)){
  effects_confounder_H=list(
    betaGU=c(rep(0,times=(1-invalid1)*N_non_null),rep(1,invalid1*N_non_null),rep(1,N_H_only+N_null))*betaGU_fixed_directional,
    betaGH=betaGH_total[,1:dtotal[i]],
    betaGX=betaGX,
    betaHX=betaHX_total[1:dtotal[i]],
    betaGY=c(rep(0,times=(1-invalid1)*N_non_null),rep(1,invalid1*N_non_null),rep(1,N_H_only+N_null))*betaGY_fixed_directional,
    betaHY=betaHY_total[1:dtotal[i]],
    betaUX=betaUXtotal[i],
    betaUY=betaUY,
    betaUH=betaUHtotal[[i]],
    theta=theta
  )
  save(effects_confounder_H,file=paste("confounder_Hdim=",dtotal[i],"_directional_correlated_pleiotropy30.RData",sep=""))
}

# 50% invalid
 
# Directional and correlated pleiotropy
for(i in 1:length(dtotal)){
  effects_confounder_H=list(
    betaGU=c(rep(0,times=(1-invalid2)*N_non_null),rep(1,invalid2*N_non_null),rep(1,N_H_only+N_null))*betaGU_fixed_directional,
    betaGH=betaGH_total[,1:dtotal[i]],
    betaGX=betaGX,
    betaHX=betaHX_total[1:dtotal[i]],
    betaGY=c(rep(0,times=(1-invalid2)*N_non_null),rep(1,invalid2*N_non_null),rep(1,N_H_only+N_null))*betaGY_fixed_directional,
    betaHY=betaHY_total[1:dtotal[i]],
    betaUX=betaUXtotal[i],
    betaUY=betaUY,
    betaUH=betaUHtotal[[i]],
    theta=theta
  )
  save(effects_confounder_H,file=paste("confounder_Hdim=",dtotal[i],"_directional_correlated_pleiotropy50.RData",sep=""))
}


 