
 
setwd("/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/")


library(data.table)
DBP_use=as.data.frame(fread("GWAS_DBP.txt",sep="\t",header = TRUE))
BMI_use=as.data.frame(fread("GWAS_BMI.txt",sep="\t",header=TRUE))
SBP_use=as.data.frame(fread("GWAS_SBP.txt",sep="\t",header=TRUE))
CAD_use=as.data.frame(fread("GWAS_CAD.txt",sep="\t",header=TRUE))

 

DBP_use=DBP_use[!duplicated(DBP_use[,c("SNP")]),]
CAD_use=CAD_use[!duplicated(CAD_use[,c("SNP")]),]
SBP_use=SBP_use[!duplicated(SBP_use[,c("SNP")]),]
DBP_use=DBP_use[!duplicated(DBP_use[,c("SNP")]),]
 
BMI_use=BMI_use[order(BMI_use$SNP),]
CAD_use=CAD_use[order(CAD_use$SNP),] 
SBP_use=SBP_use[order(SBP_use$SNP),] 
DBP_use=DBP_use[order(DBP_use$SNP),] 

nrow(na.omit(BMI_use))
nrow(na.omit(CAD_use))
nrow(na.omit(SBP_use))
nrow(na.omit(DBP_use))


BMI_ALT=BMI_use$ALT
CAD_ALT=CAD_use$ALT
SBP_ALT=SBP_use$ALT
DBP_ALT=DBP_use$ALT
BMI_beta=BMI_use$beta
CAD_beta=CAD_use$beta
SBP_beta=SBP_use$beta
DBP_beta=DBP_use$beta
BMI_z=BMI_use$z
CAD_z=CAD_use$z
SBP_z=SBP_use$z
DBP_z=DBP_use$z

for(i in 1:length(CAD_ALT)){
  if(BMI_ALT[i]!=CAD_ALT[i]){
    if((BMI_ALT[i]%in%c("A","T")&CAD_ALT[i]%in%c("C","G"))|(BMI_ALT[i]%in%c("C","G")&CAD_ALT[i]%in%c("A","T"))){
      BMI_beta[i]=-BMI_beta[i]
      BMI_z[i]=-BMI_z[i]
      BMI_ALT[i]=CAD_ALT[i] 
    }
     
  }
  if(DBP_ALT[i]!=CAD_ALT[i]){
    if((SBP_ALT[i]%in%c("A","T")&CAD_ALT[i]%in%c("C","G"))|(SBP_ALT[i]%in%c("C","G")&CAD_ALT[i]%in%c("A","T"))){
      SBP_beta[i]=-SBP_beta[i]
      SBP_z[i]=-SBP_z[i]
      SBP_ALT[i]=CAD_ALT[i] 
    }
  }
  if(SBP_ALT[i]!=CAD_ALT[i]){
    if((DBP_ALT[i]%in%c("A","T")&CAD_ALT[i]%in%c("C","G"))|(DBP_ALT[i]%in%c("C","G")&CAD_ALT[i]%in%c("A","T"))){
      DBP_beta[i]=-DBP_beta[i]
      DBP_z[i]=-DBP_z[i]
      DBP_ALT[i]=CAD_ALT[i] 
    }
  }
  print(i)
}

BMI_use$beta=BMI_beta
BMI_use$ALT=BMI_ALT
BMI_use$z=BMI_z
SBP_use$beta=SBP_beta
SBP_use$ALT=SBP_ALT
SBP_use$z=SBP_z 
DBP_use$beta=DBP_beta
DBP_use$ALT=DBP_ALT
DBP_use$z=DBP_z

BMI_use=BMI_use[order(BMI_use$POS),]
CAD_use=CAD_use[order(CAD_use$POS),]
DBP_use=DBP_use[order(DBP_use$POS),]
SBP_use=SBP_use[order(SBP_use$POS),]


write.table(BMI_use,file="GWAS_BMI_harmo.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(CAD_use,file="GWAS_CAD_harmo.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE) 
write.table(SBP_use,file="GWAS_SBP_harmo.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE) 
write.table(DBP_use,file="GWAS_DBP_harmo.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


