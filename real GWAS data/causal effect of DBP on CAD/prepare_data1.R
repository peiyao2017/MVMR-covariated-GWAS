
 

 
setwd("/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/")
library(data.table)
DBP_raw=as.data.frame(fread("ieu-b-39.vcf",sep="\t",header = TRUE))
BMI_raw=as.data.frame(fread("ieu-a-2.vcf",sep="\t",header=TRUE))
SBP_raw=as.data.frame(fread("ieu-b-38.vcf",sep="\t",header=TRUE))
CAD_raw=as.data.frame(fread("ebi-a-GCST005195.vcf",sep="\t",header=TRUE))

DBP_raw=na.omit(DBP_raw)
BMI_raw=na.omit(BMI_raw)
SBP_raw=na.omit(SBP_raw)
CAD_raw=na.omit(CAD_raw)

rs_BMI=BMI_raw$ID[BMI_raw$FILTER=="PASS"]
rs_DBP=DBP_raw$ID[DBP_raw$FILTER=="PASS"]
rs_SBP=SBP_raw$ID[SBP_raw$FILTER=="PASS"]
rs_CAD=CAD_raw$ID[CAD_raw$FILTER=="PASS"]


rs=intersect(rs_BMI,rs_DBP)  
rs=intersect(rs,rs_SBP)
rs=intersect(rs,rs_CAD)
DBP_use=DBP_raw[DBP_raw$ID%in%rs,]
BMI_use=BMI_raw[BMI_raw$ID%in%rs,]
SBP_use=SBP_raw[SBP_raw$ID%in%rs,]
CAD_use=CAD_raw[CAD_raw$ID%in%rs,]
rm(DBP_raw)
rm(BMI_raw)
rm(SBP_raw)
rm(CAD_raw)

sub_SBP=SBP_use$`ieu-b-38`
sub_DBP=DBP_use$`ieu-b-39`
sub_BMI=BMI_use$`ieu-a-2`
sub_CAD=CAD_use$`EBI-a-GCST005195`

beta_CAD=numeric()
se_CAD=numeric()
beta_SBP=numeric()
se_SBP=numeric()
beta_DBP=numeric()
se_DBP=numeric()
beta_BMI=numeric()
se_BMI=numeric()

for(i in 1:length(sub_DBP)){
  beta_DBP[i]=strsplit(sub_DBP[i],split = ":")[[1]][1]
  se_DBP[i]=strsplit(sub_DBP[i],split = ":")[[1]][2]
  print(i)
}

for(i in 1:length(sub_CAD)){
  beta_CAD[i]=strsplit(sub_CAD[i],split = ":")[[1]][1]
  se_CAD[i]=strsplit(sub_CAD[i],split = ":")[[1]][2]
  print(i)
}
attention=which(is.na(as.numeric(se_CAD)))
for(i in attention){
  beta_CAD[i]=strsplit(sub_CAD[i],split = ":")[[1]][1]
  se_CAD[i]=strsplit(sub_CAD[i],split = ":")[[1]][4]
  print(i)
}
for(i in 1:length(sub_SBP)){
  beta_SBP[i]=strsplit(sub_SBP[i],split = ":")[[1]][1]
  se_SBP[i]=strsplit(sub_SBP[i],split = ":")[[1]][2]
  print(i)
}

for(i in 1:length(sub_BMI)){
  beta_BMI[i]=strsplit(sub_BMI[i],split = ":")[[1]][1]
  se_BMI[i]=strsplit(sub_BMI[i],split = ":")[[1]][2]
  print(i)
}
BMI_allele=BMI_use$ALT
DBP_allele=DBP_use$ALT
SBP_allele=SBP_use$ALT
CAD_allele=CAD_use$ALT

BMI_other=BMI_use$REF
DBP_other=DBP_use$REF
SBP_other=SBP_use$REF
CAD_other=CAD_use$REF


beta_BMI=as.numeric(beta_BMI)
se_BMI=as.numeric(se_BMI)
beta_DBP=as.numeric(beta_DBP)
se_DBP=as.numeric(se_DBP)
beta_CAD=as.numeric(beta_CAD)
se_CAD=as.numeric(se_CAD)
beta_SBP=as.numeric(beta_SBP)
se_SBP=as.numeric(se_SBP)


 

 
gwas_BMI=data.frame(CHR=BMI_use$`#CHROM`,POS=BMI_use$POS,SNP=BMI_use$ID, beta=beta_BMI, se=se_BMI, z=beta_BMI/se_BMI,p=1-pnorm(abs(beta_BMI/se_BMI)),REF=BMI_other,ALT=BMI_allele)
gwas_CAD=data.frame(CHR=CAD_use$`#CHROM`,POS=CAD_use$POS,SNP=CAD_use$ID, beta=beta_CAD, se=se_CAD, z=beta_CAD/se_CAD,p=1-pnorm(abs(beta_CAD/se_CAD)),REF=CAD_other,ALT=CAD_allele)
gwas_SBP=data.frame(CHR=SBP_use$`#CHROM`,POS=SBP_use$POS,SNP=SBP_use$ID, beta=beta_SBP, se=se_SBP, z=beta_SBP/se_SBP,p=1-pnorm(abs(beta_SBP/se_SBP)),REF=SBP_other,ALT=SBP_allele)
gwas_DBP=data.frame(CHR=DBP_use$`#CHROM`,POS=DBP_use$POS,SNP=DBP_use$ID, beta=beta_DBP, se=se_DBP, z=beta_DBP/se_DBP,p=1-pnorm(abs(beta_DBP/se_DBP)),REF=DBP_other,ALT=DBP_allele)
write.table(gwas_BMI,file="GWAS_BMI.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(gwas_CAD,file="GWAS_CAD.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE) 
write.table(gwas_SBP,file="GWAS_SBP.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE) 
write.table(gwas_DBP,file="GWAS_DBP.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)