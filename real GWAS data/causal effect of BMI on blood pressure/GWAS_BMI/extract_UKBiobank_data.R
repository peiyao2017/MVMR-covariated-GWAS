library(data.table)
library(dplyr)

#load metabolic data

x=fread(file="/panfs/jay/groups/20/panwei/shared/UKBiobankIndiv/metabolites_old/ukb49780.tab")
x %>% select(ends_with('.0.0'),f.eid) -> x
no_phe_id = which(rowSums(is.na(x))==251)
x = x[-no_phe_id,]
metIID=x$f.eid

x=cbind(metIID,x)
x=select(x,-c("f.eid"))
met=x
met=as.data.frame(met)
for(i in 2:ncol(met)) {
  met[,i][is.na(met[,i])]=mean(met[,i], na.rm = TRUE)
}

select1=c("f.eid","f.22027.0.0","f.22019.0.0","f.22001.0.0","f.22021.0.0","f.21000.0.0",
          "f.31.0.0","f.34.0.0",
          "f.22009.0.1","f.22009.0.2","f.22009.0.3","f.22009.0.4","f.22009.0.5",
          "f.22009.0.6","f.22009.0.7","f.22009.0.8","f.22009.0.9","f.22009.0.10",
          "f.21001.0.0","f.48.0.0","f.49.0.0","f.4079.0.0","f.4080.0.0")
data=fread("/panfs/jay/groups/20/panwei/shared/UKBiobankIndiv/ukb49020.tab",select=select1)
data%>%filter(is.na(f.22027.0.0)&is.na(f.22019.0.0)&f.22001.0.0==f.31.0.0&f.22021.0.0 %in% c(0,1)&f.21000.0.0 %in% c(1,1001,1002,1003)) %>%
  select(-f.22027.0.0,-f.22019.0.0,-f.22001.0.0,-f.22021.0.0,-f.21000.0.0) ->data_clean->y
FID=y$f.eid
IID=y$f.eid
BMI=y$f.21001.0.0
sex=y$f.31.0.0
age=y$f.34.0.0
GenPC1=y$f.22009.0.1
GenPC2=y$f.22009.0.2
GenPC3=y$f.22009.0.3
GenPC4=y$f.22009.0.4
GenPC5=y$f.22009.0.5
GenPC6=y$f.22009.0.6
GenPC7=y$f.22009.0.7
GenPC8=y$f.22009.0.8
GenPC9=y$f.22009.0.9
GenPC10=y$f.22009.0.10
BMI=data.frame(FID=FID,IID=IID,sex=sex,age=age,BMI=BMI,GenPC1=GenPC1,GenPC2=GenPC2,GenPC3=GenPC3,GenPC4=GenPC4,GenPC5=GenPC5,
               GenPC6=GenPC6,GenPC7=GenPC7,GenPC8=GenPC8,GenPC9=GenPC9,GenPC10=GenPC10)
for(i in 4:ncol(BMI)){
  BMI[,i][is.na(BMI[,i])]=mean(BMI[,i], na.rm = TRUE)
}

met.20PCnoGen=read.table(file="/panfs/jay/groups/20/panwei/wan01299/collider_bias/collider_bias/M2/BMI/MET_20PC_NoGene.txt",sep="\t",header=TRUE)
ID=intersect(metIID,BMI$IID)
ID=intersect(ID,met.20PCnoGen$IID)
BMI=BMI[BMI$IID%in%ID,]
met1=met[met$metIID%in%ID,]
met.20PCnoGen=met.20PCnoGen[met.20PCnoGen$IID%in%ID,]





met.X=met1[,-c(1)]
colnames(met.X)=NULL
met.pca = prcomp(met.X,scale=TRUE)
met.20PC = met.pca$x[,1:20]
namesMPC=numeric()
namesMPCnoGen=numeric()
for(i in 1:20){
  namesMPC[i]=paste("metPC",i,sep="")
  namesMPCnoGen[i]=paste("metPCnoGen",i,sep="")
}
colnames(met.20PC)=namesMPC
colnames(met.20PCnoGen)=c("FID","IID",namesMPCnoGen)
BMItotal=cbind(BMI,met.20PC,met.20PCnoGen[,-c(1,2)])

write.table(BMItotal,file="/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/gwas_BMI/BMI.txt",sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
write.table(data.frame(ID,ID),file="/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/gwas_BMI/BMI_ID.txt",sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)


select1=c("f.eid","f.22027.0.0","f.22019.0.0","f.22001.0.0","f.22021.0.0","f.21000.0.0",
          "f.31.0.0","f.34.0.0",
          "f.22009.0.1","f.22009.0.2","f.22009.0.3","f.22009.0.4","f.22009.0.5",
          "f.22009.0.6","f.22009.0.7","f.22009.0.8","f.22009.0.9","f.22009.0.10",
          "f.21001.0.0","f.48.0.0","f.49.0.0","f.4079.0.0","f.4080.0.0")
data=fread("/panfs/jay/groups/20/panwei/shared/UKBiobankIndiv/ukb49020.tab",select=select1)
data%>%filter(is.na(f.22027.0.0)&is.na(f.22019.0.0)&f.22001.0.0==f.31.0.0&f.22021.0.0 %in% c(0,1)&f.21000.0.0 %in% c(1,1001,1002,1003)) %>%
  select(-f.22027.0.0,-f.22019.0.0,-f.22001.0.0,-f.22021.0.0,-f.21000.0.0) ->data_clean->y
FID=y$f.eid
IID=y$f.eid

sex=y$f.31.0.0
age=y$f.34.0.0

SBP=y$f.4080.0.0
DBP=y$f.4079.0.0
GenPC1=y$f.22009.0.1
GenPC2=y$f.22009.0.2
GenPC3=y$f.22009.0.3
GenPC4=y$f.22009.0.4
GenPC5=y$f.22009.0.5
GenPC6=y$f.22009.0.6
GenPC7=y$f.22009.0.7
GenPC8=y$f.22009.0.8
GenPC9=y$f.22009.0.9
GenPC10=y$f.22009.0.10
BP=data.frame(FID=FID,IID=IID,sex=sex,age=age,SBP=SBP,DBP=DBP ,GenPC1=GenPC1,GenPC2=GenPC2,GenPC3=GenPC3,GenPC4=GenPC4,GenPC5=GenPC5,
              GenPC6=GenPC6,GenPC7=GenPC7,GenPC8=GenPC8,GenPC9=GenPC9,GenPC10=GenPC10)
for(i in 4:ncol(BP)){
  BP[,i][is.na(BP[,i])]=mean(BP[,i], na.rm = TRUE)
}
BMI_ID=read.table(file="/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/gwas_BMI/BMI_ID.txt",sep="\t",header=FALSE)
BP=BP[!BP$IID%in%BMI_ID[,1],]
BP_ID=BP$IID
write.table(BP,file="/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/gwas_BP/BP.txt",sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
write.table(data.frame(BP_ID,BP_ID),file="/panfs/jay/groups/20/panwei/wan01299/paper2/real_data/gwas_BP/BP_ID.txt",sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
