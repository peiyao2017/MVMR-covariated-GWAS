 library(data.table)
 library(TwoSampleMR)
 library(ieugwasr)
bmi=read.table("D:/paper2/real_data/SBP_CAD/GWAS_BMI_harmo.txt",header = TRUE,sep="\t")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
bmi=data.table(bmi)
bmi=bmi[order(bmi$CHR),]
data_clump=data.frame(chr_name=bmi$CHR,chrom_start=bmi$POS,SNP=bmi$SNP) 
indSNPs=clump_data(dat=data_clump,pop="EUR",clump_r2 = 0.001,clump_kb = 300)
ind_snps=indSNPs$SNP
ind_snps=list(ind_snps)
save(ind_snps,file="D:/paper2/real_data/SBP_CAD/ind_snps.RData")