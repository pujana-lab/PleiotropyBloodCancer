# Reference list of SNPs to generate random set of SNPs (requirements: r2<0.8 and SNV)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(data.table) # ‘1.14.8’

# 1. List of SNPs r2<0.8 in chromosomes 1-22
# Download list of SNP with r2<0.8 
# https://data.broadinstitute.org/mpg/snpsnap/database_download.html 
SNP.EUR <- fread(file="Data/ld0.8_collection.tab")
# Filter only chromosomes 1-22
SNP.EUR.1_22 <- SNP.EUR %>% filter(!(grepl("23:",snpID)))
# Save
write.table(SNP.EUR.1_22 %>% pull(snpID),file="Data/SNP_EUR_1_22.vcf",quote=FALSE,row.names = FALSE)
# Run next command in bash line command
# sort SNP_EUR_1_22.vcf > SNP_EUR_1_22.sort.vcf

# 2. List of SNPs only SNV
# Download SNP info (REF, ALT and position)
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
# Filter SNP: only SNV
# Run next commands in bash line command
# gunzip -c All_20180423.vcf.gz | grep -v '#' | cut -f1,2,4,5 | awk '$3~/^[ACGT]$/' | awk '$4~/^[ACGT]$/' | awk {'print $1":"$2'} > All_20180423_onlycoord_SNV.vcf
# sort All_20180423_onlycoord_SNV.vcf > All_20180423_onlycoord_SNV.sort.vcf

# 3. SNPs only SNV, r2<0.8 and chromosomes 1-22
# comm -12 SNP_EUR_1_22.sort.vcf All_20180423_onlycoord_SNV.sort.vcf > SNP_EUR_1_22.SNV.vcf
SNV.EUR <- fread(file="Data/SNP_EUR_1_22.SNV.vcf",header=FALSE,col.names=c("SNP")) %>% pull(SNP)
SNP.EUR.SNV <- SNP.EUR.1_22 %>% filter(snpID %in% SNV.EUR)

save(SNP.EUR.SNV,file="RData/SNP_EUR_proc.RData")