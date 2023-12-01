# Graph showing the percentage of variants (SNPs) mapped to RNYs (± 50 kb) 
# in 1,000 random sets of 8,155 SNPs (European MAF > 0.01 and r 2 < 0.8) 
# and the observed percentage in the blood trait–cancer pleiotropy set (6.6%; 270/4,093)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(data.table) # ‘1.14.8’

# Number of distinct pleiotropic SNPs
load(file="RData/leadSNP005.RData")
nSNPlead <- dim(leadSNP005 %>% distinct(snpid))[1]

# Reference list of SNPs to generate random set of SNPs (requirements: r2<0.8 and SNV)
# 1. List of SNPs r2<0.8 in chromosomes 1-22
# Download list of SNP with r2<0.8 
# https://data.broadinstitute.org/mpg/snpsnap/database_download.html 
SNP.EUR <- fread(file=file.path(pathdata,"BroadInstitute","ld0.8_collection.tab"))
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

# Generate 1000 random list of nSNPlead SNPs
SNP.random <- list()
for(r in seq(1,1000)) {
  print (paste0("Random",r))
  SNPs <- sample(SNP.EUR.SNV %>% pull(snpID),11000)
  
  SNPs.split <- data.frame(snpID=SNPs,as.data.frame(do.call(rbind,stringr::str_split(SNPs, ":"))))
  colnames(SNPs.split) <- c("snpID","chr","pos")
  
  SNPs.EUR.sel <- SNP.EUR.SNV %>% filter(snpID %in% SNPs) %>% left_join(SNPs.split,by=c("snpID"))
  SNPs.EUR.sel <- SNPs.EUR.sel %>% mutate(pos=as.integer(pos),chr=as.integer(chr))
  
  SNPs.EUR.final <- NULL
  for (chrom in seq(1,22)) {
    print(chrom)
    SNPs.EUR.chr <- SNPs.EUR.sel %>% filter(chr==chrom)
    for (i in seq(1,nrow(SNPs.EUR.chr)-1)) {
      if (SNPs.EUR.chr$friends_ld08[i] == 0 | is.na(SNPs.EUR.chr$friends_ld08[i])) {
        SNPs.EUR.final <- rbind(SNPs.EUR.final,SNPs.EUR.chr[i,])
      } else {
        LD <- FALSE
        for (j in seq(i+1,nrow(SNPs.EUR.chr))) {
          if ((SNPs.EUR.chr$pos[j] >= SNPs.EUR.chr$loci_upstream[i]) & (SNPs.EUR.chr$pos[j] <= SNPs.EUR.chr$loci_downstream[i])) {
            LD <- TRUE
          }
        }
        if (!LD) {
          SNPs.EUR.final <- rbind(SNPs.EUR.final,SNPs.EUR.chr[i,])
        }
      }
    }
  }
  SNPs <- sample(SNPs.EUR.final %>% pull(snpID),nSNPlead)
  SNP.random[[paste0("random",r)]] <- SNPs.EUR.final %>% filter(snpID %in% SNPs)
}

# Test number SNPs within RNY window
# YRNA + RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))

# Pleiotropic SNPs next to YRNA (+/- 50kb)
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
n.yrna <- length(unique(yrna.snp %>% pull(ensembl_gene_id)))
n.snp <- length(unique(yrna.snp %>% pull(snpid)))

# SNPs next to YRNA for each random SNP set
nSNPYRNArandom <- lapply(SNP.random,countSNPYRNAGRCh37)
countSNPYRNAGRCh37 <- function(df) {
  yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>%
    left_join(df ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>%
    filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>%
    filter(!(chromosome_name == 6 & between(pos,29500000,33500000)))
  n.snp <- length(unique(yrna.snp %>% pull(snpID)))
  return(n.snp)
}

nSNPYRNArandom <- unlist(nSNPYRNArandom)
pSNPYRNApleio <- length(unique(yrna.snp %>% pull(snpid)))/nSNPlead
pSNPYRNArandom <- nSNPYRNArandom/nSNPlead

# Plot
postscript(file="Output/Figure5c.ps")
data.frame(n=pSNPYRNArandom) %>% ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white") +
  theme_classic() +
  xlim(0.015,0.07) +
  ylab("Number of random SNP sets") +
  xlab("Proportion of SNP with YRNA") +
  geom_segment(aes(x=pSNPYRNApleio,y=200,xend=pSNPYRNApleio,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off()
