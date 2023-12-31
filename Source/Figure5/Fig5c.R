# Graph showing the percentage of variants (SNPs) mapped to RNYs (± 50 kb) 
# in 1,000 random sets of 8,155 SNPs (European MAF > 0.01 and r 2 < 0.8) 
# and the observed percentage in the blood trait–cancer pleiotropy set (6.6%; 270/4,093)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(data.table) # ‘1.14.8’
library(stringr) # ‘1.5.0’
library(ggplot2) # ‘3.4.3’

## Pleitropic SNPs mapped to YRNA window
# read pleiotropic SNP list (pleio_loci.tsv). Filter to conjfdr < 0.05
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# number unique pleiotropic SNPs
nSNPlead <- dim(leadSNP005 %>% distinct(snpid))[1]
# YRNA + RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
# pleiotropic SNPs next to YRNA (+/- 50kb)
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
n.yrna <- length(unique(yrna.snp %>% pull(ensembl_gene_id)))
n.snp <- length(unique(yrna.snp %>% pull(snpid)))

## Test number random SNPs mapped to RNY window
# read reference list of SNP with requirements r2<0.8 and SNV (referenceSNPrandom.R)
load(file="RData/SNP_EUR_proc.RData")
# generate 1000 random SNPs list
SNP.random <- list()
for(r in seq(1,1000)) {
  print (paste0("Random",r))
  SNPs <- sample(SNP.EUR.SNV %>% pull(snpID),11000)
  SNPs.split <- data.frame(snpID=SNPs,as.data.frame(do.call(rbind,str_split(SNPs, ":"))))
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
# SNPs mapped to YRNA window for each random SNP set
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
# plot
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
