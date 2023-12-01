# Graph showing the percentage of variants (SNPs) mapped to RNYs (± 50 kb) 
# in 1,000 random sets of 3,847 SNPs (no filter criteria) and the observed 
# percentage in the GWAS catalog of cancer risk variants (3.7%; 144/3,845)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(biomaRt) # ‘2.56.1’
library(stringr) # ‘1.5.0’
library(ggplot2) # ‘3.4.3’

## GWAS cancer related SNPs mapped to YRNA
# GWAS cancer related SNP (GRCh38)
snp_list <- readRDS("Data/GWAS_cancer_related_snps.rds")
snp_list.clean <- snp_list %>% filter(CHR_ID %in% seq(1,22)) %>% distinct(SNP,CHR_ID,CHR_POS) %>% transmute(snpID=SNP,chr=CHR_ID,pos=as.integer(CHR_POS))

# transform from hg38 to hg19
# glimpse(snp_list.clean)
ensemblGRCh37.SNPs <- useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snpGRCh37.ensembl <- getBM(attributes = c("refsnp_id","chr_name","chrom_start"), 
                           filters = c("snp_filter","chr_name"), values = list(snp_list.clean %>% pull(snpID),snp_list.clean %>% pull(chr)), 
                           mart=ensemblGRCh37.SNPs)
snp_list.GRCh37 <- snpGRCh37.ensembl %>% rename(snpID=refsnp_id,chr=chr_name,pos=chrom_start)
totalSNP <- dim(snp_list.GRCh37)[1]

# Annotated YRNA BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                                 attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                                "start_position","end_position","gene_biotype","strand"), 
                                 mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))

# Map SNP to YRNA window (+/- 50kb)
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(snp_list.GRCh37 ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) %>%
  left_join(snp_list,by=c("snpID"="SNP"),relationship="many-to-many") %>% 
  transmute(`Cancer Study`=`DISEASE/TRAIT`,`SNP ID`=snpID,Chromosome=chromosome_name,`Position (bp)`=pos,
            `Ensembl ID`=ensembl_gene_id,`Gene start (bp)`=start_position,`Gene end (bp)`=end_position,
            `Pubmed ID GWAS`=PUBMEDID,`First author`=`FIRST AUTHOR`,Journal=JOURNAL,`Genetic context`=CONTEXT,
            `Risk alelle frequency`=`RISK ALLELE FREQUENCY`,`Cancer risk p value`=`P-VALUE`,`OR or beta effect`=`OR or BETA`)

# Number SNP mapped to YRNA window (+/- 50kb)
yrna.snp.gwas <- length(unique(biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(snp_list.GRCh37 ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) %>% pull(snpID)) )

# Test number random SNPs mapped to RNY window 
# Read reference list of SNP with requirements r2<0.8 and SNV (from referenceSNPrandom.R)
load(file="RData/SNP_EUR_proc.RData")

# Generate 1000 random SNPs list (totalSNP)
SNP.random <- list()
for(r in seq(1,20)) {
  print(r)
  SNPs <- sample(SNP.EUR.SNV %>% pull(snpID),totalSNP)
  
  SNPs.split <- data.frame(snpID=SNPs,as.data.frame(do.call(rbind,str_split(SNPs, ":"))))
  colnames(SNPs.split) <- c("snpID","chr","pos")
  
  SNPs.EUR.sel <- SNP.EUR.SNV %>% filter(snpID %in% SNPs) %>% left_join(SNPs.split,by=c("snpID"))
  SNPs.EUR.sel <- SNPs.EUR.sel %>% mutate(pos=as.integer(pos),chr=as.integer(chr))
  SNP.random[[paste0("random",r)]] <- SNPs.EUR.sel
}

# Number SNPs mapped to YRNA window
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
pSNPYRNApleio <- yrna.snp.gwas/totalSNP
pSNPYRNArandom <- nSNPYRNArandom/totalSNP

mean(pSNPYRNArandom >= pSNPYRNApleio)

# Plot
postscript(file="Output/Figure5e.ps")
data.frame(n=pSNPYRNArandom) %>% ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white") +
  theme_classic() + 
  xlim(0.015,0.04) +
  ylab("Number of random SNP sets") +
  xlab("Proportion of SNP with YRNA") +
  geom_segment(aes(x=pSNPYRNApleio,y=100,xend=pSNPYRNApleio,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off()  
