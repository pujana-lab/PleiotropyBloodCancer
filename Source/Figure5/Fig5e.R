# README Test random SNPs GWAS cancer (Roderic) + SNV (REF 1b and ALT 1b). Proportions of SNP within RNY window ####
rm(list=ls())
gc()
pathR <- c("/media/mapardo/SeagateBasic/procure/RData/pleiotStudy")
pathExt <- c("/media/mapardo/SeagateBasic/procure/dataExt")

box::use(dplyr[...])
box::use(tidyr[...])

# 
load(file=file.path(pathR,"biomart37YRNA_RNY816.RData"))
# GWAS cancer related SNP (GRCh38)
snp_list <- readRDS("data/GWAS_cancer_related_snps.rds")
snp_list.clean <- snp_list %>% filter(CHR_ID %in% seq(1,22)) %>% distinct(SNP,CHR_ID,CHR_POS) %>% transmute(snpID=SNP,chr=CHR_ID,pos=as.integer(CHR_POS))
# 3847

# transform from hg38 to hg19
glimpse(snp_list.clean)
ensemblGRCh37.SNPs <- biomaRt::useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snpGRCh37.ensembl <- biomaRt::getBM(attributes = c("refsnp_id","chr_name","chrom_start"), 
                                    filters = c("snp_filter","chr_name"), values = list(snp_list.clean %>% pull(snpID),snp_list.clean %>% pull(chr)), 
                                    mart=ensemblGRCh37.SNPs)
snp_list.GRCh37 <- snpGRCh37.ensembl %>% rename(snpID=refsnp_id,chr=chr_name,pos=chrom_start)
save(snp_list.GRCh37,file=file.path(pathR,"SNPGWAScancerGRCh37.RData"))
# 3845
glimpse(snp_list.GRCh37)
glimpse(snp_list)
glimpse(biomart37.yrnarny)
#
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(snp_list.GRCh37 ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) %>%
  left_join(snp_list,by=c("snpID"="SNP"),relationship="many-to-many") %>% 
  transmute(`Cancer Study`=`DISEASE/TRAIT`,`SNP ID`=snpID,Chromosome=chromosome_name,`Position (bp)`=pos,
            `Ensembl ID`=ensembl_gene_id,`Gene start (bp)`=start_position,`Gene end (bp)`=end_position,
            `Pubmed ID GWAS`=PUBMEDID,`First author`=`FIRST AUTHOR`,Journal=JOURNAL,`Genetic context`=CONTEXT,
            `Risk alelle frequency`=`RISK ALLELE FREQUENCY`,`Cancer risk p value`=`P-VALUE`,`OR or beta effect`=`OR or BETA`)
glimpse(yrna.snp)
write.table(yrna.snp,file=file.path("output/18.testRandomYRNA_GWAS/YRNA_mapped_SNPs_GWAScancer.csv"),sep=",",row.names = FALSE)

# 
yrna.snp.gwas <- countSNPYRNAGRCh37(snp_list.GRCh37)
save(yrna.snp.gwas,file=file.path(pathR,"YRNASNPGWAS.RData"))

countSNPYRNAGRCh37 <- function(df) {
  yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
    left_join(df ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
    filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
    filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
  n.snp <- length(unique(yrna.snp %>% pull(snpID)))
  return(n.snp)
}

# random SNP sets
load(file=file.path(pathR,"SNP_EUR_proc.RData"))

SNP.random <- list()
for(r in seq(1,1000)) {
  print(r)
  SNPs <- sample(SNP.EUR.SNV %>% pull(snpID),3845)
  
  SNPs.split <- data.frame(snpID=SNPs,as.data.frame(do.call(rbind,stringr::str_split(SNPs, ":"))))
  colnames(SNPs.split) <- c("snpID","chr","pos")
  
  SNPs.EUR.sel <- SNP.EUR.SNV %>% filter(snpID %in% SNPs) %>% left_join(SNPs.split,by=c("snpID"))
  SNPs.EUR.sel <- SNPs.EUR.sel %>% mutate(pos=as.integer(pos),chr=as.integer(chr))
  SNP.random[[paste0("random",r)]] <- SNPs.EUR.sel
}
save(SNP.random,file=file.path(pathR,"SNPrandomGWASSNP.RData"))

# test number SNPs within RNY window 
load(file=file.path(pathR,"SNPrandomGWASSNP.RData"))
load(file=file.path(pathR,"biomart37YRNA_RNY816.RData"))
load(file=file.path(pathR,"YRNASNPGWAS.RData"))

#load(file=file.path(pathdata,"R","YRNASNPpleio816.RData"))

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
nSNPYRNArandom
pSNPYRNApleio <- yrna.snp.gwas/3845
pSNPYRNArandom <- nSNPYRNArandom/3845

mean(pSNPYRNArandom >= pSNPYRNApleio)

postscript(file="output/18.testRandomYRNA_GWAS/YRNArandomTestGWASSNP.ps")
data.frame(n=pSNPYRNArandom) %>% ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white") +
  theme_classic() + 
  xlim(0.015,0.04) +
  ylab("Number of random SNP sets") +
  xlab("Proportion of SNP with YRNA") +
  geom_segment(aes(x=pSNPYRNApleio,y=100,xend=pSNPYRNApleio,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off()  



## OLD
# relations SNP vs YRNA GRCh38
load(file=file.path(pathR,"biomart38YRNA_RNY771.RData"))
# GWAS cancer related SNP (GRCh38)
snp_list <- readRDS("data/GWAS_cancer_related_snps.rds")
snp_list.clean <- snp_list %>% filter(CHR_ID %in% seq(1,22)) %>% distinct(SNP,CHR_ID,CHR_POS) %>% transmute(snpID=SNP,chr=CHR_ID,pos=as.integer(CHR_POS))
# 3847
yrna.snp.gwas <- countSNPYRNAGRCh38(snp_list.clean)

countSNPYRNAGRCh38 <- function(df) {
  yrna.snp <- biomart38.yrnarny %>%
    left_join(df ,by=c("chromosome_name"="chr"),relationship="many-to-many") %>%
    filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>%
    filter(!(chromosome_name == 6 & between(pos,29500000,33500000)))
  return(yrna.snp)
}

save(yrna.snp.gwas,file=file.path(pathR,"YRNASNPGWAS.RData"))


# random SNP sets
load(file=file.path(pathR,"SNP_EUR_proc.RData"))

SNP.random <- list()
for(r in seq(1,2000)) {
  #for(r in seq(1,10)) {
  print(r)
  SNPs <- sample(SNP.EUR.SNV %>% pull(snpID),3847)
  
  SNPs.split <- data.frame(snpID=SNPs,as.data.frame(do.call(rbind,stringr::str_split(SNPs, ":"))))
  colnames(SNPs.split) <- c("snpID","chr","pos")
  
  SNPs.EUR.sel <- SNP.EUR.SNV %>% filter(snpID %in% SNPs) %>% left_join(SNPs.split,by=c("snpID"))
  SNPs.EUR.sel <- SNPs.EUR.sel %>% mutate(pos=as.integer(pos),chr=as.integer(chr))
  SNP.random[[paste0("random",r)]] <- SNPs.EUR.sel
}
save(SNP.random,file=file.path(pathR,"SNPrandomGWASSNP2000.RData"))

# test number SNPs within RNY window
load(file=file.path(pathR,"SNPrandomGWASSNP2000.RData"))
load(file=file.path(pathR,"biomart37YRNA_RNY816.RData"))
load(file=file.path(pathR,"YRNASNPGWAS.RData"))
#load(file=file.path(pathR,"YRNASNPpleio816.RData"))

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
nSNPYRNArandom
pSNPYRNApleio <- length(unique(yrna.snp.gwas %>% pull(snpID)))/3847
pSNPYRNArandom <- nSNPYRNArandom/3847

sum(pSNPYRNArandom >= pSNPYRNApleio)

postscript(file="output/18.testRandomYRNA_GWAS/YRNArandomTestGWASSNP2000.ps")
data.frame(n=pSNPYRNArandom) %>% ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white") +
  theme_classic() +
  xlim(0.015,0.04) +
  ylab("Number of random SNP sets") +
  xlab("Proportion of SNP with YRNA") +
  geom_segment(aes(x=pSNPYRNApleio,y=100,xend=pSNPYRNApleio,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off()

