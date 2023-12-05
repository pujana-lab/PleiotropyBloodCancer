# Graphs showing the number of variants (SNPs) identified as pleiotropic in RNYs (± 50 kb) and 
# correlated (European r 2 > 0.4, left panel; r 2 > 0.8, right panel) with SLE GWAS catalog variants, 
# and compared with the results of equivalent 1,000 random variant sets (European MAF > 0.01)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(readxl) # ‘1.4.3’
library(biomaRt) # ‘2.56.1’
library(LDlinkR) # ‘1.3.0’
# PLINK v1.90b6.26 64-bit (2 Apr 2022)  www.cog-genomics.org/plink/1.9/

### 1. Check how many Sjorgen and SLE SNPs have r2 >0.4 or r2 > 0.8 with pleiotropic SNPs ####
## Pleiotropic SNPs nearby YRNA (270)
# RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
# read pleiotropic SNPs
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  dplyr::rename(cancer_trait=trait1,
                blood_trait=trait2,
                pval_neoplasm=pval_trait1,
                pval_immunological=pval_trait2,
                zscore_neoplasm=zscore_trait1,
                zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# assignment of RNY to nearby SNP
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% dplyr::rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
snp.yrna.pleio <- yrna.snp %>% distinct(snpid,chromosome_name,pos) %>% transmute(snpid,chr=chromosome_name,pos)
write.table(snp.yrna.pleio %>% pull(snpid),"Output/SNPpleio.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

## Read Sjorgen SNPs (48)
rsid.sjogren <- read_excel(path="Data/Fig7b/GWAS-Catalog_Sjögren-Lupus_112922.xlsx", 
                                   sheet="Sjögren",col_names=TRUE) %>% dplyr::select(CHR_ID,SNPS)
# https://www.ebi.ac.uk/gwas/
ensemblGRCh37.SNPs <- useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snpGRCh37.sjorgren <- getBM(attributes = c("refsnp_id","chr_name","chrom_start"), 
                                     filters = c("snp_filter","chr_name"), values = list(rsid.sjogren %>% pull(SNPS),rsid.sjogren %>% pull(CHR_ID)), 
                                     mart=ensemblGRCh37.SNPs)
snp.sjorgen <- snpGRCh37.sjorgren %>% transmute(snpid=refsnp_id,chr=chr_name,pos=chrom_start)
## Read SLE SNPs + clean SNP info (917)
rsid.lupus <- read_excel(path="Data/Fig7b/GWAS-Catalog_Sjögren-Lupus_112922.xlsx", 
                                 sheet="systemic lupus erythematosus",col_names=TRUE) %>% dplyr::select(CHR_ID,SNPS)
# https://www.ebi.ac.uk/gwas/
# get position for SNP with CHR_ID info
ensemblGRCh37.SNPs <- useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snpGRCh37.lupus.chr <- getBM(attributes = c("refsnp_id","chr_name","chrom_start"),
                                      filters = c("snp_filter","chr_name"),
                                      values = list(rsid.lupus %>% filter(!(is.na(CHR_ID))) %>% pull(SNPS),rsid.lupus %>% filter(!(is.na(CHR_ID))) %>% pull(CHR_ID)),
                                      mart=ensemblGRCh37.SNPs)
# some SNPs with chromosome and info position in CHR_ID column (chr6:106564236, chr7:50307334, ...)
snpGRCh37.lupus.nochr <- getBM(attributes = c("refsnp_id","chr_name","chrom_start"), 
                                        filters = c("snp_filter"), 
                                        values = list(rsid.lupus %>% filter(is.na(CHR_ID)) %>% pull(SNPS)), 
                                        mart=ensemblGRCh37.SNPs) %>% filter(chr_name %in% seq(1,22))
# some SNPs with chromosome and info position in SNPS column (6:32449301:32449301, 1:173191475:173191475, ...)
snp.coord <- gsub("chr","",unlist(lapply(strsplit(rsid.lupus %>% filter(grepl("chr",SNPS)) %>% pull(SNPS),":"),
                                         function(x) paste0(x[1],":",x[2],":",x[2]))))
snpGRCh37.lupus.norsid5 <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","allele"), 
                                          filters = c("chromosomal_region"), 
                                          values = snp.coord[1:5], 
                                          mart=ensemblGRCh37.SNPs)
snpGRCh37.lupus.norsid10 <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","allele"), 
                                           filters = c("chromosomal_region"), 
                                           values = snp.coord[6:10], 
                                           mart=ensemblGRCh37.SNPs)
snpGRCh37.lupus.norsid14 <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","allele"), 
                                           filters = c("chromosomal_region"), 
                                           values = snp.coord[11:14], 
                                           mart=ensemblGRCh37.SNPs)
snpGRCh37.lupus.norsid <- rbind(snpGRCh37.lupus.norsid5,snpGRCh37.lupus.norsid10,snpGRCh37.lupus.norsid14) %>% 
  filter(!grepl("[ACGT]{2,}",allele)) %>% dplyr::select(-allele)
# merge all SNPs
snp.lupus <- rbind(snpGRCh37.lupus.chr,
                   snpGRCh37.lupus.nochr,
                   snpGRCh37.lupus.norsid) %>% transmute(snpid=refsnp_id,chr=chr_name,pos=chrom_start) %>% distinct()

## Calculate r2 and check for r2>0.4 or r2>0.8
snp <- rbind(data.frame(snp.yrna.pleio,set="yrnapleio"),data.frame(snp.sjorgen,set="sjorgen"),data.frame(snp.lupus,set="lupus"))
# calculate LD r2
LDres <- list()
for (chrom in seq(1,22)) {
  print(paste0("chr",chrom))
  snp.rsid <- unique(snp %>% filter(chr==chrom) %>% pull(snpid))
  LDres[[paste0("chr",chrom)]] <- LDlinkR::LDmatrix(snp.rsid,pop="EUR",r2d="r2",token="e0b6745d1a68",genome_build = "grch37")
}
LDtotal <- NULL
for (chrom in paste0("chr",seq(1,22))) {
  LDtotal <- rbind(LDtotal,LDres[[chrom]] %>% rename(snp1=RS_number) %>% pivot_longer(-c(snp1),names_to = "snp2", values_to = "r2")) %>%
    filter(snp1!=snp2)
}
# r2>0.8
SNP.LD.08 <- LDtotal %>% filter(r2>0.8) %>% left_join(snp %>% select(snpid,set),by=c("snp1"="snpid"),relationship="many-to-many") %>% rename(set1=set) %>% 
  left_join(snp %>% select(snpid,set),by=c("snp2"="snpid"),relationship="many-to-many") %>% rename(set2=set) %>% filter(set1!=set2) %>%
  filter(!((set1=="lupus" & set2=="sjorgen") | (set1=="sjorgen" & set2=="lupus"))) %>% filter(set1 != "yrnapleio")
ld08 <- length(unique(SNP.LD.08 %>% pull(snp2)))
# 8
write.table(SNP.LD.08,file="Data/SNP_LD_lupus_sjorgen_08.csv",sep=",",row.names = FALSE)
# r2>0.4
SNP.LD.04 <- LDtotal %>% filter(r2>0.4) %>% left_join(snp %>% select(snpid,set),by=c("snp1"="snpid"),relationship="many-to-many") %>% rename(set1=set) %>% 
  left_join(snp %>% select(snpid,set),by=c("snp2"="snpid"),relationship="many-to-many") %>% rename(set2=set) %>% filter(set1!=set2) %>%
  filter(!((set1=="lupus" & set2=="sjorgen") | (set1=="sjorgen" & set2=="lupus"))) %>% filter(set1 != "yrnapleio")
ld04 <- length(unique(SNP.LD.04 %>% pull(snp2)))
# 17
write.table(SNP.LD.04,file="Data/SNP_LD_lupus_sjorgen_04.csv",sep=",",row.names = FALSE)

### 2. Generate random test same size lupus catalog (917) ####
# 1. Download 1000G data 
# https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# 2. Prepare 1000G data for PLINK
#     1000GEURplink.sh
# 3. Build list SNP to generate random SNPs sets (pleiotropic SNPs were removed)
#     plink --bfile EUR_phase3_autosomes --recode --out EUR_phase3_autosomes_snp
#     cut -f 2 EUR_phase3_autosomes_snp.map > snps.map
#     grep -v -w -f SNPpleio.txt snps.map > snpsNOpleio.map 
# 4. Build 1000 randon set 917 SNP
#     randomSetSNPlupus.sh (shuf -n 917 snpsNOpleio.map > randomSetXXXX.txt)

### 3. Calculate how many SNPs in random sets have LD > 0.4 or LD > 0.8  ####
# pleiotropic SNPs nearby RNY (270)
snp.yrna.pleio <- data.frame(snp.yrna.pleio,coord=paste0(snp.yrna.pleio$chr,":",snp.yrna.pleio$pos))
snp.pleio <- data.frame(coord=paste(snp.yrna.pleio$chr,snp.yrna.pleio$pos,sep=":"),snp.yrna.pleio %>% dplyr::select(snpid))
# add info from SNPsnap to pleiotropic SNPs nearby RNY. 
# SNP out from this range (loci_upstream and loci_downstream) it is unlikely to have LD > 0.4 or LD > 0.8
# randomSetXXXX.txt ---> LDrandomXXXX_08.csv/LDrandomXXXX_04.csv
LD08 <- fread("Data/ld0.8_collection_rsID_loci.tab")
LD04 <- fread("Data/ld0.4_collection_rsID_loci.tab")
# https://data.broadinstitute.org/mpg/snpsnap/database_download.html
snp.pleio.loci08 <- snp.pleio %>% left_join(LD08,by=c("snpid"="rsID")) %>% 
  separate(coord,c("chr","posPleio")) %>% dplyr::rename(snpidPleio=snpid)
snp.pleio.loci04 <- snp.pleio %>% left_join(LD04,by=c("snpid"="rsID")) %>% 
  separate(coord,c("chr","posPleio")) %>% dplyr::rename(snpidPleio=snpid)
for (rn in seq(1,10)) {
  print(date())
  print(paste0("Random",rn))
  random <- fread(file=file.path("Output/Fig7b/randomSetLupus",paste0("randomSet",rn,".txt")),
                              col.names = c("chr","snpidRandom","posRandom"),header=FALSE) %>% mutate(chr=as.character(chr))
  # select only SNPs in region [loci_upstream,loci_downstream] for plink running. Reduce time in calculation.
  snp.pleio.random08 <- snp.pleio.loci08 %>% 
    inner_join(random,by=c("chr"), relationship="many-to-many") %>%
    filter(posRandom > loci_upstream & posRandom < loci_downstream) %>% dplyr::select(chr,snpidPleio,snpidRandom)
  snp.pleio.random04 <- snp.pleio.loci04 %>% 
    inner_join(random,by=c("chr"), relationship="many-to-many") %>%
    filter(posRandom > loci_upstream & posRandom < loci_downstream) %>% dplyr::select(chr,snpidPleio,snpidRandom)
  # calculate LD with plink
  print("LD > 0.8")
  outputfile <- file.path("Output/Fig7b/LDrandom",paste0("LDrandom",rn,"_08.csv"))
  if (file.exists(outputfile)) file.remove(outputfile)
  if (nrow(snp.pleio.random08) > 0) {
    for (i in seq(1:nrow(snp.pleio.random08))) {
      plinkfile <- file.path("Data/1000Gplink",paste0("EUR_phase3_chr",snp.pleio.random08[i,"chr"]))
      snp1 <- snp.pleio.random08[i,"snpidPleio"]
      snp2 <- snp.pleio.random08[i,"snpidRandom"]
      bash <- paste0("| grep 'R-sq' | sed 's/   R-sq = /",snp1," ",snp2," /' | cut -d' ' -f1,2,3 | sed 's/ /,/g'")
      plinkrun <- paste("plink --bfile",plinkfile,"--ld",snp1,snp2,bash,">>",outputfile)
      system(plinkrun)
    } 
  } else {
    system(paste0("echo rsXXXX rsXXXX 0.000 > ",outputfile))
  }
  # calculate LD with plink
  print("LD > 0.4")
  outputfile <- file.path("Output/Fig7b/LDrandom",paste0("LDrandom",rn,"_04.csv"))
  if (file.exists(outputfile)) file.remove(outputfile)
  if (nrow(snp.pleio.random04) > 0) {
    for (i in seq(1:nrow(snp.pleio.random04))) {
      plinkfile <- file.path("Data/1000Gplink",paste0("EUR_phase3_chr",snp.pleio.random04[i,"chr"]))
      snp1 <- snp.pleio.random04[i,"snpidPleio"]
      snp2 <- snp.pleio.random04[i,"snpidRandom"]
      bash <- paste0("| grep 'R-sq' | sed 's/   R-sq = /",snp1," ",snp2," /' | cut -d' ' -f1,2,3 | sed 's/ /,/g'")
      plinkrun <- paste("plink --bfile",plinkfile,"--ld",snp1,snp2,bash,">>",outputfile)
      system(plinkrun)
    }
  }
}

###  4. Plot results ####
LD08 <- list()
LD04 <- list()
# SNP pleio near YRNA that correlates with SNP random
for (rn in seq(1,1000)) {
  LD08[[rn]] <- fread(file=file.path("Output/Fig7b/LDrandom",paste0("LDrandom",rn,"_08.csv")),col.names=c("snp1","snp2","r2")) %>% 
    group_by(snp1,snp2) %>% summarise(r2=max(r2)) %>% ungroup() %>% filter(r2>0.8)
  LD04[[rn]] <- fread(file=file.path("Output/Fig7b/LDrandom",paste0("LDrandom",rn,"_04.csv")),col.names=c("snp1","snp2","r2")) %>% 
    group_by(snp1,snp2) %>% summarise(r2=max(r2)) %>% ungroup() %>% filter(r2>0.4)
}
postscript(file="Output/Fig7b/Figure7bLeft.ps")
data.frame(n=unlist(lapply(LD04,function(x) length(unique(x$snp2))))) %>% 
  ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white",bins=18) +
  theme_classic() + 
  ylab("Number of random SNP sets") +
  xlab("Number of SNP pairs r2 > 0.4") +
  geom_segment(aes(x=ld04,y=100,xend=ld04,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off() 
postscript(file="Output/Fig7b/Figure7bRight.ps")
data.frame(n=unlist(lapply(LD08,function(x) length(unique(x$snp2))))) %>%
  ggplot(aes(x=n)) +
  geom_histogram(color="black",fill="white",bins=9) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,8,2)) +
  ylab("Number of random SNP sets") +
  xlab("Number of SNP pairs r2 > 0.8") +
  geom_segment(aes(x=ld08,y=100,xend=ld08,yend=0),color="red",arrow = arrow(length = unit(0.25, "cm"))) +
  theme(text=element_text(size=18))
dev.off() 
