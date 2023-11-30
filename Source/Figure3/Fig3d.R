# Histograms showing the percentage of pleiotropic variants identified as eQTL in whole blood (left panel) or 
# immortalized lymphocytes (right panel) of the corresponding candidate genes (X-axis)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’

# SNP pleio near (+/-50kb) target genes (ATM, CHRHR1, HLAs, APOBEC3, ATM, TERT) ####
load(file=file.path(pathR,"leadSNP005.RData"))

genes <- c("TERT","HLA-B","HLA-DPB1","HLA-DQA2","HLA-DQB1","CRHR1","ATM")

ensemblGRCh37.genes <- biomaRt::useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
genes.info <- biomaRt::getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"), 
                             values=genes,mart=ensemblGRCh37.genes)
genes.info <- genes.info %>% filter(chromosome_name %in% seq(1,22)) %>% mutate(chromosome_name=as.integer(chromosome_name))
genes.info <- genes.info %>% mutate(startSNP=start_position-50000,endSNP=end_position+50000)

# Chromosome 22: 39,348,746-39,359,188 (APOBEC3A) Chromosome 22: 39,493,229-39,500,072 (APOBEC3H)
genes.info <- rbind(genes.info,
                    data.frame(hgnc_symbol="APOBEC3",chromosome_name=22,start_position=39348746,end_position=39500072,startSNP=39348746-50000,endSNP=39500072+50000))
glimpse(genes.info)

genes.SNP <- leadSNP005 %>% 
  inner_join(genes.info %>% select(hgnc_symbol,chromosome_name,startSNP,endSNP),by=c("chrnum"="chromosome_name"),relationship="many-to-many") %>% 
  filter(chrpos>startSNP & chrpos<endSNP) %>% select(snpid,chrnum,chrpos,ends_with("trait"),starts_with("zscore"),hgnc_symbol)
# glimpse(genes.SNP %>% distinct(snpid))
save(genes.SNP,file=file.path(pathR,"genesSNP.RData"))

# load(file=file.path(pathR,"genesSNP.RData"))
# load(file="/media/mapardo/SeagateBasic/procure/data/pleiotStudy/R/genesSNP.RData") ###

# rebuild genes.SNP: 1.convert SNP pleio to GRCh38 + 2.add ensgID for target genes ####
box::use(rtracklayer)

load(file=file.path(pathR,"genesSNP.RData"))

chainhg19toHg38 <- rtracklayer$import.chain(file.path(pathExt,"UCSC","hg19ToHg38.over.chain"))
df <- genes.SNP  %>% mutate(chrnum=paste0("chr",chrnum))
grhg19 <- GenomicRanges::GRanges(names=df %>% pull(snpid),seqnames=df %>% pull(chrnum),
                                 ranges=IRanges::IRanges(start=df %>% pull(chrpos),
                                                         end=df %>% pull(chrpos)))
grHg38 <- rtracklayer$liftOver(grhg19, chainhg19toHg38)
SNPsGRCh38 <- as.data.frame(grHg38) %>% distinct(seqnames,start,names)

# save SNP(GRCh38 coord) + target genes (ATM, CHRHR1, HLAs, ...)
genes.SNP.GRCh38 <- genes.SNP %>% left_join(SNPsGRCh38,by=c("snpid"="names")) %>% 
  mutate(variant_id=paste(seqnames,start,sep="_")) %>% select(-seqnames,-start,-chrnum,-chrpos) %>% distinct(hgnc_symbol,variant_id) %>%
  rename(target=hgnc_symbol)
unique(genes.SNP.GRCh38 %>% pull(target)) 
glimpse(genes.SNP.GRCh38)

# assign HGNC symbol & ensg ID to target genes
gene.target <- data.frame(target=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",rep("APOBEC3",7)),
                          gene=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",paste0("APOBEC3",c("A","B","C","D","F","G","H"))))
ensemblGRCh38.genes <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes.name <- biomaRt::getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","ensembl_gene_id","chromosome_name"), 
                             values=gene.target %>% pull(gene),mart=ensemblGRCh38.genes) %>% filter(chromosome_name %in% seq(1,22))
gene.target <- gene.target %>% left_join(genes.name %>% select(-chromosome_name),by=c("gene"="hgnc_symbol"))

genes.SNP.GRCh38 <- genes.SNP.GRCh38 %>% left_join(gene.target,by=c("target"),relationship="many-to-many") %>% rename(gene_id=ensembl_gene_id)
save(genes.SNP.GRCh38,file=file.path(pathR,"genesSNPGRCh38.RData"))

# Proportions Lymph ####
load(file=file.path(pathR,"genesSNPGRCh38.RData"))
load(file=file.path(pathR,"SNPgrpGRCh38All.RData")) ## rescatado de /media/mapardo/SeagateBasic/procure/data/pleiotStudy/R
load(file=file.path(pathR,"eQTLsLymph.RData"))
glimpse(SNPs.grp.gene)
glimpse(genes.SNP.GRCh38)
glimpse(eQTLs.lymp)

SNP.prop <- genes.SNP.GRCh38 %>% left_join(eQTLs.lymp,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(target,eQTL) 
glimpse(SNP.prop)

SNP.grp.prop <- SNPs.grp.gene %>% left_join(eQTLs.lymp,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(eQTL) %>% mutate(target="Random")
glimpse(SNP.grp.prop)

SNP.prop <- rbind(SNP.prop,SNP.grp.prop) %>% mutate(eQTL=factor(eQTL,levels=c("NO","Slope < 0","Slope > 0"))) 

totalSNP <- SNP.prop %>% group_by(target) %>% summarise(total=sum(n)) 

SNP.prop <- SNP.prop %>% left_join(totalSNP,by=c("target")) %>% mutate(porc=n/total)
glimpse(SNP.prop)

gene.level <- c(SNP.prop %>% filter(eQTL=="NO" & target!="Random") %>% arrange(porc) %>% pull(target),"Random")

SNP.prop <- SNP.prop %>% mutate(target=factor(target,levels=gene.level))

# plot Lymph
postscript(file="output/10.eQTLsRegion/SNPgenesGTExLymp.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",size=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()  

# Test two-proportions z-Test Lymph
# H0: PA <= PB
# HA: PA > PB
SNP.prop.random <- SNP.prop %>% filter(target == "Random") %>% select(-porc,-target) %>% rename(nRandom=n,totalRandom=total)
SNP.prop.gene <- SNP.prop %>% filter(target != "Random") 
data.test <- SNP.prop.gene %>% select(-porc) %>% left_join(SNP.prop.random,by=c("eQTL")) %>% rename(nA=n,nTA=total,nB=nRandom,nTB=totalRandom) 
glimpse(data.test)

PropTest <- data.test %>% mutate(pHatA=nA/nTA) %>% mutate(pHatB=nB/nTB) %>% 
  mutate(pHat=(nA+nB)/(nTA+nTB)) %>% mutate(pHatnum=pHat*(1-pHat)) %>%
  mutate(SE=sqrt((pHatnum/nTA)+(pHatnum/nTB))) %>% mutate(zscore=(pHatA-pHatB)/SE) %>%
  mutate(pval=pnorm(zscore,lower.tail=FALSE))
PropTest <- data.frame(PropTest,padj=p.adjust(PropTest %>% pull(pval),method="BH")) %>% mutate(sign = case_when(
  padj < 0.05 ~ "*",
  TRUE ~ ""
))
glimpse(PropTest)

write.table(PropTest,file="output/10.eQTLsRegion/SNPgenesGTExLymp.tsv",sep="\t",row.names = FALSE)

# Proportions Wblood
load(file=file.path(pathR,"genesSNPGRCh38.RData"))
load(file=file.path(pathR,"SNPgrpGRCh38All.RData"))
load(file=file.path(pathR,"eQTLsWBlood.RData"))

SNP.prop <- genes.SNP.GRCh38 %>% left_join(eQTLs.wblood,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(target,eQTL) 
# glimpse(SNP.prop)

SNP.grp.prop <- SNPs.grp.gene %>% left_join(eQTLs.wblood,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(eQTL) %>% mutate(target="Random")
glimpse(SNP.grp.prop)

SNP.prop <- rbind(SNP.prop,SNP.grp.prop) %>% mutate(eQTL=factor(eQTL,levels=c("NO","Slope < 0","Slope > 0"))) 

totalSNP <- SNP.prop %>% group_by(target) %>% summarise(total=sum(n)) 

SNP.prop <- SNP.prop %>% left_join(totalSNP,by=c("target")) %>% mutate(porc=n/total)
glimpse(SNP.prop)

# gene.level <- c(SNP.prop %>% filter(eQTL=="NO" & target!="Random") %>% arrange(porc) %>% pull(target),"Random")

gene.level <- c(SNP.prop %>% filter(eQTL!="NO" & target!="Random") %>% 
                  group_by(target) %>% summarise(pTotal=sum(porc)) %>% arrange(desc(pTotal)) %>% pull(target),
                SNP.prop %>% filter(eQTL=="NO" & porc==1) %>% pull(target),"Random")

SNP.prop <- SNP.prop %>% mutate(target=factor(target,levels=gene.level))

# plot WBlood
postscript(file="output/10.eQTLsRegion/SNPgenesGTExWBlood.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",size=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()  

# Test two-proportions z-Test WBlood
# H0: PA <= PB
# HA: PA > PB
SNP.prop.random <- SNP.prop %>% filter(target == "Random") %>% select(-porc,-target) %>% rename(nRandom=n,totalRandom=total)
SNP.prop.gene <- SNP.prop %>% filter(target != "Random") 
data.test <- SNP.prop.gene %>% select(-porc) %>% left_join(SNP.prop.random,by=c("eQTL")) %>% rename(nA=n,nTA=total,nB=nRandom,nTB=totalRandom) 
glimpse(data.test)

PropTest <- data.test %>% mutate(pHatA=nA/nTA) %>% mutate(pHatB=nB/nTB) %>% 
  mutate(pHat=(nA+nB)/(nTA+nTB)) %>% mutate(pHatnum=pHat*(1-pHat)) %>%
  mutate(SE=sqrt((pHatnum/nTA)+(pHatnum/nTB))) %>% mutate(zscore=(pHatA-pHatB)/SE) %>%
  mutate(pval=pnorm(zscore,lower.tail=FALSE)) 
PropTest <- data.frame(PropTest,padj=p.adjust(PropTest %>% pull(pval),method="BH")) %>% mutate(sign = case_when(
  padj < 0.05 ~ "*",
  TRUE ~ ""
))
glimpse(PropTest)
write.table(PropTest,file="output/10.eQTLsRegion/SNPgenesGTExWBlood.tsv",sep="\t",row.names = FALSE)


# 1000 random genes
# pathaux <- c("/home/mapardo/procure/pleiotStudy")
# load(file=file.path(pathaux,"data","randSNPs1000Genes","LinkDeseq","SNPSelFinal.RData")) ### DONDE ESTAN????
# load(file=file.path(pathaux,"data","randSNPs1000Genes","LinkDeseq","ldSNPs1000Group.RData")) ### DONDE ESTAN???
# SNPs.grp.gene <- data.frame(gene=names(unlist(SNPfinal)),snpid=unlist(SNPfinal)) %>% mutate(gene=substr(gene,1,15)) %>% 
#   left_join(ldSNPs.subgrps %>% select(gene,grp) %>% mutate(gene=substr(gene,1,15)),by=c("gene"))

# ensemblGRCh38.SNPs <- biomaRt::useEnsembl(biomart = "snps",dataset = "hsapiens_snp")
# SNPs.GRCh38 <- biomaRt::getBM(filters= "snp_filter",
#                                                 attributes = c("refsnp_id","chr_name","chrom_start","chrom_end","allele"),
#                                                 values=unique(SNPs.grp.gene %>% pull(snpid)),mart=ensemblGRCh38.SNPs)
# SNPs.GRCh38 <- SNPs.GRCh38 %>% filter(chr_name %in% seq(1,22))
# save(SNPs.GRCh38,file=file.path(pathdata,"R","SNPsGrpGRCh38.RData"))        
# glimpse(SNPs.GRCh38)

# load(file=file.path(pathdata,"R","SNPsGrpGRCh38.RData"))
# SNPs.grp.gene <- SNPs.grp.gene %>% inner_join(SNPs.GRCh38 %>% select(refsnp_id,chr_name,chrom_start),by=c("snpid"="refsnp_id")) %>% 
#   mutate(variant_id=paste0("chr",chr_name,"_",chrom_start)) %>% select(gene,variant_id) %>% rename(gene_id=gene)
# glimpse(SNPs.grp.gene)
# save(SNPs.grp.gene,file=file.path(pathdata,"R","SNPgrpGRCh38All.RData"))