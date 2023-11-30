# Histograms showing the percentage of pleiotropic variants identified as eQTL in whole blood (left panel) or 
# immortalized lymphocytes (right panel) of the corresponding candidate genes (X-axis)

# TODO README Z-Test proportion SNP pleio selected genes in GTEx lymph and wblood (SNP random genes as background) FIGURE 3C

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(biomaRt) # ‘2.56.1’
library(rtracklayer) # ‘1.60.1’

load(file="RData/leadSNP005.RData")
load(file="RData/SNPgrpGRCh38All.RData") ## TODO rescatado de /media/mapardo/SeagateBasic/procure/data/pleiotStudy/R
load(file="RData/eQTLsLymph.RData") 
load(file="RData/eQTLsWBlood.RData")

## SNP pleio near (+/-50kb) target genes
# ATM, CHRHR1, HLAs, APOBEC3, ATM, TERT
genes <- c("TERT","HLA-B","HLA-DPB1","HLA-DQA2","HLA-DQB1","CRHR1","ATM")

ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
genes.info <- getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"), 
                             values=genes,mart=ensemblGRCh37.genes) %>% filter(chromosome_name %in% seq(1,22)) %>% 
  mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  mutate(startSNP=start_position-50000,endSNP=end_position+50000)
# Chromosome 22: 39,348,746-39,359,188 (APOBEC3A) Chromosome 22: 39,493,229-39,500,072 (APOBEC3H)
genes.info <- rbind(genes.info,
                    data.frame(hgnc_symbol="APOBEC3",chromosome_name=22,start_position=39348746,end_position=39500072,startSNP=39348746-50000,endSNP=39500072+50000))
# SNP pleio
genes.SNP <- leadSNP005 %>% 
  inner_join(genes.info %>% dplyr::select(hgnc_symbol,chromosome_name,startSNP,endSNP),by=c("chrnum"="chromosome_name"),relationship="many-to-many") %>% 
  filter(chrpos>startSNP & chrpos<endSNP) %>% dplyr::select(snpid,chrnum,chrpos,ends_with("trait"),starts_with("zscore"),hgnc_symbol)

## rebuild genes.SNP: 1.convert SNP pleio to GRCh38 + 2.add ensgID for target genes
chainhg19toHg38 <- import.chain("Data/hg19ToHg38.over.chain")
df <- genes.SNP  %>% mutate(chrnum=paste0("chr",chrnum))
grhg19 <- GenomicRanges::GRanges(names=df %>% pull(snpid),seqnames=df %>% pull(chrnum),
                                 ranges=IRanges::IRanges(start=df %>% pull(chrpos),
                                                         end=df %>% pull(chrpos)))
grHg38 <- liftOver(grhg19, chainhg19toHg38)
SNPsGRCh38 <- as.data.frame(grHg38) %>% distinct(seqnames,start,names)
# SNP(GRCh38 coord) + target genes (ATM, CHRHR1, HLAs, ...)
genes.SNP.GRCh38 <- genes.SNP %>% left_join(SNPsGRCh38,by=c("snpid"="names")) %>% 
  mutate(variant_id=paste(seqnames,start,sep="_")) %>% dplyr::select(-seqnames,-start,-chrnum,-chrpos) %>% distinct(hgnc_symbol,variant_id) %>%
  dplyr::rename(target=hgnc_symbol)

## assign HGNC symbol & ensg ID to target genes
gene.target <- data.frame(target=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",rep("APOBEC3",7)),
                          gene=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",paste0("APOBEC3",c("A","B","C","D","F","G","H"))))
ensemblGRCh38.genes <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes.name <- getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","ensembl_gene_id","chromosome_name"), 
                             values=gene.target %>% pull(gene),mart=ensemblGRCh38.genes) %>% filter(chromosome_name %in% seq(1,22))
gene.target <- gene.target %>% left_join(genes.name %>% dplyr::select(-chromosome_name),by=c("gene"="hgnc_symbol"))
genes.SNP.GRCh38 <- genes.SNP.GRCh38 %>% left_join(gene.target,by=c("target"),relationship="many-to-many") %>% dplyr::rename(gene_id=ensembl_gene_id)


# Proportions Lymph
# TODO SNPs.grp.gene 
glimpse(genes.SNP.GRCh38)
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
postscript(file="Output/Figure3dRight.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",size=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()  

# Proportions Wblood
SNP.prop <- genes.SNP.GRCh38 %>% left_join(eQTLs.wblood,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(target,eQTL) 

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

gene.level <- c(SNP.prop %>% filter(eQTL!="NO" & target!="Random") %>% 
                  group_by(target) %>% summarise(pTotal=sum(porc)) %>% arrange(desc(pTotal)) %>% pull(target),
                SNP.prop %>% filter(eQTL=="NO" & porc==1) %>% pull(target),"Random")

SNP.prop <- SNP.prop %>% mutate(target=factor(target,levels=gene.level))

# plot WBlood
postscript(file="Output/Figure3dWBlood.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",size=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()  
