# Violin plot of the expression level of the pleiotropic and non-pleiotropic 
# RNY signatures in blood plasma from cancer patients and healthy individuals, 
# as indicated on the X-axis 

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(GSVA) # ‘1.48.3’
library(purrr) # ‘1.0.2’
library(ggsignif) # ‘0.6.4’
library(data.table) # ‘1.14.8’
library(biomaRt) # ‘2.56.1’

##  Read expression data for GSE71008 (exRNA_Patel_tushar)
# https://exrna-atlas.org/exat/datasets/EXR-TPATE1OqELFf-AN#EXRTPATE1OqELFfAN
# Data info
bios.files <- list.files(path="Data/Patel_Tushar/biosampleMetadata",pattern="*metadata.tsv")
bios.info <- NULL
for (file.name in bios.files) {
  bio <- fread(file.path("Data/Patel_Tushar/biosampleMetadata",file.name))
  bios.info <- rbind(bios.info,as.data.frame(bio[bio$"#property" %in% c("- Donor ID","- Name","-- Disease Type"),]) %>% pull(value))
}
colnames(bios.info) <- c("sample","file","tumor")
bios.info <- as.data.frame(bios.info)

donor.files <- list.files(path="Data/Patel_Tushar/donorMetadata",pattern="*metadata.tsv")
donor.info <- NULL
for (file.name in donor.files) {
  donor.info <- rbind(donor.info,as.data.frame(fread(file.path("Data/Patel_Tushar/donorMetadata",file.name))) %>% pull(value))
}
colnames(donor.info) <- c("file","type","status","sex","age")
donor.info <- as.data.frame(donor.info) %>% mutate(age = as.integer(gsub(" y","",age)))

info <- bios.info %>% left_join(donor.info,by=c("file"))

# Counts data
sample.counts <- list()
sample.stats <- NULL
count <- c(0)
ensemblID.list <- NULL
for (id in info$sample) {
  count <- count + 1
  print(count)
  print(id)
  stats <- fread(file.path("Data/Patel_Tushar/datafiles",paste("sample",id,"fastq.stats",sep="_")),skip=c("#"))
  stats <- data.frame(sample=id,stats)
  sample.stats <- rbind(sample.stats,stats)
  counts <- fread(file.path("Data/Patel_Tushar/datafiles",paste("sample",id,"fastq",sep="_"),"readCounts_gencode_sense.txt")) %>%
    separate(ReferenceID,c("ensemblID","biotype","gene"),c(":")) %>% mutate(ensemblIDv=ensemblID) %>% 
    mutate(ensemblID=substr(ensemblIDv,1,15))
  sample.counts[[id]] <- data.frame(sample=id,antisense=FALSE,counts)
  counts <- fread(file.path("Data/Patel_Tushar/datafiles",paste("sample",id,"fastq",sep="_"),"readCounts_gencode_antisense.txt")) %>%
    separate(ReferenceID,c("ensemblID","biotype","gene"),c(":")) %>% mutate(ensemblIDv=ensemblID) %>% 
    mutate(ensemblID=substr(ensemblIDv,1,15))
  if (dim(counts)[1] != 0) {
    counts <- data.frame(sample=id,antisense=TRUE,counts)
    sample.counts[[id]] <- rbind(sample.counts[[id]],counts)
  }
  ensemblID.list <- c(ensemblID.list,sample.counts[[id]] %>% pull(ensemblID))
}
ensemblID.list <- unique(ensemblID.list)

# Ensembl gene ID for Ensembl transcript ID annotation
ensembl.genes <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
biomart.exrna <- getBM(filters="ensembl_transcript_id",values=ensemblID.list,
                                attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name","chromosome_name",
                                               "start_position","end_position","transcript_length",
                                               "gene_biotype","strand"), 
                                mart=ensembl.genes) %>% filter(chromosome_name %in% seq(1,22)) 

# Build expression matrix (genes/samples)
rnaseq.matrix <- NULL
for (id in names(sample.counts)) {
  ReadCount <- sample.stats %>% filter(sample==id & Stage=="input") %>% pull(ReadCount)
  rnaseq.matrix <- cbind(rnaseq.matrix,biomart.exrna %>% left_join(sample.counts[[id]] %>% filter(antisense==FALSE) %>% dplyr::select(ensemblID,uniqueReadCount) ,
                                                                   by=c("ensembl_transcript_id"="ensemblID")) %>% 
                           replace(is.na(.),0) %>% mutate(uniqueReadNorm=log10(uniqueReadCount/((ReadCount/1000000)*transcript_length)+1)) %>% 
                           dplyr::select(ensembl_transcript_id,uniqueReadNorm) %>% pull(uniqueReadNorm))
}
rownames(rnaseq.matrix) <- biomart.exrna %>% pull(ensembl_transcript_id)
colnames(rnaseq.matrix) <- names(sample.counts)

## Build pleiotropic and non-pleiotropic RNY signatures
load(file="RData/leadSNP005.RData")
# RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
# Assignment of RNY to nearby SNP
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
# Group RNY: pleiotropic and non-pleiotropic
yrna.pleio <- yrna.snp %>% distinct(ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,gene_biotype,strand)
yrna.nonpleio <- biomart37.yrnarny %>% filter(!(ensembl_gene_id %in% yrna.pleio$ensembl_gene_id))
yrna.class <- list(pleio=yrna.pleio$ensembl_gene_id,nonpleio=yrna.nonpleio$ensembl_gene_id)
genesign <- list()
genesign[["without"]] <- biomart.exrna %>% filter(ensembl_gene_id %in% yrna.class[["nonpleio"]]) %>% pull(ensembl_transcript_id)
genesign[["with"]] <- biomart.exrna %>% filter(ensembl_gene_id %in% yrna.class[["pleio"]]) %>% pull(ensembl_transcript_id)

## Calculate signatures
ssGSEAs <- gsva(as.matrix(rnaseq.matrix), genesign,
                     method = "ssgsea",
                     kcdf = "Gaussian",
                     ssgsea.norm = T)
ssGSEAs.result <- as.data.frame(ssGSEAs) %>% mutate(SNPs=rownames(ssGSEAs)) %>% 
  pivot_longer(-SNPs,names_to="sample",values_to="sign") %>% left_join(info,by=c("sample")) %>%
  mutate(type=factor(type,levels=c("Control","Experimental"))) %>% 
  mutate(SNPs=factor(SNPs,levels=c("with","without"))) %>%
  mutate(tumor=factor(tumor,levels=c("Colon Carcinoma","Pancreatic Carcinoma","Prostate Carcinoma","Healthy Control")))

## Violin plot Colon Carcinoma/Pancreatic Carcinoma/Prostate Carcinoma/HealthyControl + Wilcoxon Test
wt.results <- ssGSEAs.result %>% nest(d=-c(tumor)) %>% mutate(wt=map(d,~ wilcox.test(sign ~ SNPs,data=.x))) %>%
  mutate(pval=map(wt,~ .x$p.value)) %>% unnest(pval) %>% transmute(feature=tumor,pval)

comp <- list(paste("Colon Carcinoma",c("with","without"),sep=":"),
             paste("Pancreatic Carcinoma",c("with","without"),sep=":"),
             paste("Prostate Carcinoma",c("with","without"),sep=":"),
             paste("Healthy Control",c("with","without"),sep=":")) 
annot <- c("p < 2.2e-16","p = 0.002","p < 2.2e-16", "p = 1.5e-16")
postscript(file="Output/Figure6i.ps")
ssGSEAs.result %>%
  ggplot(aes(x=tumor:SNPs,y=sign)) + 
  geom_violin(aes(color=SNPs),trim=FALSE) +
  scale_color_manual(values=c("red","grey")) +
  scale_x_discrete(labels=c("Colon/+SNPs","Colon/-SNPs","Pancreatic/+SNPs","Pancreatic/-SNPs",
                            "Prostate/+SNPs","Prostate/-SNPs","Control/+SNPs","Control/-SNPs"))+
  geom_jitter(aes(color=SNPs),shape=16,position=position_jitter(0.1),size=0.5) +
  geom_boxplot(width=0.05,outlier.alpha=0) +
  theme_classic() +
  xlab("") +
  ylab("RNY signature") +
  geom_signif(comparisons=comp,annotations=annot,textsize=4,margin_top=c(0.2))
dev.off()
