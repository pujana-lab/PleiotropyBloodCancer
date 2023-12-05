# Histogram showing the relative contribution of pleiotropic variants (%; Y-axis) 
# in RNY-containing loci (± 50 kb centered on each variant) across cancer studies (X-axis)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(biomaRt) # ‘2.56.1’
library(ggplot2) # ‘3.4.3’

# read pleiotropic SNP list (pleio_loci.tsv). Filter to conjfdr < 0.05
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# read metadata
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
# cancer label 
cancer.label <- fread(file="Data/cancerLabel.csv")
# cancer levels
cancer.levels <- c("Thyroid","Rectum","Prostate","Pancreas","Oropharyngeal","NHL","Melanoma","Lung","Leukemia","Kidney",
                   "Gastroesophageal","Endometrium","Colon","Cervix","Bladder","Ovary","BRCA2 OC","BRCA1 OC","BC#2","BRCA2 BC",
                   "BRCA1 BC","BRCA1 TNBC","BC#1 TNBC","BC#1 HER2+","BC#1 LumB/HER2-","BC#1 LumB","BC#1 LumA","BC#1")

# RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                                 attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                                "start_position","end_position","gene_biotype","strand"), 
                                 mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
# assignment of RNY to nearby SNP
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
n.yrna <- length(unique(yrna.snp %>% pull(ensembl_gene_id)))
n.snp <- length(unique(yrna.snp %>% pull(snpid)))
# percentage SNPs mapped to YRNA
SNP <- NULL
for (cancer in unique(leadSNP005 %>% pull(cancer_trait))) {
  n <- length(unique(leadSNP005 %>% filter(cancer_trait == cancer) %>% pull(snpid)))
  n <- c(n,length(unique(yrna.snp %>% filter(cancer_trait == cancer) %>% pull(snpid))))
  SNP <- rbind(SNP,n)
}
n <- length(unique(leadSNP005 %>% pull(snpid)))
n <- c(n,length(unique(yrna.snp %>% pull(snpid))))
SNP <- rbind(SNP,n)
rownames(SNP) <- c(unique(leadSNP005 %>% pull(cancer_trait)),"All")
colnames(SNP) <- c("total","YRNA")
numSNP.yrna.cancer <- as.data.frame(SNP) %>% mutate(cancer_trait=rownames(SNP)) %>% 
  pivot_longer(-cancer_trait,names_to = "type",values_to = "SNPs")
percSNP.yrna.cancer <- numSNP.yrna.cancer %>% 
  left_join(numSNP.yrna.cancer %>% filter(type == "total") %>% dplyr::select(-type) %>% rename(total=SNPs),by=c("cancer_trait")) %>%
  filter(type != "total") %>% mutate(perc=SNPs/total*100)
percSNP.yrna.cancer <- rbind(percSNP.yrna.cancer,percSNP.yrna.cancer %>% mutate(perc = 100-perc) %>% mutate(type="Others"))
percSNP.yrna.cancer <- percSNP.yrna.cancer %>% left_join(metadata %>% dplyr::select(id,Phenotype),by=c("cancer_trait"="id")) %>% 
  left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% dplyr::select(-Phenotype) %>% rename(Phenotype=LabelNew)  %>% 
  mutate(type=factor(type,levels=c("Others","YRNA"))) %>% mutate_if(is.character,coalesce,"All") %>%
  mutate(Phenotype=factor(Phenotype, levels=c(rev(cancer.levels),"All")))
# plot
postscript(file="Output/Figure5a.ps")
percSNP.yrna.cancer %>%
  ggplot(aes(x=Phenotype,y=perc,fill=type)) +
  geom_bar(stat="identity",color="black",linewidth=0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=70,hjust=1)) + 
  xlab("") +
  ylab("Percentage pleiotropic SNPs") +
  theme(axis.text.x=element_text(size=10,angle=45,hjust=1,vjust=1),axis.text.y=element_text(size=10),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) + 
  scale_fill_manual(values=c("orange","white"),name="",labels=c("YRNA","Others"),limits=c("YRNA","Others")) 
dev.off()
