# Venn diagram showing the overlap between mouse gene orthologs that, when mutated, 
# cause immune system alterations (MP:0005387; "immune system phenotype") 
# and the pleiotropic gene set (all cancers included)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’
library(readxl) # ‘1.4.3’
library(ggVennDiagram) # ‘1.2.3’

load(file=file.path("RData/leadSNP005Gene.RData"))

# mammalian phenotype data immu and infl (excel MAPujana)
mouse.pheno.immu <- as.data.frame(read_excel(path="Data/Mouse-phenotypes_immune-inflammation.xlsx", 
                                                     sheet="Immune",skip=0,col_names=c("MP","info1","info2"),col_types=rep("text",3))) %>% select(MP)

mouse.pheno.infl <- as.data.frame(read_excel(path="Data/Mouse-phenotypes_immune-inflammation.xlsx", 
                                                     sheet="Inflammation",skip=0,col_names=c("MP","info1","info2"),col_types=rep("text",3))) %>% select(MP)

mamma.pheno <- fread("Data/VOC_MammalianPhenotype.rpt", header=FALSE,col.names=c("MP","term","description")) %>% as.data.frame()
# http://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt

mamma.pheno.immu <- mamma.pheno %>% filter(grepl("immune",term) | grepl("immune",description))
mamma.pheno.infl <- mamma.pheno %>% filter(grepl("inflammation|inflammatory",term) | grepl("inflammation|inflammatory",description))

# Human/Mouse linked genes for immune and inflamation term MGI
# read human/mouse gene link
gene.link <- data.table::fread("Data/HMD_HumanPhenotype.rpt") %>% as.data.frame()
# http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
colnames(gene.link) <- c("geneNAME","geneEntrez","geneName","MGIid","linkMP","F")
gene.link <- gene.link %>% select(geneNAME,linkMP) %>% mutate(linkMP=gsub(" ","",linkMP)) %>% separate_rows(linkMP,sep=",")

# Genes linkados con immu/infl
gene.link.immu <- gene.link %>% filter(linkMP %in% mamma.pheno.immu$MP)
gene.link.infl <- gene.link %>% filter(linkMP %in% mamma.pheno.infl$MP)

# Gene link to SNPs pleio
pleio.gene.immu <- leadSNP005.gene %>% filter(gene %in% unique(gene.link.immu$geneNAME))
pleio.gene.infl <- leadSNP005.gene %>% filter(gene %in% unique(gene.link.infl$geneNAME))

# Fisher's exact test
Total <- length(unique(gene.link$geneNAME))
group1 <- length(unique(leadSNP005.gene %>% filter(gene %in% unique(gene.link$geneNAME)) %>% pull(gene)))
group2 <- length(unique(gene.link.immu$geneNAME))
Overlap <- length(unique(pleio.gene.immu$gene))

res.fisher <- fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, Total-group2-group1 +Overlap), 2, 2), alternative='greater')

res.fisher$p.value
# 3.518338e-08
res.fisher$estimate
# odds ratio 
# 1.453306

# Venn diagram
venn.data <- list(Pleio = unique(leadSNP005.gene %>% filter(gene %in% unique(gene.link$geneNAME)) %>% pull(gene)),
                  Immune = unique(gene.link.immu$geneNAME))

postscript(file="Output/Figure4e.ps")
ggVennDiagram(venn.data, label_alpha = 0,category.names = c("Pleiotropic","Immune"),label="count") +
  scale_fill_continuous(low="white",high="white") + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) +
  ggtitle("OR = 1.45, p = 3.5e-8")
dev.off()
