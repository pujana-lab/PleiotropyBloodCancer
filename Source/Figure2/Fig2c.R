# Histogram depicting the distribution of classes of genetic elements across the identified pleiotropic loci and cancer studies

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’
library(biomaRt) # ‘2.56.1’

# association genes/pleiotropic SNPs (gene2SNP.R)
load(file="RData/leadSNP005Gene.RData")
# read metadata
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
# cancer label 
cancer.label <- fread(file="Data/cancerLabel.csv")

# cancer traits order
cancer.sort <- leadSNP005.gene %>% count(cancer_trait,gene) %>% filter(!is.na(gene)) %>% count(cancer_trait) %>% 
  dplyr::rename(numGenes = n) %>% left_join(metadata %>% dplyr::select(id,Phenotype),by = c("cancer_trait" = "id"))  %>% dplyr::select(Phenotype,numGenes) %>% 
  arrange(desc(numGenes)) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  dplyr::rename(oldLabel=Phenotype) %>% dplyr::rename(Phenotype=LabelNew) %>% pull(Phenotype)
# gene type statistics
pleio.cancer.gene <- leadSNP005.gene %>% count(cancer_trait,gene,gene_type) %>% filter(!is.na(gene)) %>% 
  count(cancer_trait,gene_type) %>% dplyr::rename(numGenes = n) %>% left_join(metadata %>% dplyr::select(id,Phenotype),by = c("cancer_trait" = "id"))  %>%
  dplyr::select(Phenotype,gene_type,numGenes) %>% mutate(gene_type = case_when(
    gene_type == "protein_coding" ~ "Protein coding",
    grepl("pseudogene",gene_type) ~ "Pseudogene",
    gene_type == "protein_coding" ~ "Protein coding",
    gene_type == "antisense" ~ "Antisense",
    gene_type == "processed_transcript" ~ "Processed transcript",
    gene_type == "lincRNA" ~ "linc-RNA",
    gene_type == "sense_overlapping" ~ "Sense overlapping",
    gene_type == "miRNA" ~ "mi-RNA",
    TRUE ~ NA_character_)) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  dplyr::rename(oldLabel=Phenotype) %>% dplyr::rename(Phenotype=LabelNew)
# biotype order
level.biotype <- pleio.cancer.gene %>% group_by(gene_type) %>% summarize(total=sum(numGenes)) %>% arrange(desc(total)) %>% pull(gene_type)
# plot
postscript(file=file.path("Output/Figure2c.ps"))
pleio.cancer.gene %>%
  mutate(gene_type = factor(gene_type,levels=level.biotype)) %>%
  mutate(Phenotype = factor(Phenotype,levels=cancer.sort)) %>% 
  ggplot(aes(x=Phenotype, y=numGenes, fill=gene_type)) +
  geom_bar(stat="identity") +
  scale_fill_discrete(name="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=70,hjust=1)) + 
  xlab("") +
  ylab("Shared loci (n) conjFDR < 0.05")
dev.off()
