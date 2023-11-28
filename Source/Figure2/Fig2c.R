rm(list=ls())
gc()

library(dplyr) # [1] ‘1.1.3’
library(tidyr) # [1] ‘1.3.0’
library(ggplot2) # [1] ‘3.4.3’

load(file="Data/leadSNP005Gene.RData")
load(file="Data/metadata.RData")
load(file="Data/CancerLabel.RData")
load(file="Data/BloodLabel.RData")

cancer.sort <- leadSNP005.gene %>% count(cancer_trait,gene) %>% filter(!is.na(gene)) %>% count(cancer_trait) %>% 
  rename(numGenes = n) %>% left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id"))  %>% select(Phenotype,numGenes) %>% 
  arrange(desc(numGenes)) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew) %>% pull(Phenotype)

pleio.cancer.gene <- leadSNP005.gene %>% count(cancer_trait,gene,gene_type) %>% filter(!is.na(gene)) %>% 
  count(cancer_trait,gene_type) %>% rename(numGenes = n) %>% left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id"))  %>%
  select(Phenotype,gene_type,numGenes) %>% mutate(gene_type = case_when(
    gene_type == "protein_coding" ~ "Protein coding",
    grepl("pseudogene",gene_type) ~ "Pseudogene",
    gene_type == "protein_coding" ~ "Protein coding",
    gene_type == "antisense" ~ "Antisense",
    gene_type == "processed_transcript" ~ "Processed transcript",
    gene_type == "lincRNA" ~ "linc-RNA",
    gene_type == "sense_overlapping" ~ "Sense overlapping",
    gene_type == "miRNA" ~ "mi-RNA",
    TRUE ~ NA_character_)) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew)

level.biotype <- pleio.cancer.gene %>% group_by(gene_type) %>% summarize(total=sum(numGenes)) %>% arrange(desc(total)) %>% pull(gene_type)

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