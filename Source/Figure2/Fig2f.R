# Histogram depicting the distribution of classes of genetic elements across the identified pleiotropic loci and blood traits

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’

# association genes/pleiotropic SNPs (gene2SNP.R)
load(file="RData/leadSNP005Gene.RData")
# read metadata
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
# blood label 
blood.label <- fread(file="Data/bloodLabel.csv")

# blood traits order
blood.sort <- leadSNP005.gene %>% left_join(metadata %>% select(id,Phenotype),by = c("blood_trait" = "id")) %>% 
  count(Phenotype,gene) %>% filter(!is.na(gene)) %>% count(Phenotype) %>% rename(numGenes = n) %>% 
  select(Phenotype,numGenes) %>% arrange(desc(numGenes)) %>% left_join(blood.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew) %>% pull(Phenotype)
# gene type statistics
pleio.blood.gene <- leadSNP005.gene %>% left_join(metadata %>% select(id,Phenotype),by = c("blood_trait" = "id")) %>% 
  count(Phenotype,gene,gene_type) %>% filter(!is.na(gene)) %>% count(Phenotype,gene_type) %>% 
  rename(numGenes = n) %>% mutate(gene_type = case_when(
    gene_type == "protein_coding" ~ "Protein coding",
    grepl("pseudogene",gene_type) ~ "Pseudogene",
    gene_type == "protein_coding" ~ "Protein coding",
    gene_type == "antisense" ~ "Antisense",
    gene_type == "processed_transcript" ~ "Processed transcript",
    gene_type == "lincRNA" ~ "linc-RNA",
    gene_type == "sense_overlapping" ~ "Sense overlapping",
    gene_type == "miRNA" ~ "mi-RNA",
    TRUE ~ NA_character_)) %>% left_join(blood.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew)
# biotype order
level.biotype <- pleio.blood.gene %>% group_by(gene_type) %>% summarize(total=sum(numGenes)) %>% arrange(desc(total)) %>% pull(gene_type)
# plot
postscript(file="Output/Figure2f.ps")
pleio.blood.gene %>%
  mutate(gene_type = factor(gene_type,levels=level.biotype)) %>% 
  mutate(Phenotype = factor(Phenotype,levels=blood.sort)) %>% 
  ggplot(aes(x=Phenotype, y=numGenes, fill=gene_type)) +
  geom_bar(stat="identity") +
  scale_fill_discrete(name="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=70,hjust=1)) + 
  xlab("") +
  ylab("Shared loci (n) conjFDR < 0.05")
dev.off()
