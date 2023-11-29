# Plot depicting the relationship between the number of individuals in each GWAS analyzed and the number of identified pleiotropic variants

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(ggrepel) # ‘0.9.3’

load(file="RData/leadSNP005.RData")
load(file="RData/metadata.RData")
load(file="RData/CancerLabel.RData")
load(file="RData/BloodLabel.RData")

pleio.cancer.loci <- leadSNP005 %>% group_by(cancer_trait,blood_trait,locusnum,chrnum) %>% 
  summarize(start = min(chrpos), end = max(chrpos), cFDRmin = min(conjfdr), cFDRmax = max(conjfdr)) %>%
  ungroup() %>% arrange(cancer_trait,chrnum,start) %>% 
  mutate(last_start = lag(start), last_end = lag(end), last_chrnum = lag(chrnum), last_cancer_trait = lag(cancer_trait)) %>%
  filter(is.na(last_end) | last_cancer_trait != cancer_trait | chrnum != last_chrnum | start > last_end ) %>%
  count(cancer_trait) %>% rename(numLoci = n)

pleio.cancer.snp <- leadSNP005 %>% count(cancer_trait,snpid) %>% count(cancer_trait) %>% rename(numSNPs = n) %>% 
  left_join(pleio.cancer.loci,by=c("cancer_trait")) %>%
  left_join(metadata %>% select(id,Phenotype,Cases),by = c("cancer_trait" = "id"))  %>%
  select(Phenotype,numSNPs,numLoci,Cases) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew)

postscript(file=file.path("Output/Figure2d.ps"))
pleio.cancer.snp  %>%
  mutate(lognumSNPs = log10(numSNPs),logCases=log10(Cases)) %>%
  ggplot(aes(x=logCases, y=lognumSNPs, label=Phenotype)) +
  geom_point() + geom_text_repel() +
  theme_classic() + theme(text = element_text(size=14)) +
  ylab("Shared SNPs conjFDR < 0.05 (log10(n))") + xlab("Individuals (log10(n))")
dev.off()