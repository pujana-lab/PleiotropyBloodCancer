rm(list=ls())
gc()

library(dplyr) # [1] ‘1.1.3’
library(tidyr) # [1] ‘1.3.0’
library(ggplot2) # [1] ‘3.4.3’

load(file="Data/leadSNP005.RData")
load(file="Data/metadata.RData")
load(file="Data/CancerLabel.RData")

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

postscript(file=file.path("Output/Figure2b.ps"))
pleio.cancer.snp  %>% arrange(desc(numSNPs)) %>% 
  mutate(Phenotype = factor(Phenotype,levels=pleio.cancer.snp %>% arrange(desc(numSNPs)) %>% pull(Phenotype))) %>% 
  mutate(lognumSNPs = log10(numSNPs)) %>%
  ggplot(aes(x=Phenotype, y=numSNPs/1000, fill=Cases/1000)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=numLoci),vjust=-0.5) +
  geom_col() +
  viridis::scale_fill_viridis(name=bquote("Number cases"~ (x10^3)),option="D",limits = c(0,150)) +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle=70,hjust=1,size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=14)) + 
  xlab("") +
  ylab(bquote('Shared SNPs '~(nx10^3)~'conjFDR < 0.05'))
dev.off()