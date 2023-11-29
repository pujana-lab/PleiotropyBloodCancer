# Histogram depicting the number of variants shared by blood traits and cancer risk

library(dplyr) # [1] ‘1.1.3’
library(tidyr) # [1] ‘1.3.0’
library(ggplot2) # [1] ‘3.4.3’

load(file="Data/leadSNP005.RData")
load(file="Data/metadata.RData")
load(file="Data/CancerLabel.RData")
load(file="Data/BloodLabel.RData")

pleio.blood.loci <- leadSNP005 %>% mutate(blood_trait = substr(blood_trait,1,nchar(blood_trait)-2)) %>% 
  group_by(blood_trait,cancer_trait,locusnum,chrnum) %>% 
  summarize(start = min(chrpos), end = max(chrpos), cFDRmin = min(conjfdr), cFDRmax = max(conjfdr)) %>%
  ungroup() %>% arrange(blood_trait,chrnum,start) %>% 
  mutate(last_start = lag(start), last_end = lag(end), last_chrnum = lag(chrnum), last_blood_trait = lag(blood_trait)) %>%
  filter(is.na(last_end) | last_blood_trait != blood_trait | chrnum != last_chrnum | start > last_end ) %>%
  count(blood_trait) %>% rename(numLoci = n)
glimpse(pleio.blood.loci)

pleio.blood.snp <- leadSNP005 %>% count(blood_trait,snpid) %>% count(blood_trait) %>% rename(numSNPs = n) %>% 
  left_join(metadata %>% select(id,Phenotype,Total),by = c("blood_trait" = "id"))  %>% filter(grepl("_b",blood_trait)) %>%
  mutate(Phenotype=gsub("_b","",Phenotype)) %>% select(Phenotype,numSNPs,Total) %>% left_join(blood.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(oldLabel=Phenotype) %>% rename(Phenotype=LabelNew) %>% left_join(pleio.blood.loci,by=c("Phenotype"="blood_trait"))

postscript(file=file.path("Output/Figure2e.ps"))
pleio.blood.snp %>% arrange(desc(numSNPs)) %>% 
  mutate(Phenotype = factor(Phenotype,levels=pleio.blood.snp %>% arrange(desc(numSNPs)) %>% pull(Phenotype))) %>% 
  mutate(lognumSNPs = log10(numSNPs)) %>%
  ggplot(aes(x=Phenotype, y=numSNPs/1000)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=numLoci),vjust=-0.5) +
  geom_col() +
  theme_classic() +
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle=70,hjust=1,size=13),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=14)) + 
  xlab("") +
  ylab(bquote('Shared SNPs '~(nx10^3)~'conjFDR < 0.05'))
dev.off()