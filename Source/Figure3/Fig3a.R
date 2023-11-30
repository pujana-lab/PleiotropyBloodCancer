# Pie charts showing the contribution of pleotropic variants in telomere length-associated gene loci across the cancer studies

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(RColorBrewer) # ‘1.1.3’

load(file="RData/telomereGene005.RData")

cancer.blood <- telomereSNP005 %>% select(V12,V13) %>% rename(cancer=V12,blood=V13) %>% mutate(cancer=case_when(
  cancer == "Luminal A breast cancer" ~ "BC#1 LumA",
  cancer == "Breast cancer overall (BCAC)" ~ "BC#1",
  cancer == "Breast cancer overall (UKBB)" ~ "BC#2",
  cancer == "TNBC" ~ "BC#1 TNBC",
  cancer == "HER2 breast cancer" ~ "BC#1 HER2",
  cancer == "Luminal B - HER2 negative breast cancer" ~ "BC#1 LumB/HER2-",
  cancer == "Cervical" ~ "Cervix",
  TRUE ~ cancer
)) %>% mutate(cancer=factor(cancer,levels=c("BC#1","BC#2","BC#1 LumA","BC#1 LumB/HER2-","BC#1 HER2","BC#1 TNBC","BRCA1 breast cancer","BRCA1 TNBC",
                                            "Bladder","Cervix","Melanoma","NHL","Prostate","Rectum")))

# piechart Telomere SNP cancer
cancer.prop <- cancer.blood %>% count(cancer) %>% 
  mutate(total=sum(n)) %>% mutate(prop=n/total) %>% select(cancer,n,prop)

pie.colors <- c(brewer.pal(8,"Dark2"),
                brewer.pal(11,"Set3")[seq(1,11)])

postscript(file="Output/Figure3a.ps")
cancer.blood %>% count(cancer) %>%
  ggplot(aes(x="",y=n,fill=cancer)) +
  geom_bar(stat="identity",width=1) +
  coord_polar(theta="y") + 
  scale_fill_manual(values=pie.colors) +
  theme(axis.line=element_blank(), axis.ticks = element_blank(),panel.border = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(), axis.text = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(),legend.title = element_blank())
dev.off()
