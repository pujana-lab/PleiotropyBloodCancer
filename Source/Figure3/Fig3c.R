# Histogram showing the relative contribution of pleiotropic variants in defined genomic regions or hotspots (bottom inset) across the cancer studies

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’

load(file="RData/leadSNP005.RData")
load(file="RData/metadata.RData")
load(file="RData/CancerLabel.RData")
load(file="RData/CancerLevels.RData")
load(file="RData/regionSNP.RData")

snp.region <- list()
for (r in regionSNPs %>% pull(region)) {
  c <- regionSNPs %>% filter(region==r) %>% pull(chr)
  s <- regionSNPs %>% filter(region==r) %>% pull(start)
  e <- regionSNPs %>% filter(region==r) %>% pull(end)
  if (r %in% c("TERT","HLA")) {
    snp.region[[r]] <- unique(leadSNP005 %>% filter(chrnum==c,between(chrpos,s,e)) %>% pull(snpid))
  } else {
    snp.region[[r]] <- unique(leadSNP005 %>% filter(chrnum==c,between(chrpos,s-50000,e+50000)) %>% pull(snpid))
  }
}

SNP <- NULL
for (cancer in unique(leadSNP005 %>% pull(cancer_trait))) {
  n <- length(unique(leadSNP005 %>% filter(cancer_trait == cancer) %>% pull(snpid)))
  for (r in names(snp.region)) {
    n <- c(n,sum(snp.region[[r]] %in% unique(leadSNP005 %>% filter(cancer_trait == cancer) %>% pull(snpid))))
  }
  SNP <- rbind(SNP,n)
}
rownames(SNP) <- unique(leadSNP005 %>% pull(cancer_trait))
colnames(SNP) <- c("total",names(snp.region))

numSNP.region.cancer <- as.data.frame(SNP) %>% mutate(cancer_trait=rownames(SNP)) %>% 
  pivot_longer(-cancer_trait,names_to = "region",values_to = "SNPs")

numSNP.totalSNP <- numSNP.region.cancer %>% filter(region=="total") %>% select(-region) %>% rename(totalSNPs=SNPs)
numSNP.region.cancer <- numSNP.region.cancer %>% filter(region!="total")
numSNP.inRegion <- numSNP.region.cancer %>% group_by(cancer_trait) %>% summarise(inSNPs=sum(SNPs))
numSNP.totalSNP <- numSNP.totalSNP %>% left_join(numSNP.inRegion,by=c("cancer_trait")) %>% mutate(outSNPs=totalSNPs-inSNPs)

numSNP.region.cancer <- rbind(numSNP.region.cancer,
                              data.frame(region="Others",numSNP.totalSNP) %>% rename(SNPs=outSNPs) %>% select(cancer_trait,region,SNPs)) 
numSNP.region.cancer <- numSNP.region.cancer %>% left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id")) %>% 
  left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(cancer_trait_old=cancer_trait) %>% rename(cancer_trait=LabelNew)

numSNP.totalSNP <- numSNP.totalSNP %>% left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id")) %>% 
  left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(cancer_trait_old=cancer_trait) %>% rename(cancer_trait=LabelNew)

region.levels <- c("Others",rev(regionSNPs %>% arrange(chr,start) %>% pull(label)))
percSNP.region.cancer <- numSNP.region.cancer %>% left_join(numSNP.totalSNP %>% select(cancer_trait,totalSNPs),by=c("cancer_trait")) %>%
  mutate(perc=SNPs/totalSNPs*100) %>% select(cancer_trait,region,perc) %>% 
  mutate(cancer_trait=factor(cancer_trait, levels=rev(cancer.levels))) %>% 
  left_join(regionSNPs %>% select(region,label),by=c("region")) %>%
  mutate(label=case_when(is.na(label) ~ "Others",
                         TRUE ~ label)) %>% 
  mutate(label=factor(label,levels=region.levels))

coul <- RColorBrewer::brewer.pal(12, "Paired")
coul <-c(coul,"#F0027F","white")

postscript(file="Output/Figure3c.ps")
percSNP.region.cancer %>%
  ggplot(aes(x=cancer_trait,y=perc,fill=label)) +
  geom_bar(stat="identity",color="black",size=0.1) + ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) + 
  scale_fill_manual(values=coul,name="",limits=rev(region.levels))
dev.off()