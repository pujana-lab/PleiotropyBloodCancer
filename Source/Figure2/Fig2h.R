# Heatmap showing the overrepresentation and underrepresentation of blood trait-shared variants for each cancer study

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(purrr) # ‘1.0.2’
source("Source/ownFunc.R")

load(file="RData/leadSNP005.RData")
load(file="RData/metadata.RData")
load(file="RData/CancerLabel.RData")
load(file="RData/BloodLabel.RData")
load(file="RData/CancerLevels.RData")
load(file="RData/BloodLevels.RData")

cancer <- unique(leadSNP005 %>% pull(cancer_trait))
blood.cancer <- expand.grid(cancer,blood.levels)
colnames(blood.cancer) <- c("cancer_trait","blood_trait")

# Fischer test
nSNP.cancer.blood <- leadSNP005 %>% left_join(metadata %>% select(id,Phenotype),by = c("blood_trait" = "id")) %>%
  left_join(blood.label,by=c("Phenotype"="LabelInit")) %>% mutate(LabelNew = case_when(
    grepl("BASO",LabelNew) ~ "BASO",
    grepl("EO",LabelNew) ~ "EO",
    grepl("HLSR",LabelNew) ~ "HLSR",
    grepl("LYMPH",LabelNew) ~ "LYMPH",
    grepl("MONO",LabelNew) ~ "MONO",
    grepl("NEUT",LabelNew) ~ "NEUT",
    grepl("RET",LabelNew) ~ "RET",
    TRUE ~ LabelNew
  )) %>% distinct(snpid,cancer_trait,LabelNew) %>% rename(blood_trait=LabelNew) %>%
  count(cancer_trait,blood_trait)

dfc <- nSNP.cancer.blood %>% rename(subgroup=cancer_trait,state=blood_trait)

res.test <- dfContin(dfc) %>% nest(data.contin=c(yy,ny,yn,nn)) %>%
  mutate(m.contin=map(data.contin,~ matrix(c(.x$yy,.x$ny,.x$yn,.x$nn),ncol=2))) %>%
  mutate(tfisher=map(m.contin,~ fisher.test(.x,alternative="two.sided"))) %>%
  mutate(pval=map(tfisher,~.x$p.value)) %>% unnest(pval) %>% select(-data.contin,-m.contin,-tfisher)
df.fisher <- data.frame(res.test,padj=p.adjust(res.test %>% pull(pval),method = "BH")) %>% mutate(signif=case_when(
  padj < 0.1 ~ "*",
  TRUE ~ ""
)) %>% transmute(cancer_trait=subgroup,blood_trait=state,pvalue=pval,padj,signif)

# Chi2 Independence Test
nSNP.cancer.blood.w <- nSNP.cancer.blood %>% pivot_wider(names_from = "cancer_trait", values_from = "n", values_fill=0)
chi2.df <- nSNP.cancer.blood.w %>% select(-blood_trait) %>% as.data.frame()
rownames(chi2.df) <- nSNP.cancer.blood.w %>% pull(blood_trait)
chisq <- chisq.test(chi2.df)

chisq.resid <- data.frame(blood_trait=rownames(chisq$residuals),chisq$residuals) %>% pivot_longer(-blood_trait,names_to="cancer_trait",values_to="resid")

# Heatmap
df.plot <- df.fisher %>% left_join(chisq.resid,by=c("blood_trait","cancer_trait")) %>% 
  left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id")) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>%
  rename(cancer_trait_old=cancer_trait) %>% rename(cancer_trait=LabelNew) %>% select(cancer_trait,blood_trait,pvalue,padj,signif,resid)

postscript(file="Output/Figure2h.ps")
df.plot %>% mutate(cancer_trait=factor(cancer_trait,levels=cancer.levels)) %>%
  ggplot(aes(x = blood_trait , y = cancer_trait, fill = resid)) + 
  geom_tile() +
  geom_point(aes(alpha=signif),shape=0,size=8,stroke=0.4,show.legend = FALSE) +
  scale_fill_gradient2(name="Chi2 resid",low="blue",high="red") +
  xlab("") +
  ylab("") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,size=8,hjust=1,vjust=1),axis.text.y = element_text(size=8)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(size=0.5,colour="black"),
        plot.margin = unit(c(0, 0, 0, 0), "points")) +
  coord_equal()
dev.off()
