# Pie charts showing the contribution of each blood trait to each cancer risk study based on the number of shared variants

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(cowplot) # ‘1.1.1’
library(grid) # ‘4.3.2’

# read pleiotropic SNP list (pleio_loci.tsv). Filter to conjfdr < 0.05
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# read metadata
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
# cancer label 
cancer.label <- fread(file="Data/cancerLabel.csv")
# blood label
blood.label <- fread(file="Data/bloodLabel.csv")
# cancer levels
cancer.levels <- c("Thyroid","Rectum","Prostate","Pancreas","Oropharyngeal","NHL","Melanoma","Lung","Leukemia","Kidney",
                   "Gastroesophageal","Endometrium","Colon","Cervix","Bladder","Ovary","BRCA2 OC","BRCA1 OC","BC#2","BRCA2 BC",
                   "BRCA1 BC","BRCA1 TNBC","BC#1 TNBC","BC#1 HER2+","BC#1 LumB/HER2-","BC#1 LumB","BC#1 LumA","BC#1")
# blood levels
blood.levels <- c("HCT","HGB","MCHC","NRBCD","RBC#","RDW","MCV","HLSR#","HLSR%","IRF","RET#","RET%","MPV","PCT","PDW","PLT#",
                  "BASO#","BASO%","EO#","EO%","LYMPH#","LYMPH%","MONO#","MONO%","NEUT#","NEUT%","WBC#")

# number shared variants each cancer trait/blood trait
nSNP.cancer.blood <- leadSNP005 %>% count(cancer_trait,blood_trait) %>%
  left_join(metadata %>% select(id,Phenotype),by = c("blood_trait" = "id"))  %>%
  left_join(blood.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(blood_trait_old=blood_trait) %>% rename(blood_trait=LabelNew) %>% select(-Phenotype) %>% 
  left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id")) %>% left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(cancer_trait_old=cancer_trait) %>% rename(cancer_trait=LabelNew) %>% select(cancer_trait,blood_trait,n) %>% 
  pivot_wider(names_from="cancer_trait",values_from="n",values_fill = 0) %>% pivot_longer(-blood_trait,names_to="cancer_trait",values_to="n") 
# plot
pie.colors <- c(RColorBrewer::brewer.pal(9,"Reds")[seq(3,9)],
                RColorBrewer::brewer.pal(9,"YlOrRd")[seq(1,5)],
                RColorBrewer::brewer.pal(9,"RdPu")[seq(5,8)],
                RColorBrewer::brewer.pal(9,"RdPu")[c(1,3)],
                RColorBrewer::brewer.pal(9,"Paired")[c(1,2)],
                RColorBrewer::brewer.pal(9,"Paired")[c(3,4)],
                RColorBrewer::brewer.pal(9,"BrBG")[c(2,1)],
                RColorBrewer::brewer.pal(9,"BrBG")[c(7,9)],
                c("grey50"))
plots <- list()
for (cancer in rev(cancer.levels)) {
  plots[[cancer]] <- nSNP.cancer.blood %>% filter(cancer_trait==cancer) %>%
    ggplot(aes(x="",y=n,fill=blood_trait)) +
    ggtitle(cancer) +
    geom_bar(stat="identity",width=1) +
    coord_polar(theta="y") + 
    scale_fill_manual(values=pie.colors) +
    theme(axis.line=element_blank(), axis.ticks = element_blank(),panel.border = element_blank(),
          panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(), axis.text = element_blank(),
          panel.background = element_blank(),legend.position = "none", axis.title = element_blank())
}
postscript(file="Output/Figure2g.ps",paper="special",width=10,height=20,onefile=FALSE,horizontal=FALSE)
ggpubr::ggarrange(plotlist=plots,ncol = 4, nrow = 7)
dev.off()
# plot legend in a separated file
plot <- nSNP.cancer.blood %>% filter(cancer_trait==cancer) %>%
  ggplot(aes(x="",y=n,fill=blood_trait)) +
  geom_bar(stat="identity",width=1) +
  coord_polar(theta="y") + 
  scale_fill_manual(values=pie.colors) +
  theme(axis.line=element_blank(), axis.ticks = element_blank(),panel.border = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor = element_blank(), axis.text = element_blank(),
        panel.background = element_blank(), axis.title = element_blank())
legend <- get_legend(plot)
postscript(file="Output/Figure2gLegend.ps")
grid.newpage()
grid.draw(legend)
dev.off()
