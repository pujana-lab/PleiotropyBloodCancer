# Genomic diagram showing the relative position of the pleiotropic variants (dots) across human chromosomes 1–22 (X-axis) and cancer-risk studies (Y-axis)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’

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
# number bases at each chr (start and end for each chr in an acumulative way)
chrlen <- c(249250621,243199373,198022430,191154276,180915260,171115067,
            159138663,146364022,141213431,135534747,135006516,133851895,	
            115169878,107349540,102531392,90354753,81195210,78077248,
            59128983,63025520,48129895,51304566)
names(chrlen) <- seq(1:22)
chrend <- cumsum(chrlen)
chrstart <- c(0,chrend[1:21]) + 1 
chrcum <- as.data.frame(cbind(names(chrlen),chrstart,chrend))
colnames(chrcum) <- c("chrnum","chrstart","chrend")
chrcum <- chrcum %>% mutate(chrnum = as.integer(chrnum), chrstart = as.double(chrstart), chrend = as.double(chrend))

# asignment number 1..28 for cancer trait (Support representation at y-axis)
cancer.num <- as.data.frame(cbind(cancer.levels,seq(1:length(cancer.levels))))
colnames(cancer.num) <- c("cancer_trait","cancernum")
cancer.num <- cancer.num %>% mutate(cancernum = as.integer(cancernum))

# asignment number 1..27 for blod trait
blood.levels <- blood.label %>% pull(LabelNew)
blood.num <- data.frame(blood_trait=blood.label %>% pull(LabelNew),bloodnum=seq(1:dim(blood.label)[1]))

# plot
gap <- c(10000000) # gap between chr for avoiding overlap
pleio.chrcum <- leadSNP005 %>% left_join(metadata %>% select(id,Phenotype),by = c("cancer_trait" = "id"))  %>%
  left_join(cancer.label,by=c("Phenotype"="LabelInit")) %>% 
  rename(cancer_trait_old=cancer_trait) %>% rename(cancer_trait=LabelNew) %>% # change cancer trait id
  mutate(blood_trait = substr(blood_trait,1,nchar(blood_trait)-2)) %>% # modify blood trait id
  inner_join(cancer.num,by = c("cancer_trait")) %>%
  inner_join(blood.num,by = c("blood_trait")) %>% 
  mutate(posy = (cancernum-1)*length(blood.levels) + bloodnum + (cancernum-1)) %>% # calculate y position
  inner_join(chrcum, by = c("chrnum" = "chrnum")) %>% 
  mutate(chrstartgap = chrstart + (chrnum-1)*gap) %>% # add gap between chr (chrstart y chrend)
  mutate(chrendgap = chrend + (chrnum-1)*gap) %>% mutate(chrcum = chrpos + chrstartgap)  # calculate position of SNPs
postscript(file="Output/Figure3b.ps",height=10,width=30)
pleio.chrcum %>% ggplot(aes(x=chrcum,y=posy,color=cancer_trait)) +
  geom_hline(yintercept = seq(1,(length(cancer.levels)-1))*(length(blood.levels)+1), color = "lightblue", linewidth=0.1) +
  scale_y_continuous(label=cancer.levels,
                     breaks=length(blood.levels)/2 + seq(0,length(cancer.levels)-1)*length(blood.levels) + seq(0,length(cancer.levels)-1),
                     expand=c(0,3),position="right") +
  geom_vline(xintercept = gap/2 + chrcum$chrend + (seq(1,22)-1)*gap, color = "lightblue", linewidth=0.1) +
  scale_x_continuous(label=seq(1:22),
                     breaks= chrcum %>% mutate(labpos = chrstart + (chrnum-1)*gap + (chrend-chrstart)/2) %>% pull(labpos),
                     expand=c(0,gap/2),position="top") +
  geom_point(show.legend=FALSE, size = 0.1) +
  labs(x = "Chromosome", y = NULL) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_rect(colour="black",fill=NA,size=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 10, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(colour=c("black")))
dev.off()
