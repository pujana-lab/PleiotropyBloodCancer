# Phylogenetic analysis of RNY sequences from the human genome

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(ggdendro) # ‘0.1.23’
library(data.table) # ‘1.14.8’

# read genetic correlation
gencorr <- fread("Data/stats_gencorr.tsv") %>% as.data.frame()
# read metadata
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
# process data
gencorr <- gencorr %>% left_join(metadata %>% select(id,Phenotype),by = c("p1" = "id")) %>% rename(cancer=Phenotype)
gencorr <- gencorr %>% left_join(metadata %>% select(id,Phenotype),by = c("p2" = "id")) %>% rename(blood=Phenotype)
gencorr <- gencorr %>% mutate(p2 = gsub("_b","",p2)) %>% mutate(p2 = gsub("_f","",p2)) %>% mutate(p2 = gsub("_m","",p2)) %>%
  mutate(logp = -log10(p))

gencorr[gencorr$p1 == 'brca2_oc' & gencorr$p2 == 'NRBCD','rg'] <- NA_real_

gencorr <- gencorr %>% mutate(sig = case_when(
  p < 0.05 ~ "*",
  TRUE ~ ""))

# Cluster Blood
blood <- sort(unique(gencorr$blood))
cancer <- sort(unique(gencorr$cancer))
gencorr.blood <- gencorr %>% arrange(blood,cancer) %>% mutate(rg = case_when(
  is.na(rg) ~ 0,
  TRUE ~ rg))# Se ordena primero blood y luego cancer
split.blood <- split(gencorr.blood$rg,cancer)
matrix.blood <- as.matrix(bind_cols(split.blood))
rownames(matrix.blood) <- blood
d.blood = dist(matrix.blood,method="euclidean")
h.blood = hclust(d.blood, method ="complete")
dendro.blood <- as.dendrogram(h.blood)

# Cluster Cancer
gencorr.cancer <- gencorr %>% arrange(cancer,blood) %>% mutate(rg = case_when(
  is.na(rg) ~ 0,
  TRUE ~ rg))# Se ordena primero cancer
split.cancer <- split(gencorr.cancer$rg,blood)
matrix.cancer <- as.matrix(bind_cols(split.cancer))
rownames(matrix.cancer) <- cancer
d.cancer = dist(matrix.cancer,method="euclidean")
h.cancer = hclust(d.cancer, method ="complete")
dendro.cancer <- as.dendrogram(h.cancer)

cancer.traits <- cancer[order.dendrogram(dendro.cancer)]
cancer.traits <- as.data.frame(cbind(cancer.traits,seq(1,length(cancer.traits))))
colnames(cancer.traits) <- c("cancer","ncancer")
gencorr <- gencorr %>% left_join(cancer.traits,by = c("cancer" = "cancer")) %>% mutate(ncancer = as.integer(ncancer))

blood.traits <- blood[order.dendrogram(dendro.blood)]
blood.traits <- as.data.frame(cbind(blood.traits,seq(1,length(blood.traits))))
colnames(blood.traits) <- c("blood","nblood")
gencorr <- gencorr %>% left_join(blood.traits,by = c("blood" = "blood")) %>% mutate(nblood = as.integer(nblood))

glimpse(gencorr)

# Heatmap
heatmap.gencorr <- gencorr %>% ggplot(aes(x = nblood , y = ncancer, color = rg, size = logp, label=sig)) +
  geom_point(shape=15) +
  scale_size(name="-log10(p)",range = c(0,4.8),breaks = c(0.5,1,1.5,2.0,2.5,3.0)) +
  scale_color_gradient2(name="Genetic Corr.",limits=c(-0.3,0.3),
                        low="blue",high="red",na.value="white",
                        oob=squish,breaks =c(-0.2,0,0.2)) +
  guides(color=guide_colourbar(ticks=FALSE)) +
  theme(legend.position = c(1.6, 0.5)) +
  xlab("") +
  ylab("") +
  geom_text(size=3,color="black") +
  scale_x_continuous(breaks = seq(1, dim(blood.traits)[1],1) ,minor_breaks = 0.5 + seq(1, dim(blood.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(blood.traits)[1]),labels = blood.traits$blood,
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(1, dim(cancer.traits)[1],1), minor_breaks = 0.5 + seq(1, dim(cancer.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(cancer.traits)[1]),labels = cancer.traits$cancer,
                     expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90,size=7,hjust=1),axis.text.y = element_text(size=7)) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black",size = 0.1,linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()

dendro.blood.plot <- ggdendro::ggdendrogram(dendro.blood) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), text = element_text(size=0.01),
        plot.margin=unit(c(0,0,0,0), 'cm'))

dendro.cancer.plot <- ggdendro::ggdendrogram(dendro.cancer, rotate = TRUE) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), text = element_text(size=0.01),
        plot.margin=unit(c(0,0,0,0), 'cm'))

postscript("Output/FigureS2_1.ps")
grid::grid.newpage()
print(heatmap.gencorr, vp = grid::viewport(x = 0.35, y = 0.416, width = 0.8, height = 0.9))
print(dendro.cancer.plot, vp = grid::viewport( x = 0.738, y = 0.5525, width = 0.15, height = 0.637))
print(dendro.blood.plot, vp = grid::viewport( x = 0.46, y = 0.929, width = 0.425, height = 0.146))
dev.off()

# Genetic correlation Cancer Traits
gencorr.cancer <- fread("Data/stats_gencorr_cancer_fdr.tsv") %>% as.data.frame()
gencorr.cancer <- gencorr.cancer %>% left_join(metadata %>% select(id,Phenotype),by = c("p1" = "id")) %>% rename(cancer1=Phenotype)
gencorr.cancer <- gencorr.cancer %>% left_join(metadata %>% select(id,Phenotype),by = c("p2" = "id")) %>% rename(cancer2=Phenotype)
gencorr.cancer <- gencorr.cancer %>% mutate(FDR=case_when(
  FDR==0 ~ 1e-300,
  TRUE ~ FDR
)) %>% mutate(logFDR = -log10(FDR)) %>% mutate(p=case_when(
  p==0 ~ 1e-300,
  TRUE ~ p
)) %>% mutate(logp = -log10(p))

gencorr.cancer <- gencorr.cancer %>% mutate(sigFDR = case_when(
  FDR < 0.05 ~ "*",
  TRUE ~ "")) %>% mutate(sigp = case_when(
    p < 0.05 ~ "*",
    TRUE ~ ""))

cancer.traits <- data.frame(unique(gencorr.cancer %>% arrange(cancer1) %>% pull(cancer1)),seq(1:28))
colnames(cancer.traits) <- c("cancer","ncancer")
gencorr.cancer <- gencorr.cancer %>% left_join(cancer.traits,by = c("cancer1" = "cancer")) %>% 
  mutate(ncancer1 = as.integer(ncancer)) %>% select(-ncancer)
gencorr.cancer <- gencorr.cancer %>% left_join(cancer.traits,by = c("cancer2" = "cancer")) %>% 
  mutate(ncancer2 = as.integer(ncancer)) %>% select(-ncancer)

postscript(file="Output/FigureS2_2.ps")
gencorr.cancer %>% mutate(logFDR = case_when(
  logFDR > 5 ~ 5,
  TRUE ~ logFDR
)) %>% ggplot(aes(x = ncancer1 , y = ncancer2, color = rg, size = logFDR, label=sigFDR)) +
  geom_point(shape=15) +
  scale_size(name="-log10(FDR)",range = c(0,5),breaks = c(1,2,3,4,5), labels = c("1","2","3","4",">5")) +
  scale_color_gradient2(name="Genetic Corr.",limits=c(-1,1),
                        low="blue",high="red",na.value="white",
                        oob=scales::squish,breaks =c(-1,0,1)) +
  guides(color=guide_colourbar(ticks=FALSE)) +
  xlab("") +
  ylab("") +
  geom_text(size=3,color="black") +
  scale_x_continuous(breaks = seq(1, dim(cancer.traits)[1],1) ,minor_breaks = 0.5 + seq(1, dim(cancer.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(cancer.traits)[1]),labels = cancer.traits$cancer,
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(1, dim(cancer.traits)[1],1), minor_breaks = 0.5 + seq(1, dim(cancer.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(cancer.traits)[1]),labels = cancer.traits$cancer,
                     expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90,size=7,hjust=1),axis.text.y = element_text(size=7)) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black",size = 0.1,linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()
dev.off()

# Genetic correlation Blood Traits
gencorr.blood <- fread("Data/stats_gencorr_blood_fdr.tsv") %>% as.data.frame()
gencorr.blood <- gencorr.blood %>% left_join(metadata %>% select(id,Phenotype),by = c("p1" = "id")) %>% rename(blood1=Phenotype)
gencorr.blood <- gencorr.blood %>% left_join(metadata %>% select(id,Phenotype),by = c("p2" = "id")) %>% rename(blood2=Phenotype)
gencorr.blood <- gencorr.blood %>% mutate(FDR=case_when(
  FDR==0 ~ 1e-300,
  TRUE ~ FDR
)) %>% mutate(logFDR = -log10(FDR)) %>% mutate(p=case_when(
  p==0 ~ 1e-300,
  TRUE ~ p
)) %>% mutate(logp = -log10(p))

gencorr.blood <- gencorr.blood %>% mutate(sigFDR = case_when(
  FDR < 0.05 ~ "*",
  TRUE ~ "")) %>% mutate(sigp = case_when(
    p < 0.05 ~ "*",
    TRUE ~ ""))

blood.traits <- data.frame(unique(gencorr.blood %>% arrange(blood1) %>% pull(blood1)),seq(1:27))
colnames(blood.traits) <- c("blood","nblood")
gencorr.blood <- gencorr.blood %>% left_join(blood.traits,by = c("blood1" = "blood")) %>% 
  mutate(nblood1 = as.integer(nblood)) %>% select(-nblood)
gencorr.blood <- gencorr.blood %>% left_join(blood.traits,by = c("blood2" = "blood")) %>% 
  mutate(nblood2 = as.integer(nblood)) %>% select(-nblood)

postscript("Output/FigureS2_3.ps")
gencorr.blood %>% mutate(logFDR = case_when(
  logFDR > 5 ~ 5,
  TRUE ~ logFDR
)) %>% ggplot(aes(x = nblood1 , y = nblood2, color = rg, size = logFDR, label=sigFDR)) +
  geom_point(shape=15) +
  scale_size(name="-log10(FDR)",range = c(0,5),breaks = c(1,2,3,4,5), labels = c("1","2","3","4",">5")) +
  scale_color_gradient2(name="Genetic Corr.",limits=c(-1,1),
                        low="blue",high="red",na.value="white",
                        oob=scales::squish,breaks =c(-1,0,1)) +
  guides(color=guide_colourbar(ticks=FALSE)) +
  xlab("") +
  ylab("") +
  geom_text(size=3,color="black") +
  scale_x_continuous(breaks = seq(1, dim(blood.traits)[1],1) ,minor_breaks = 0.5 + seq(1, dim(blood.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(blood.traits)[1]),labels = blood.traits$blood,
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(1, dim(blood.traits)[1],1), minor_breaks = 0.5 + seq(1, dim(blood.traits)[1]-1, 1),
                     limits =c(0.5, 0.5 + dim(blood.traits)[1]),labels = blood.traits$blood,
                     expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90,size=7,hjust=1),axis.text.y = element_text(size=7)) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black",size = 0.1,linetype = "solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()
dev.off()