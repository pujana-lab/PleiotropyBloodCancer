
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

metadata <- data.table::fread("./data/yrna_rnaSeq_samples_metadata.csv")
countsData <- data.table::fread("./data/exceRpt_gencode_ReadsPerMillion.txt") %>% as.data.frame()

countsData <- data.table::fread("./data/exceRpt_miRNA_ReadsPerMillion.txt") %>% as.data.frame()
mirna_list <- c("hsa-miR-150-5p", "hsa-miR-21-5p", "hsa-miR-122-5p", "hsa-miR-16-5p")
countsData <- countsData %>% dplyr::filter(V1 %in% mirna_list)
rownames(countsData) <- countsData$V1
countsData <- countsData %>% dplyr::select(-V1)
countsData <- log2(countsData+1)
countsData <- as.data.frame(scale(t(countsData)))

countsData$Sample_ID <- rownames(countsData)
all.data <- merge(countsData, metadata, by = "Sample_ID")

#######

brca.data <- all.data %>% dplyr::filter(Mut %in% c("BRCA1", "BRCA2"))
brca.data <- brca.data %>% dplyr::select(-ng_uL, -Mut, -Sample_ID)

brca.data.long <- gather(brca.data, miRNA, log2RPM, 
                         `hsa-miR-21-5p`:`hsa-miR-150-5p`)

pdf("Suppl8a.pdf")
ggboxplot(brca.data.long, x = "Phenotype", y = "log2RPM",
          xlab = "Phenotype", ylab = "log2(RPM+1)",
          title = "BRCA1/2 mutated samples miRNA expression affected vs unaffected",
          color = "Phenotype", palette = "jco",
          add = "jitter") + 
  facet_wrap(~miRNA) +
  stat_compare_means(method = "wilcox.test", label.x.npc = 0.4, label.y.npc = 0.95)
dev.off()

#######

hepa.data <- all.data %>% dplyr::filter(Phenotype %in% c("case", "control"))
hepa.data <- hepa.data %>% dplyr::select(-ng_uL, -Mut, -Sample_ID)

hepa.data.long <- gather(hepa.data, miRNA, log2RPM, 
                         `hsa-miR-21-5p`:`hsa-miR-150-5p`)

pdf("Suppl8b.pdf")
ggboxplot(hepa.data.long, x = "Phenotype", y = "log2RPM",
          xlab = "Phenotype", ylab = "log2(RPM+1)",
          title = "Heparin samples miRNA expression case vs control",
          color = "Phenotype", palette = "jco",
          add = "jitter") + 
  facet_wrap(~miRNA) +
  stat_compare_means(method = "wilcox.test", label.x.npc = 0.4, label.y.npc = 0.95)
dev.off()



