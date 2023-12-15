library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

all_tcga_normal <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "HNSC", 
                     "KIRC", "LUAD", "LUSC", "PAAD", 
                     "PRAD", "SKCM", "STAD", "THCA", "UCEC")

ensembl_ids <- c("ENSG00000116747")
gene_names <- c("RO60")

# Function to process RNA-seq data for each TCGA normal dataset
process_rnaseq_data <- function(tcga_name, ensembl_id, gene_name) {
  rnaseq <- readRDS(paste0("./yrna/data/FPKM_UQ.df/TCGA-", tcga_name, ".rds"))
  rnaseq <- as.data.frame(t(rnaseq))
  rnaseq$Sample_ID <- rownames(rnaseq)
  rnaseq <- rnaseq %>% select(Sample_ID, ensembl_id)
  colnames(rnaseq)[2] <- gene_name
  rnaseq <- rnaseq[substr(rnaseq$Sample_ID, 14, 15) == 11,]
  rnaseq
}

# Process and combine RNA-seq data from all TCGA normal datasets
all.rnaseq <- do.call(rbind, lapply(all_tcga_normal, process_rnaseq_data, 
                                    ensembl_id = ensembl_ids, gene_name = gene_names))
all.rnaseq[,gene_names] <- log2(all.rnaseq[,gene_names] + 1)

# Read and merge signature scores
signature_scores <- readRDS("./yrna_mirror_2/data/ssGSEA/normal/Pancancer.rds")
all.data <- merge(all.rnaseq, signature_scores, by = "Sample_ID")

# Prepare data for plotting
all.data <- all.data %>%
  gather(key = "gene", value = "expression", gene_names) %>%
  gather(key = "pleiotropy", value = "score", YRNA_pleio:YRNA_nonpleio, factor_key = F) %>%
  mutate(pleiotropy = ifelse(pleiotropy == "YRNA_pleio", "Yes", "No")) %>%
  filter(expression != 0, gene == "RO60") %>%
  mutate(pleiotropy = factor(pleiotropy, levels = c("Yes", "No")))

# Fit linear models based on Pleiotropic genes (Yes/No)
fits <- all.data %>%
  group_by(pleiotropy) %>%
  do(model = lm(score ~ expression, data = .)) %>%
  summarise(intercept = model$coefficients[[1]], slope = model$coefficients[[2]], gene = "RO60")
fits$pleiotropy = c("Yes", "No")  

title_y <- expression(paste(italic("RNY"), " signature"))
title_x <- expression(paste(italic("RO60"), " (log2 FPKM-UQ)"))
title_normal <- "TCGA normal tissue (n = 593)"

# Generate plot
pdf("Fig7a.pdf")
ggplot(all.data, aes(x = expression, y = score, col = pleiotropy)) +
  geom_point() +
  geom_abline(data = filter(fits, pleiotropy == "Yes"), aes(intercept = intercept, slope = slope), color = "#FF0000", size = 1.5) +
  geom_abline(data = filter(fits, pleiotropy == "No"), aes(intercept = intercept, slope = slope), color = "#4D4D4D", size = 1.5) +
  labs(title = title_normal, x = title_x, y = title_y) +
  theme_classic() +
  scale_color_manual(values = c("#FF0000", "#B2B2B2")) +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = 1, cor.coef.name = "PCC")
dev.off()

