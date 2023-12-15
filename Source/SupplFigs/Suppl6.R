
library(dplyr)
library(tidyr)

library(ComplexHeatmap)
library(circlize)
library(rtracklayer)

library(GSVA)
library(doParallel)
library(ggplot2)
library(ggpubr)

nonpleio_genes <- scan("./data/New_nonpleioRNYs_699.txt", as.character())
pleio_genes <- scan("./data/New_pleioRNYs_118.txt", as.character())

all_rny <- data.table::fread("./data/yrna_seq/RNYs_found.txt")
metadata <- data.table::fread("./data/yrna_seq/yrna_rnaSeq_samples.csv")
countsData <- data.table::fread("./data/yrna_seq/exceRpt_gencode_ReadsPerMillion.txt") %>% as.data.frame()
hamin <- data.frame(Sample_ID = colnames(countsData))
allData <- merge(hamin, metadata, by = "Sample_ID")
rny.genes <- all_rny$transcript_name
countsData <- countsData %>% dplyr::filter(V1 %in% rny.genes)
rownames(countsData) <- gsub("\\:.*","",(countsData$V1))
countsData <- countsData %>% dplyr::select(-V1)
countsData <- log2(countsData+1)
#countsData <- as.data.frame(scale(t(countsData)))
countsData <- as.data.frame((t(countsData)))

###

sel_genes <- data.frame(gene_name = colnames(countsData))
all_rny$gene_name <- gsub("\\:.*","",(all_rny$transcript_name))
found_rny <- merge(sel_genes, all_rny, by = "gene_name")
found_rny$pleio <- ""
found_rny$pleio[found_rny$gene_id %in% nonpleio_genes] <- "non-pleio"
found_rny$pleio[found_rny$gene_id %in% pleio_genes] <- "pleio"
found_rny.pleio <- found_rny %>% dplyr::filter(pleio == "pleio")
found_rny.nonpleio <- found_rny %>% dplyr::filter(pleio == "non-pleio")
found_rny.unknown <- found_rny %>% dplyr::filter(pleio == "")
###

countsData <- countsData %>% dplyr::select(-found_rny.unknown$gene_name)
all.data <- merge(countsData, metadata, by.x=0, by.y="Sample_ID")
rny_vec <- rep("", (ncol(countsData)))
rny_vec[(colnames(countsData) %in% found_rny.pleio$gene_name)] = "pleio"
rny_vec[(colnames(countsData) %in% found_rny.nonpleio$gene_name)] = "non-pleio"

countsData.hep <- countsData[grepl("HEPARIN", rownames(countsData), fixed=T),]
countsData.hep <- as.data.frame(scale(countsData.hep))
countsData.hep <- countsData.hep[,!unname(is.na(countsData.hep)[1,])]
rny_vec_hep <- rep("", (sum((colnames(countsData.hep) %in% found_rny.nonpleio$gene_name))+
                          sum((colnames(countsData.hep) %in% found_rny.pleio$gene_name))))
rny_vec_hep[(colnames(countsData.hep) %in% found_rny.pleio$gene_name)] = "pleio"
rny_vec_hep[(colnames(countsData.hep) %in% found_rny.nonpleio$gene_name)] = "non-pleio"
all.data.hep <- all.data[grepl("HEPARIN", all.data$Row.names, fixed = T),]

ha = HeatmapAnnotation(
  Pleiotropic = rny_vec_hep,
  col = list(Pleiotropic = c("pleio" = "black", 
                             "non-pleio" = "white")
  )
)

ha_row = rowAnnotation(
  Pheno = all.data.hep$Phenotype,
  col = list(Pheno = c("case" = "grey", 
                       "control" = "white"),
             gp = gpar(col = "black")
  )
)

#table(all.data$Phenotype)
#table(all.data$Mut)
col_fun = colorRamp2(c(-2, 0, 2), c("purple", "white", "orange"))

pdf("Suppl6.pdf", height=6, width=8)
Heatmap(as.matrix(countsData.hep), name = "log2 (RPM+1)", 
        column_title = "Y-RNA expression (Heparin samples)",
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        show_column_names = F, 
        show_row_names = T, 
        top_annotation = ha,
        right_annotation = ha_row,
        col = col_fun
)
dev.off()





