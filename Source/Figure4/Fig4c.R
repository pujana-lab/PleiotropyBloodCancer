# Uniform Manifold Approximation and Projection (UMAP) of the pleiotropic gene signature expression (score indicated in inset) 
# in the bone marrow single-cell RNA sequencing profiles

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’
library(biomaRt) # ‘2.56.1’
library(Seurat) # ‘4.3.0.1’

# background genes DICE (select genes TPM > 1)
DICE.tpm <- fread("Data/mean_tpm_merged.csv") %>% as.data.frame()
# https: https://dice-database.org/download/mean_tpm_merged.csv
DICE.tpm.mat <- as.matrix(DICE.tpm %>% dplyr::select(-gene))
rownames(DICE.tpm.mat) <- DICE.tpm %>% pull(gene)
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37" ,dataset = "hsapiens_gene_ensembl")
DICE.biomart37.genes <- getBM(filters="ensembl_gene_id",values=rownames(DICE.tpm.mat),
                                       attributes = c("hgnc_symbol","ensembl_gene_id","gene_biotype"),
                                       mart=ensemblGRCh37.genes)
DICE.tpm.expr <- rownames(DICE.tpm.mat)[(apply(DICE.tpm.mat > 1,1,all))]
DICE.tpm.genes <- DICE.biomart37.genes %>% filter(ensembl_gene_id %in% DICE.tpm.expr) %>% pull(hgnc_symbol)
# read WTA_projected (select genes DICE TPM>1)
scRNA.data <- readRDS("Data/WTA_projected.rds")
# https://figshare.com/articles/dataset/Expression_of_97_surface_markers_and_RNA_transcriptome_wide_in_13165_cells_from_a_healthy_young_bone_marrow_donor/13397987
scRNA.data.DICE <- subset(x=scRNA.data,features = DICE.tpm.genes)
WTAgenesDICE <- rownames(scRNA.data.DICE)
# association genes/pleiotropic SNPs (gene2SNP.R)
load(file="RData/leadSNP005Gene.RData")

# select genes associated to pleiotropic SNPs within WTA selection + random gene set
genes.leadSNP <- unique(leadSNP005.gene %>% filter(!is.na(gene) & gene_biotype == "protein_coding") %>% pull(gene))
genes.leadSNP.DICE <- genes.leadSNP[genes.leadSNP %in% WTAgenesDICE]
gene.set <- list()
gene.set[["leadSNP"]] <- genes.leadSNP.DICE
for (i in seq(1,100)) {
  gene.set[[paste0("Random",i)]] <-sample(WTAgenesDICE,1088,replace=FALSE)  
}
# run score calculation
wta.leadSNP <- AddModuleScore(scRNA.data.DICE,
                                      features = gene.set,
                                      name="WTA")
# UMAP plot + scores
postscript(file="Output/Figure4c.ps")
FeaturePlot(wta.leadSNP,
                    features = "WTA1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
dev.off()
