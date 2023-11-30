# 

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’

# Seurat scRNAseq data WTA ####
# https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/
# background genes DICE
DICE.tpm <- data.table::fread(file.path(pathExt,"DICE","mean_tpm_merged.csv")) %>% as.data.frame()
DICE.tpm.mat <- as.matrix(DICE.tpm %>% select(-gene))
rownames(DICE.tpm.mat) <- DICE.tpm %>% pull(gene)
save(DICE.tpm.mat,file=file.path(pathR,"DICE_TPMmatrix.RData"))
ensemblGRCh37.genes <- biomaRt::useEnsembl(biomart = "genes", version="GRCh37" ,dataset = "hsapiens_gene_ensembl")
DICE.biomart37.genes <- biomaRt::getBM(filters="ensembl_gene_id",values=rownames(DICE.tpm.mat),
                                       attributes = c("hgnc_symbol","ensembl_gene_id","gene_biotype"),
                                       mart=ensemblGRCh37.genes)
save(DICE.biomart37.genes,file=file.path(pathR,"DICEGRCh37genes.RData"))

# select genes TPM > 1
load(file=file.path(pathR,"DICE_TPMmatrix.RData"))
load(file=file.path(pathR,"DICEGRCh37genes.RData"))

DICE.tpm.expr <- rownames(DICE.tpm.mat)[(apply(DICE.tpm.mat > 1,1,all))]
DICE.tpm.genes <- DICE.biomart37.genes %>% filter(ensembl_gene_id %in% DICE.tpm.expr) %>% pull(hgnc_symbol)
length(DICE.tpm.genes)
# TPM > 1
# 10871
save(DICE.tpm.genes,file=file.path(pathR,"DICEgenesExpr_TPM1.RData"))

# select genes DICE(TPM>1) within WAT
load(file=file.path(pathR,"DICEgenesExpr_TPM1.RData"))
# https://figshare.com/articles/dataset/Expression_of_97_surface_markers_and_RNA_transcriptome_wide_in_13165_cells_from_a_healthy_young_bone_marrow_donor/13397987
scRNA.data <- readRDS(file.path(pathExt,"FigShare","WTA_projected.rds"))
scRNA.data.DICE <- subset(x=scRNA.data,features = DICE.tpm.genes)
save(scRNA.data.DICE,file=file.path(pathR,"WTAdataDICE_TPM1.RData"))
WTAgenesDICE <- rownames(scRNA.data.DICE)
length(WTAgenesDICE)
# TPM > 5
# 6422
# TPM > 1
# 8574
save(WTAgenesDICE,file=file.path(pathR,"WTAgenesDICE_TPM1.RData"))

# select genes lead SNP within WTA selection + random gene set
load(file=file.path(pathR,"WTAdataDICE_TPM1.RData"))
load(file=file.path(pathR,"WTAgenesDICE_TPM1.RData"))
load(file=file.path(pathR,"leadSNP005Gene.RData"))
genes.leadSNP <- unique(leadSNP005.gene %>% filter(!is.na(gene) & gene_biotype == "protein_coding") %>% pull(gene))
# length(genes)
# 2097
genes.leadSNP.DICE <- genes.leadSNP[genes.leadSNP %in% WTAgenesDICE]
length(genes.leadSNP.DICE)
# TPM > 1
# 640

gene.set <- list()
gene.set[["leadSNP"]] <- genes.leadSNP.DICE
for (i in seq(1,100)) {
  gene.set[[paste0("Random",i)]] <-sample(WTAgenesDICE,1088,replace=FALSE)  
}
head(lapply(gene.set,length))

# run score calculation
wta.leadSNP <- Seurat::AddModuleScore(scRNA.data.DICE,
                                      features = gene.set,
                                      name="WTA")

score.WTA <- wta.leadSNP[[]]
save(score.WTA,file=file.path(pathR,"ScoreWTADICE_TPM1.RData"))

# UMAP plot + scores
load(file=file.path(pathR,"ScoreWTADICE_TPM1.RData"))
glimpse(score.WTA)

postscript(file="output/14.scRNAimmu/scRNASeq_SeuratDICE_TPM1.ps")
Seurat::FeaturePlot(wta.leadSNP,
                    features = "WTA1", label = TRUE, repel = TRUE) +
  ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
dev.off()



# violin plot
load(file=file.path(pathdata,"R","ScoreWTADICE_TPM1.RData"))
score.leadSNP <- score.WTA %>% select(Prediction_Healthy,WTA1)
score.random <- score.WTA %>% select(Prediction_Healthy,starts_with("WTA")) %>% select(-WTA1)

cell.levels <- score.leadSNP %>% group_by(Prediction_Healthy) %>% summarise(m=median(WTA1)) %>% arrange(desc(m)) %>% pull(Prediction_Healthy)

random.cell <- score.random %>%
  pivot_longer(-c(Prediction_Healthy),names_to="random",values_to="score") %>% group_by(Prediction_Healthy) %>%
  summarise(m=mean(score))

postscript(file="output/14.scRNAimmu/WTA_violinDICE_TPM1.ps",paper="special",width=40,height=20,onefile=FALSE,horizontal=FALSE)
score.leadSNP %>% mutate(Prediction_Healthy=factor(Prediction_Healthy,levels=cell.levels)) %>% ggplot(aes(x=Prediction_Healthy,y=WTA1)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1,outlier.alpha=0) +
  geom_jitter(shape=16,color="grey",position=position_jitter(0.1),size=0.5) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 60, hjust=1),text=element_text(size=26)) +
  xlab("") +
  ylab("Score") +
  geom_point(color="red",data=random.cell,aes(x=Prediction_Healthy,y=m))
dev.off()

# One sample t-test (mean expression vs mean random (100))
test.res <- score.WTA %>% select(Prediction_Healthy,starts_with("WTA")) %>% pivot_longer(-c(Prediction_Healthy),names_to = "random",values_to = "score") %>% 
  group_by(Prediction_Healthy,random) %>% summarise(avg=mean(score)) %>% pivot_wider(names_from = "random",values_from = "avg") %>% 
  nest(random = -c(Prediction_Healthy,WTA1)) %>% mutate(ttest = purrr::map2(random,WTA1,~ t.test(t(.x)[,1],mu=.y))) %>% 
  mutate(pval=purrr::map(ttest,~ .x$p.value)) %>% select(Prediction_Healthy,pval) %>% unnest(pval) %>% arrange(pval)

test.res.padj <- data.frame(test.res,padj=p.adjust(test.res$pval,method="fdr"))
# glimpse(test.res.padj)
write.table(test.res.padj,file="output/14.scRNAimmu/testRandom.csv",sep=",",row.names = FALSE)

