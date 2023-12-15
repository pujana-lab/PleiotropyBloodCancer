# Forest plot showing the OR and 95% CI of the overlap between the pleiotropic gene set and 
# hematopoiesis gene modules, depicted by the corresponding master regulators (Y-axis)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(readxl) # ‘1.4.3’
library(purrr) # ‘1.0.2’
library(questionr) # ‘0.7.8’
library(biomaRt) # ‘2.56.1’

# association genes/pleiotropic SNPs (gene2SNP.R)
load(file="RData/leadSNP005Gene.RData")
# hematopoiesis data
d <- read_xlsx("Data/GeneModules_Hematopiesis_Verten_Ind1+2_sets_071923.xlsx",sheet="Merged1-2")
modules <- lapply(lapply(d, as.vector),function(x) {x[!is.na(x)]})

# ensembl gene info from biomaRt ---> GRCh38 <--- ---> only canonical chromosomes <---
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes.biomart.all <- getBM(attributes=c('hgnc_symbol','hgnc_id','gene_biotype',
                                                'ensembl_gene_id','ensembl_gene_id_version','chromosome_name'),mart = ensembl) 
genes.biomart <- genes.biomart.all %>% filter(chromosome_name %in% c(seq(1,22),"MT","X","Y"))
pc <- genes.biomart %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
modules.clean <- lapply(modules,function(x){x[x %in% pc]})
pleio <- unique(leadSNP005.gene %>% pull(gene))
pleio.clean <- pleio[pleio %in% pc]
# contingency table
noverlap <- unlist(lapply(modules.clean,function(x){sum(x%in%pleio.clean)}))
nmodules <- unlist(lapply(modules,length))
npleio <-length(pleio.clean)
npc <-length(pc)
yn <- nmodules - noverlap
ny <- npleio - noverlap
nn <- npc - yn - ny - noverlap
df <- data.frame(module=names(modules),yy=noverlap,ny=ny,yn=yn,nn=nn)
# odd ratio
res <- df %>% nest(d=-c(module)) %>% mutate(ORres=map(d,~ odds.ratio(matrix(unlist(as.vector(.x)),nrow=2),level=0.95))) %>%
  mutate(OR=map(ORres,~ .x$OR)) %>% unnest(OR) %>% mutate(LL=map(ORres,~ .x$`2.5 %`)) %>% unnest(LL) %>%
  mutate(UL=map(ORres,~ .x$`97.5 %`)) %>% unnest(UL) %>% dplyr::select(-d, -ORres)
module.levels <- res %>% arrange(OR) %>% pull(module)
# plot
cairo_ps(file="Output/Figure4b.ps",width=10,height=10,onefile=FALSE,fallback_resolution = 600)
res %>% filter(OR != 0) %>% mutate(module=factor(module,levels=module.levels)) %>% 
  ggplot(aes(y = module, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", linewidth = 1, alpha = 0.5) +
  coord_trans(x="log10") +
  xlim(0.001,25) +
  scale_x_continuous(breaks=c(0.1,10)) +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))
dev.off()
