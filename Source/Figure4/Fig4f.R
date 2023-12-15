# Venn diagrams showing the overlap between mouse gene orthologs linked to myeloid cell alterations (phenotypes are indicated) 
# and the pleiotropic gene set (all cancers included)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(purrr) # ‘1.0.2’
library(readxl) # ‘1.4.3’
library(data.table) # ‘1.14.8’
source("Source/ownFunc.R")

# association genes/pleiotropic SNPs (gene2SNP.R)
load(file="RData/leadSNP005Gene.RData")
# read genes link to myeloid cell alterations
pheno <- lapply(as.list(read_excel(path="Data/Myeloid-mouse-genes_phenotypes-5.xlsx", sheet="Sheet1",skip=0)),function(x) x[!is.na(x)])
# database MGI symbol MGI <-> HGNC
human.mouse.mgi <- fread("Data/kk/HOM_MouseHumanSequence.rpt",header=TRUE) %>% select("DB Class Key","Common Organism Name",Symbol)
# https://www.informatics.jax.org/downloads/reports/index.html#homology
colnames(human.mouse.mgi) <- c("DBClassKey","Name","Symbol")
human.mgi <- human.mouse.mgi %>% filter(Name == "human") %>% rename(hgncSymbolHuman=Symbol)
mouse.mgi <- human.mouse.mgi %>% filter(Name == "mouse, laboratory") %>% rename(mgiSymbolMouse=Symbol)
human.mouse.mgi <- mouse.mgi %>% inner_join(human.mgi,by=c("DBClassKey"),relationship="many-to-many") %>% select("hgncSymbolHuman","mgiSymbolMouse")

# overlap pleiotropic genes vs myeloid cell alterations genes
leadSNP005.mgi <- leadSNP005.gene %>% inner_join(human.mouse.mgi,by=c("gene"="hgncSymbolHuman"),relationship="many-to-many") %>%
  transmute(snpid,chrnum,chrpos,A1,A2,cancer_trait,blood_trait,mgi=mgiSymbolMouse)
pheno.clean <- lapply(pheno,function(x) x[x %in% human.mouse.mgi$mgi])
pheno.leadSNP005 <- lapply(pheno.clean,function(x) x[x %in% leadSNP005.mgi$mgi])
df.group <- data.frame(pheno=names(pheno.clean),total=length(unique(human.mouse.mgi$mgiSymbolMouse)),
                       group1=length(unique(leadSNP005.mgi %>% pull(mgi))),
                       group2=unlist(lapply(pheno.clean,length)),
                       overlap=unlist(lapply(pheno.leadSNP005,length))) 
# Fisher's exact test
df.cont <- df.group %>% mutate(yy=overlap,ny=group2-overlap,yn=group1-overlap,nn=total-group2-group1+overlap) %>%
  select(pheno,yy,ny,yn,nn)
df.res <- df.cont %>% pivot_longer(-c(pheno),names_to = "type",values_to = "n") %>% nest(d=-c(pheno)) %>%
  mutate(mat=map(d,~ matrix(.x %>% pull(n),2,2))) %>% mutate(ft=map(mat, ~ fisher.test(.x,alternative="greater"))) %>%
  mutate(pval=map(ft,~ .x$p.value)) %>% unnest(pval) %>% mutate(OR=map(ft,~ .x$estimate)) %>% unnest(OR) %>%
  select(pheno,pval,OR) %>% left_join(df.group,by=c("pheno")) %>% left_join(df.cont,by=c("pheno"))
# plot
data.venn <- df.res 
lead.mgi <- unique(leadSNP005.mgi$mgi)
data.venn$pleio.gene <- list(lead.mgi,lead.mgi,lead.mgi,lead.mgi,lead.mgi)
data.venn$pheno.gene <- pheno.clean
plots <- data.venn %>% nest(d=c(pval,OR,pleio.gene,pheno.gene)) %>% mutate(plots=map2(d,pheno,vennDiagramPlot)) %>% pull(plots)
cairo_ps(file="Output/Figure4f.ps",width=10,height=15,onefile=FALSE,fallback_resolution = 600)
ggpubr::ggarrange(plotlist=plots,ncol = 2, nrow = 3)
dev.off()