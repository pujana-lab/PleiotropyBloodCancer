# Phylogenetic analysis of RNY sequences from the human genome

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(Biostrings) # ‘2.68.1’
library(msa) # ‘1.32.0’
library(ape) # ‘5.7.1’
library(ggtree) # ‘3.8.2’
library(biomaRt) # ‘2.56.1’

## Download gene sequences from biomart (YRNA816.fa)
# RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
# write list RNY ensembl ID
write.table(biomart37.yrnarny %>% pull(ensembl_gene_id),file="Data/FigS6/YRNA37ensgID816.txt",row.names = FALSE,quote = FALSE,col.names = FALSE)
# search sequences in biomart GRCh37 database
# YRNA37ensgID816.txt ---> http://grch37.ensembl.org/biomart ---> YRNA816.fa
# YRNA816 and YRNA37ensgID816.txt files available at Data directory

## Build phylogenetic tree
# multiple sample alignment (ClustalW) + distance matrix (APE)
# readDNAStringSet() was used to read FASTA format file 
seq <- readDNAStringSet("Data/FigS6/YRNA816.fa")
# all samples were aligned to the same length by ClustalW algorithm
multAlign <- msaClustalW(seq)
# store multiple sequence alignments as a DNAbin object 
nbin <- as.DNAbin(multAlign)
# The msaplot() command provided by ggtree package and ggplot2 package was used to 
# demonstrate the aligned sequences with the phylogenetic tree [11]. 
YRNA.labels <- data.frame(YRNA=labels(nbin),label=paste0("YRNA_",seq(1,length(labels(nbin)))))
rownames(nbin) <- paste0("YRNA_",seq(1,length(labels(nbin))))
an <- as.alignment(nbin)  # converting DNAbin to alignment format
nm <- as.matrix(an)         # converting alignment to matrix
nbinmat <- as.matrix(labels(nbin)) # extraction of the sample names
dnbin <- dist.dna(nbin, model = "raw",pairwise.deletion=TRUE) # computing distance by ape package with K80 model derived by Kimura (1980)

## plot phylogenetic tree
# select RNY nearby pleiotropic SNPs to plot in red color
load(file="RData/leadSNP005.RData")
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% dplyr::rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
yrna.pleio <- data.frame(YRNA=unique(yrna.snp %>% pull(ensembl_gene_id)),SNP="Pleio")
rny <- biomart37.yrnarny %>% filter(external_gene_name != "Y_RNA") %>% dplyr::select(ensembl_gene_id,external_gene_name)
YRNA.plot <- YRNA.labels %>% left_join(yrna.pleio,by=c("YRNA")) %>% 
  mutate_if(is.character,coalesce,"NoPleio") %>% left_join(rny,by=c("YRNA"="ensembl_gene_id")) %>% mutate(external_gene_name=case_when(
    is.na(external_gene_name) ~ YRNA,
    TRUE ~ external_gene_name
  ))

# generate plot
tree<-nj(dnbin)
YRNA.pleio <- list(Pleio=which(YRNA.plot$SNP == "Pleio" & YRNA.plot$external_gene_name != "RNY4"),
                   NoPleio=which(YRNA.plot$SNP == "NoPleio" & YRNA.plot$external_gene_name != "RNY4"),
                   RNAY4=which(YRNA.plot$external_gene_name == "RNY4"))
tree$tip.label <- YRNA.plot$external_gene_name
branches <- YRNA.pleio
tree <- groupOTU(tree,branches)
ggt <- ggplot(tree) + ggtree::geom_tree() + ggtree::theme_tree2() + ggtree::geom_tiplab(aes(color=group),size=3,family="sans") +
  scale_color_manual(values=c("black","red","green"))
postscript(file="Output/FigureS6.ps",family="sans",
           paper="special",width=30,height=80,onefile=FALSE,horizontal=FALSE)
ggt
dev.off()
