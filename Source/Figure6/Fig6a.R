# Density distribution of the pleiotropic SNPs identified nearby (± 50 kb) RNY TSSs

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(biomaRt) # ‘2.56.1’

# RNY annotated BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))

# Assignment of RNY to nearby SNP
load(file="RData/leadSNP005.RData")

yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 

# Density plot
yrna.snp.count <- yrna.snp %>% distinct(ensembl_gene_id,snpid) %>% count(ensembl_gene_id) %>% 
  mutate(total.yrna=length(unique(yrna.snp %>% pull(ensembl_gene_id))))

yrna.snp.dens <- yrna.snp %>% distinct(ensembl_gene_id,snpid,dist) %>% left_join(yrna.snp.count,by=c("ensembl_gene_id")) %>% mutate(s=1/(n*total.yrna))

postscript(file="Output/Figure6a.ps")
yrna.snp.dens %>%
  ggplot(aes(x=dist,weight=s)) +
  geom_density(show.legend=FALSE)+
  ylab("Density pleiotropic SNPs") +
  xlab("5' RNY gene/pseudogene") +
  scale_x_continuous(breaks = c(-60000,-40000,-20000,0,20000,40000,60000),limits = c(-65000,65000)) +
  scale_y_continuous(breaks=c(c(0,2,4,6,8,10,12)*1e-6),labels=c(0,2,4,6,8,10,12),limits=c(0,15)*1e-6) +
  theme_classic() +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18)) +
  geom_vline(xintercept=0,linetype="dotted") +
  geom_vline(xintercept=-50000) +
  geom_vline(xintercept=50000) +
  geom_hline(yintercept=0,colour="grey")
dev.off()
