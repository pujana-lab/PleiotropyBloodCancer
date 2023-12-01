# Histogram showing the distribution of identified RNA repeat elements across the pleiotropic loci

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(data.table) # ‘1.14.8’
library(biomaRt) # ‘2.56.1’

# Unique pleitropic SNPs
load(file=file.path("RData/leadSNP005.RData"))
uniqSNP <- leadSNP005 %>% distinct(snpid,chrnum,chrpos)

# Download RepeatMasker
# https://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz
rmasker <- fread(file="Data/hg19.fa.out",fill=TRUE,skip=3,header=FALSE)
rmasker.coord.type <- rmasker[,paste0("V",c(5,6,7,11))]
colnames(rmasker.coord.type) <- c("chrom","pos.init","pos.end","type")

# Select rmasker (chr1...chr22) y calculate position at the midpoint (pos.init, pos.end)
rmasker.pos.type <- rmasker.coord.type %>% filter(chrom %in% paste0("chr",seq(1,22))) %>% mutate(chrnum=as.integer(gsub("chr","",chrom))) %>%
  mutate(pos=round((pos.init+pos.end)/2)) %>% dplyr::select(chrnum,pos,type)

# Number of pleiotropic SNPs within a window +/- 50kb
overl <- NULL
for (i in seq(1,819)) {
  # print(i)
  d <- rmasker.pos.type %>% inner_join(uniqSNP[seq((i-1)*5+1,i*5),],by=c("chrnum"),relationship = "many-to-many") %>% 
    filter((chrpos > (pos-50000)) & (chrpos < (pos+50000)))
  overl <- rbind(overl,d)
}

# total number of annotated YRNA BiomaRt GRCh37
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
biomart37.misc <- getBM(filters="biotype",values=c("misc_RNA"),
                        attributes = c("ensembl_gene_id","external_gene_name","chromosome_name",
                                       "start_position","end_position","gene_biotype","strand"), 
                        mart=ensemblGRCh37.genes)
biomart37.yrnarny <- biomart37.misc %>% filter((external_gene_name=="Y_RNA" | grepl("RNY",external_gene_name)) & chromosome_name %in% seq(1,22))
totalYRNA <- dim(biomart37.yrnarny)[1]

# number YRNA next to pleiotropic SNP
yrna.snp <- biomart37.yrnarny %>% mutate(chromosome_name=as.integer(chromosome_name)) %>% 
  left_join(leadSNP005 %>% rename(chr=chrnum,pos=chrpos),by=c("chromosome_name"="chr"),relationship="many-to-many") %>% 
  filter(pos>=(start_position-50000) & pos<=(end_position+50000)) %>% 
  mutate(dist_start = pos-start_position) %>% mutate(dist_end = end_position-pos) %>% mutate(dist = case_when(
    strand == 1 ~ dist_start,
    strand == -1 ~ dist_end
  )) %>% filter(!(chromosome_name == 6 & between(pos,29500000,33500000))) 
nYRNA.pleio <- length(unique(yrna.snp %>% pull(ensembl_gene_id)))

# stats by 7x selected subgroups: proportion of feature with a pleiotropic SNP 
overl.subgrp <- overl %>% filter(type %in% c("RNA","rRNA","scRNA","SINE/tRNA-RTE","snRNA","srpRNA","tRNA")) %>% distinct(type,chrnum,pos)

rmasker.count <- rmasker.coord.type %>% filter(chrom %in% paste0("chr",seq(1,22))) %>% count(type) %>% 
  filter(type %in% c("RNA","rRNA","scRNA","SINE/tRNA-RTE","snRNA","srpRNA","tRNA"))

subgrp.count <- overl.subgrp %>% count(type) %>% left_join(rmasker.count %>% rename(total=n),by=c("type")) %>%
  transmute(type, YES=n, NO=total-n,total) %>% mutate(prop=YES/total)

subgrp.prop <- rbind(subgrp.count,data.frame(type="YRNA",YES=nYRNA.pleio,NO=totalYRNA-nYRNA.pleio,total=totalYRNA,prop=nYRNA.pleio/totalYRNA))
subgrp.levels <- subgrp.prop %>% arrange(desc(prop)) %>% pull(type)

postscript(file="Output/Figure5d.ps")
subgrp.prop %>% mutate(type=factor(type,levels=subgrp.levels)) %>%
  ggplot(aes(x=type,y=prop)) +
  geom_bar(stat="identity") + ylab("proportion of feature with SNP pleio") + xlab("") +
  theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1),text=element_text(size=18),legend.title=element_blank())
dev.off()
