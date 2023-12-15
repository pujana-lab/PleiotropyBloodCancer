# Search for genes associated to pleiotropic SNPs

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(data.table) # ‘1.14.8’
library(biomaRt) # ‘2.56.1’

# read pleiotropic SNP list (pleio_loci.tsv)
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
# filter to conjfdr < 0.05
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# gene associated to SNPs
ensemblGRCh37.SNPs <- useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snptogene.ensembl <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","ensembl_gene_name","ensembl_type","refsnp_source"), 
                           filters = c("snp_filter","chr_name"), values = list(unique(leadSNP005 %>% pull(snpid)),unique(leadSNP005 %>% pull(chrnum))), 
                           mart=ensemblGRCh37.SNPs)
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
genes.ensembl.hgnc <- getBM(filters= "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"), 
                            values=unique(snptogene.ensembl$ensembl_gene_name),mart=ensemblGRCh37.genes)
snptogene.hgnc <- snptogene.ensembl %>% left_join(genes.ensembl.hgnc,by=c("ensembl_gene_name" = "ensembl_gene_id")) %>% 
  select(refsnp_id,hgnc_symbol,gene_biotype) %>% filter(!is.na(hgnc_symbol)) %>% filter(hgnc_symbol != "")
snptogene.hgnc <- snptogene.hgnc[!duplicated(snptogene.hgnc),]
leadSNP005.gene <- leadSNP005 %>% left_join(snptogene.hgnc,by = c("snpid" = "refsnp_id")) %>% rename(gene=hgnc_symbol) %>% mutate(gene_type = gene_biotype)
save(leadSNP005.gene,file="RData/leadSNP005Gene.RData")