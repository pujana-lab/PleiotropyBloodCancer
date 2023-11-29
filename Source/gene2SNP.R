# search for gene associated to SNPs

library(biomaRt) # [1] ‘2.56.1’ 

load(file="RData/leadSNP005.RData")

leadSNP005.rsid <- leadSNP005 %>% distinct(snpid,chrnum)

ensemblGRCh37.SNPs <- useEnsembl(biomart = "snps", version="GRCh37", dataset = "hsapiens_snp")
snptogene.ensembl <- NULL
for (i in seq(1,21)) {
  print(i)
  s <- 200*(i-1) + 1
  e <- 200*(i)
  l <- dim(leadSNP005.rsid)[1]
  if (e < l) {
    r <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","ensembl_gene_name","ensembl_type","refsnp_source"),
               filters = c("snp_filter","chr_name"), 
               values = list(leadSNP005.rsid[seq(s,e),] %>% pull(snpid),leadSNP005.rsid[seq(s,e),] %>% pull(chrnum)),
               mart=ensemblGRCh37.SNPs)
  } else {
    r <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","ensembl_gene_name","ensembl_type","refsnp_source"),
               filters = c("snp_filter","chr_name"), 
               values = list(leadSNP005.rsid[seq(s,l),] %>% pull(snpid),leadSNP005.rsid[seq(s,l),] %>% pull(chrnum)),
               mart=ensemblGRCh37.SNPs)
  }
  snptogene.ensembl <- rbind(snptogene.ensembl,r)
}
save(snptogene.ensembl,file="RData/SNPtoGene.RData")

# Search for HGNC symbol gene (Avoid duplications for Entrez code)
ensemblGRCh37.genes <- biomaRt::useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
genes.ensembl.hgnc <- biomaRt::getBM(filters= "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"), 
                                     values=unique(snptogene.ensembl$ensembl_gene_name),mart=ensemblGRCh37.genes)
save(genes.ensembl.hgnc,file="RData/GenesEnsemblHGNC.RData")

snptogene.hgnc <- snptogene.ensembl %>% left_join(genes.ensembl.hgnc,by=c("ensembl_gene_name" = "ensembl_gene_id"),relationship="many-to-many") %>% 
  dplyr::select(refsnp_id,hgnc_symbol,gene_biotype) %>% filter(!is.na(hgnc_symbol)) %>% filter(hgnc_symbol != "")
snptogene.hgnc <- snptogene.hgnc[!duplicated(snptogene.hgnc),]

# Add gene info to leadSNP file
leadSNP005.gene <- leadSNP005 %>% left_join(snptogene.hgnc,by = c("snpid" = "refsnp_id"),relationship="many-to-many") %>% rename(gene=hgnc_symbol) %>% mutate(gene_type = gene_biotype)
write.table(leadSNP005.gene,file="Data/leadSNPgene.tsv",sep="\t",row.names = FALSE)
save(leadSNP005.gene,file="RData/leadSNP005Gene.RData")