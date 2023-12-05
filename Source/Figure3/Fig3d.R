# Histograms showing the percentage of pleiotropic variants identified as eQTL in whole blood (left panel) or 
# immortalized lymphocytes (right panel) of the corresponding candidate genes

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(biomaRt) # ‘2.56.1’
library(rtracklayer) # ‘1.60.1’
# library(GenomicRanges) # ‘1.52.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’
library(LDlinkR) # ‘1.3.0’
# PLINK v1.90b6.26 64-bit (2 Apr 2022)  www.cog-genomics.org/plink/1.9/

## Pleiotropic SNPs nearby (+/-50kb) target genes (ATM, CHRHR1, HLAs, APOBEC3, ATM, TERT) ####
# read pleiotropic SNPs
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  dplyr::rename(cancer_trait=trait1,
                blood_trait=trait2,
                pval_neoplasm=pval_trait1,
                pval_immunological=pval_trait2,
                zscore_neoplasm=zscore_trait1,
                zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
# regions (TERT, HLA-B, HLA-DPB1, HLA-DQA2, HLA-DQB1, CRHR1, ATM)
genes <- c("TERT","HLA-B","HLA-DPB1","HLA-DQA2","HLA-DQB1","CRHR1","ATM")
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
genes.info <- getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"), 
                             values=genes,mart=ensemblGRCh37.genes)
genes.info <- genes.info %>% filter(chromosome_name %in% seq(1,22)) %>% mutate(chromosome_name=as.integer(chromosome_name))
genes.info <- genes.info %>% mutate(startSNP=start_position-50000,endSNP=end_position+50000)
# Chromosome 22: 39,348,746-39,359,188 (APOBEC3A) Chromosome 22: 39,493,229-39,500,072 (APOBEC3H)
genes.info <- rbind(genes.info,
                    data.frame(hgnc_symbol="APOBEC3",chromosome_name=22,start_position=39348746,end_position=39500072,startSNP=39348746-50000,endSNP=39500072+50000))# SNPs nearby genes
genes.SNP <- leadSNP005 %>% 
  inner_join(genes.info %>% dplyr::select(hgnc_symbol,chromosome_name,startSNP,endSNP),by=c("chrnum"="chromosome_name"),relationship="many-to-many") %>% 
  filter(chrpos>startSNP & chrpos<endSNP) %>% dplyr::select(snpid,chrnum,chrpos,ends_with("trait"),starts_with("zscore"),hgnc_symbol)
# convert SNP to GRCh38
chainhg19toHg38 <- import.chain("Data/hg19ToHg38.over.chain")
df <- genes.SNP  %>% mutate(chrnum=paste0("chr",chrnum))
grhg19 <- GenomicRanges::GRanges(names=df %>% pull(snpid),seqnames=df %>% pull(chrnum),
                                 ranges=IRanges::IRanges(start=df %>% pull(chrpos),
                                                         end=df %>% pull(chrpos)))
grHg38 <- liftOver(grhg19, chainhg19toHg38)
SNPsGRCh38 <- as.data.frame(grHg38) %>% distinct(seqnames,start,names)
genes.SNP.GRCh38 <- genes.SNP %>% left_join(SNPsGRCh38,by=c("snpid"="names")) %>% 
  mutate(variant_id=paste(seqnames,start,sep="_")) %>% dplyr::select(-seqnames,-start,-chrnum,-chrpos) %>% distinct(hgnc_symbol,variant_id) %>%
  dplyr::rename(target=hgnc_symbol)
# assign HGNC symbol and ENSGID to target genes
gene.target <- data.frame(target=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",rep("APOBEC3",7)),
                          gene=c("TERT","ATM","HLA-B","CRHR1","HLA-DQB1","HLA-DQA2","HLA-DPB1",paste0("APOBEC3",c("A","B","C","D","F","G","H"))))
ensemblGRCh38.genes <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes.name <- getBM(filters= "hgnc_symbol", attributes = c("hgnc_symbol","ensembl_gene_id","chromosome_name"), 
                             values=gene.target %>% pull(gene),mart=ensemblGRCh38.genes) %>% filter(chromosome_name %in% seq(1,22))
gene.target <- gene.target %>% left_join(genes.name %>% dplyr::select(-chromosome_name),by=c("gene"="hgnc_symbol"))
genes.SNP.GRCh38 <- genes.SNP.GRCh38 %>% left_join(gene.target,by=c("target"),relationship="many-to-many") %>% dplyr::rename(gene_id=ensembl_gene_id)

## Build random SNP set ####
# 1. Prepare 1000G data for PLINK
# Download 1000G data 
# https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# Prepare 1000G data for PLINK
#     1000GEURplink.sh
# 2. LD matrix for SNPs nearby random selected genes
# read 1000 random selected genes
random.gene <- as.data.frame(read_excel(path="Data/Fig3d/1000randomgenes_DICE_TPM-higher1.xlsx",
                                        sheet="1000randomgenes_DICE_TPM-higher",skip=0,
                                        col_names=c("gene"),col_types=c("text"))) %>% pull(gene)
# genes coordinates
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
random.gene.ensembl <- getBM(filters= "ensembl_gene_id", attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position"), 
                             values=random.gene,mart=ensemblGRCh37.genes)
# SNPs nearby genes (+/-100kb)
for (i in seq(1,nrow(random.gene.ensembl))) {
  if (random.gene.ensembl[i,"chromosome_name"] == "X") {
    chr <- c("23")
  } else {
    chr <- random.gene.ensembl[i,"chromosome_name"]
  }
  plinkfile <- file.path("Data/1000Gplink",paste0("EUR_phase3_chr",chr))
  plinkrun <- paste0("plink --bfile ",plinkfile," --chr ",chr,
                     " --from-bp ",random.gene.ensembl[i,"start_position"]-100000,
                     " --to-bp ",random.gene.ensembl[i,"end_position"]+100000," --make-bed",
                     " --out Output/Fig3d/SNPsgene",i)
  system(plinkrun)
}
# 3. Linkage disequilibrium matrix
SNPs.files <- list.files(path="Output/Fig3d",pattern="*.bim")
ldSNPs <- list()
count <- 0
for (file in SNPs.files) {
  count <- count+1
  print(count)
  if (file.info(path.file("Output/Fig3d",file))$size != 0) {
    SNPs.sel <- fread(file.path("Output/Fig3d",file),
                      col.names=c("chrNum","rsID","x","chrPos","A1","A2")) %>% 
      mutate(SNPid=paste("chr",chrNum,":",chrPos,sep="")) %>% pull(SNPid)
    if (length(SNPs.sel)>1000) {
      SNPs.sel <- sample(SNPs.sel,1000)
    }
    ldSNPs1000[[file]] <- LDlinkR::LDmatrix(SNPs.sel,pop="EUR",r2d="r2",token="e0b6745d1a68")
  }
}
# 4. Assign each gene to a group (0.1,0.25,0.5,0.75,0.90)
# min and max for mean LD are use to fix possible groups for each gene
r2ange <- NULL
cuts <- matrix(c(0.000,0.1,0.89,0.99),ncol=2,byrow=TRUE)
for (gene in names(ldSNPs1000)) {
  # transform and clean LD matrix
  ld <- ldSNPs1000[[gene]]
  rownames(ld) <- ld$RS_number
  ldmatrix <- as.matrix(ld %>% dplyr::select(-RS_number))
  SNPnotNA <- apply(ldmatrix,2,function(x) any(!(is.na(x))))
  ldmatrix <- ldmatrix[SNPnotNA,SNPnotNA]
  vmin <- 1
  vmax <- 0
  # calculate aprox min or max value that can be obtained by randomly selecting SNP
  for (i in c(1,2)) {
    # cut-off to select most properly SNP to get min (0.0,0.1) or max value (0.89,0.99)
    cutmin <- cuts[i,1]
    cutmax <- cuts[i,2]
    nSNPs <- 30
    ldlog <- ldmatrix > cutmin & ldmatrix < cutmax
    SNPsprob <- colSums(ldlog) + 0.1
    SNPsel <- colnames(ldlog)
    # only for matrix with more than 30 SNPs
    if (length(SNPsel) >= nSNPs) {
      count <- 0
      # build random SNP set (30) and calculate min and max value
      while (count < 2000) {
        count <- count+1
        SNPran <- sample(SNPsel,nSNPs,prob=SNPsprob)
        n <- length(SNPran)
        res <- (sum(ldmatrix[SNPran,SNPran])-n)/(n*n-n)
        # keep max and min value that has been calculated 
        if(res<vmin) vmin <- res
        if(res>vmax) vmax <- res
      }
    }
  }
  r2ange <- rbind(r2ange,c(vmin,vmax))
}
colnames(r2ange) <- c("min","max")
# min and max values for each matrix are use to discard groups
ldSNPs.subgrps <- data.frame(gene=names(ldSNPs1000),r2ange) %>% mutate(G010=case_when(
  min < 0.15 ~ TRUE,
  TRUE ~ FALSE,
)) %>% mutate(G025=case_when(
  min < 0.30 & max > 0.20 ~ TRUE,
  TRUE ~ FALSE,
)) %>% mutate(G050=case_when(
  max > 0.4 ~ TRUE,
  TRUE ~ FALSE
)) %>% mutate(G075=case_when(
  max > 0.65 ~ TRUE,
  TRUE ~ FALSE
)) %>% mutate(G090=case_when(
  max > 0.75 ~ TRUE,
  TRUE ~ FALSE
))
# 5. Assign group to each gene. 
# each group should have aproximately 200 genes, but groups 0.90 and 0.75 are difficult to build so requirement is relaxed
gene090 <- sample(ldSNPs.subgrps %>% filter(G090) %>% pull(gene),170)
ldSNPs.subgrps <- ldSNPs.subgrps %>% mutate(grp =case_when(
  gene %in% gene090 ~ "G090",
  TRUE ~ NA_character_))
gene075 <- sample(ldSNPs.subgrps %>% filter(is.na(grp) & G075) %>% pull(gene),165)
ldSNPs.subgrps <- ldSNPs.subgrps %>% mutate(grp =case_when(
  gene %in% gene075 ~ "G075",
  TRUE ~ grp))
gene050 <- sample(ldSNPs.subgrps %>% filter(is.na(grp) & G050) %>% pull(gene),200)
ldSNPs.subgrps <- ldSNPs.subgrps %>% mutate(grp =case_when(
  gene %in% gene050 ~ "G050",
  TRUE ~ grp))
gene025 <- sample(ldSNPs.subgrps %>% filter(is.na(grp) & G025) %>% pull(gene),200)
ldSNPs.subgrps <- ldSNPs.subgrps %>% mutate(grp =case_when(
  gene %in% gene025 ~ "G025",
  TRUE ~ grp))
gene010 <- sample(ldSNPs.subgrps %>% filter(is.na(grp) & G010) %>% pull(gene),200)
ldSNPs.subgrps <- ldSNPs.subgrps %>% mutate(grp =case_when(
  gene %in% gene010 ~ "G010",
  TRUE ~ grp))
# 6. Create set of SNPs for each gene that meet requirements
SNPfinal <- list()
r2med <- NULL
# for example, a gene in group 0.50, random selected SNPs must not have a mean LD lower than 0.35 or higher than 0.65
set.param <- matrix(c(0.1,0.05,0.15,0.000,0.1,
                      0.25,0.20,0.30,0.89,0.99,
                      0.50,0.35,0.65,0.89,0.99,
                      0.75,0.65,0.8,0.89,0.99,
                      0.90,0.75,1.0,0.89,0.99),ncol=5,byrow = TRUE)
colnames(set.param) <- c("med","min","max","cutmin","cutmax")
rownames(set.param) <- c("G010","G025","G050","G075","G090")
for (g in ldSNPs.subgrps %>% filter(!(is.na(grp))) %>% pull(gene)) {
  # transform and clean matrix
  ld <- ldSNPs1000[[g]]
  rownames(ld) <- ld$RS_number
  ldmatrix <- as.matrix(ld %>% dplyr::select(-RS_number))
  SNPnotNA <- apply(ldmatrix,2,function(x) any(!(is.na(x))))
  ldmatrix <- ldmatrix[SNPnotNA,SNPnotNA]
  # get parameters according to the assigned group
  grp <- ldSNPs.subgrps %>% filter(gene==g) %>% pull(grp)
  param <- set.param[grp,]
  limitmin <- param[2] # Requirement 1 to accept set of SNPs: mean LD for selected SNPs lower than limitmin
  limitmax <- param[3] # Requirement 2 to accept set of SNPs: mean LD for selected SNPs higher than limitmax
  cutmin <- param[4] # Not all SNPS are going to be use to generate random SNP set
  cutmax <- param[5] # cutmin and cutmax are used to determine which SNP are available
  # start with a set of 100 random selected SNPs
  nSNPs <- 100
  
  res <- 0
  vmin <- 1
  vmax <- 0
  print(c(g,param[1],lubridate::date()))
  while (((res < limitmin) | (res > limitmax))) {
    ldlog <- ldmatrix > cutmin & ldmatrix < cutmax
    SNPsprob <- colSums(ldlog) + 0.1
    SNPsel <- colnames(ldlog)
    if (length(SNPsel) >= nSNPs) {
      count <- 0
      while (((res < limitmin) | (res > limitmax)) & count < 5000) {        
        count <- count+1
        SNPran <- sample(SNPsel,nSNPs,prob=SNPsprob)
        n <- length(SNPran)
        res <- (sum(ldmatrix[SNPran,SNPran])-n)/(n*n-n)
        if(res<vmin) vmin <- res
        if(res>vmax) vmax <- res
      }
    }
    # after several loops, if requirements have not been met, number of random selected SNPs is reduced
    nSNPs <- nSNPs-5
  }
  r2med <- c(r2med,res,nSNPs+5)
  SNPfinal[[g]] <- SNPran 
}
SNPs.grp.gene <- data.frame(gene=names(unlist(SNPfinal)),snpid=unlist(SNPfinal)) %>% mutate(gene=substr(gene,1,15)) %>% 
  left_join(ldSNPs.subgrps %>% select(gene,grp) %>% mutate(gene=substr(gene,1,15)),by=c("gene"))

## Proportions Wblood ####
# read EQTLs WBlood v1 gene_id, slope, chr_pos (GRCh38/hg38)
filename <- c("Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
# https://gtexportal.org/home/datasets
eQTLs.wblood <- fread(file.path("Data",filename)) %>% 
  separate(variant_id,c("chr","pos","A1","A2","b")) %>% mutate(variant_id=paste(chr,pos,sep="_")) %>%
  mutate(gene_id=substr(gene_id,1,15)) %>% distinct(gene_id,slope,variant_id)
# proportion of pleiotropic SNPs associated to region of interest are annotated as eQTL
SNP.prop <- genes.SNP.GRCh38 %>% left_join(eQTLs.wblood,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(target,eQTL) 
# proportion of random SNPs associated annotated as eQTL
SNP.grp.prop <- SNPs.grp.gene %>% left_join(eQTLs.wblood,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(eQTL) %>% mutate(target="Random")
# plot WBlood
SNP.prop <- rbind(SNP.prop,SNP.grp.prop) %>% mutate(eQTL=factor(eQTL,levels=c("NO","Slope < 0","Slope > 0"))) 
totalSNP <- SNP.prop %>% group_by(target) %>% summarise(total=sum(n)) 
SNP.prop <- SNP.prop %>% left_join(totalSNP,by=c("target")) %>% mutate(porc=n/total)
gene.level <- c(SNP.prop %>% filter(eQTL!="NO" & target!="Random") %>% 
                  group_by(target) %>% summarise(pTotal=sum(porc)) %>% arrange(desc(pTotal)) %>% pull(target),
                SNP.prop %>% filter(eQTL=="NO" & porc==1) %>% pull(target),"Random")
SNP.prop <- SNP.prop %>% mutate(target=factor(target,levels=gene.level))
postscript(file="Output/Fig3d/Figure3bLeft.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",linewidth=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()  

## Proportions Lymph ####
# read EQTLs Lymph v1 gene_id, slope, chr_pos (GRCh38/hg38)
filename <- c("Cells_EBV-transformed_lymphocytes.v8.signif_variant_gene_pairs.txt.gz")
# https://gtexportal.org/home/datasets
eQTLs.lymp <- fread(file.path("Data/GTEx_Analysis_v8_eQTL",filename)) %>% 
  separate(variant_id,c("chr","pos","A1","A2","b")) %>% mutate(variant_id=paste(chr,pos,sep="_")) %>%
  mutate(gene_id=substr(gene_id,1,15)) %>% distinct(gene_id,slope,variant_id)
# proportion of pleiotropic SNPs associated to region of interest are annotated as eQTL
SNP.prop <- genes.SNP.GRCh38 %>% left_join(eQTLs.lymp,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(target,eQTL) 
# proportion of random SNPs associated annotated as eQTL
SNP.grp.prop <- SNPs.grp.gene %>% left_join(eQTLs.lymp,by=c("gene_id","variant_id")) %>% mutate(eQTL=case_when(
  slope > 0 ~ "Slope > 0",
  slope < 0 ~ "Slope < 0",
  is.na(slope) ~ "NO"
)) %>% count(eQTL) %>% mutate(target="Random")
# plot Lymph
SNP.prop <- rbind(SNP.prop,SNP.grp.prop) %>% mutate(eQTL=factor(eQTL,levels=c("NO","Slope < 0","Slope > 0"))) 
totalSNP <- SNP.prop %>% group_by(target) %>% summarise(total=sum(n)) 
SNP.prop <- SNP.prop %>% left_join(totalSNP,by=c("target")) %>% mutate(porc=n/total)
gene.level <- c(SNP.prop %>% filter(eQTL=="NO" & target!="Random") %>% arrange(porc) %>% pull(target),"Random")
SNP.prop <- SNP.prop %>% mutate(target=factor(target,levels=gene.level))
postscript(file="Output/Fig3d/Figure3bRight.ps")
SNP.prop %>% ggplot(aes(x=target,y=porc,fill=eQTL)) +
  geom_bar(stat="identity",color="black",linewidth=0.1) +
  ylab("Percentage pleiotropic SNPs") + xlab("") +
  theme(axis.text.x=element_text(size=13,angle=45,hjust=1),axis.text.y=element_text(size=13),legend.text=element_text(size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title=element_text(size=14)) +
  scale_fill_manual(values=c("white","deepskyblue","coral1"),name="eQTLs")
dev.off()