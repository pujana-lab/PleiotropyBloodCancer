# TODO Listado de regiones para los 1000 genes
library(readxl) # ‘1.4.3’
library(biomaRt) # ‘2.56.1’
library(data.table) #
library(LDlinkR) # ‘1.3.0’
# PLINK v1.90b6.26 64-bit (2 Apr 2022)  www.cog-genomics.org/plink/1.9/

## Prepare 1000G data for PLINK
# 1. Download 1000G data 
# https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# 2. Prepare 1000G data for PLINK
#     1000GEURplink.sh
#     splitfilechr.sh

## LD matrix for SNPs nearby random selected genes
# read 1000 random selected genes
random.gene <- as.data.frame(read_excel(path="Data/Fig3d/1000randomgenes_DICE_TPM-higher1.xlsx",
                                        sheet="1000randomgenes_DICE_TPM-higher",skip=0,
                                        col_names=c("gene"),col_types=c("text"))) %>% pull(gene)
# genes coordinates
ensemblGRCh37.genes <- useEnsembl(biomart = "genes", version="GRCh37", dataset = "hsapiens_gene_ensembl")
random.gene.ensembl <- getBM(filters= "ensembl_gene_id", attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position"), 
                                      values=random.gene,mart=ensemblGRCh37.genes)
# SNPs nearby genes (+/-100kb)
# plink --bfile dataExt/1000Gplink/EUR_phase3_chr9 --chr 9 --from-bp 4679559 --to-bp 4708398 --make-bed --out GENE
pathExt <- c("/media/mapardo/SeagateBasic/procure/dataExt")
for (i in seq(1,nrow(random.gene.ensembl))) {
  if (random.gene.ensembl[i,"chromosome_name"] == "X") {
    chr <- c("23")
  } else {
    chr <- random.gene.ensembl[i,"chromosome_name"]
  }
  #plinkfile <- file.path("Data/1000Gplink",paste0("EUR_phase3_chr",chr))
  plinkfile <- file.path(pathExt,"1000Gplink",paste0("EUR_phase3_chr",chr))
  plinkrun <- paste0("plink --bfile ",plinkfile," --chr ",chr,
                     " --from-bp ",random.gene.ensembl[i,"start_position"]-100000,
                     " --to-bp ",random.gene.ensembl[i,"end_position"]+100000," --make-bed",
                     " --out Output/Fig3d/SNPsgene",i)
  system(plinkrun)
}
# linkage disequilibrium matrix
SNPs.files <- list.files(path="Output/Fig3d",pattern="*.bim")
ldSNPs <- list()
count <- 0
# file <- SNPs.files[10]
for (file in SNPs.files) {
  count <- count+1
  print(count)
  print(lubridate::date())
  if (file.info(path.file("Output/Fig3d",file))$size != 0) {
    SNPs.sel <- fread(file.path("Output/Fig3d",file),
                                  col.names=c("chrNum","rsID","x","chrPos","A1","A2")) %>% 
      mutate(SNPid=paste("chr",chrNum,":",chrPos,sep="")) %>% pull(SNPid)
    if (length(SNPs.sel)>1000) {
      SNPs.sel <- sample(SNPs.sel,1000)
    }
    ldSNPs[[file]] <- LDlinkR::LDmatrix(SNPs.sel,pop="EUR",r2d="r2",token="e0b6745d1a68")
  }
}

names(ldSNPs)
# asignación de grupos
# Calcular los subgrupos en los que pueden participar 1000Genes (0.1,0.25,0.5,0.75,0.90)

r2ange <- NULL
cuts <- matrix(c(0.000,0.1,0.89,0.99),ncol=2,byrow=TRUE)
# for (gene in names(ldSNPs1000)) {
for (gene in names(ldSNPs)) {
  # transform and clean LD matrix
  # ld <- ldSNPs1000[[gene]]
  ld <- ldSNPs[[gene]]
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
  # print(lubridate::date())
  # print(c(gene,vmin,vmax))
  r2ange <- rbind(r2ange,c(vmin,vmax))
}
colnames(r2ange) <- c("min","max")
# min and max values for each matrix are use to discard groups
#ldSNPs.subgrps <- data.frame(gene=names(ldSNPs1000),r2ange) %>% mutate(G010=case_when(
ldSNPs.subgrps <- data.frame(gene=names(ldSNPs),r2ange) %>% mutate(G010=case_when(
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

# assign group to each gene. Each group should have 200 genes aprox. 
# Groups 0.90 and 0.75 are difficult to build so requeriment is relaxed
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
#save(ldSNPs.subgrps,file=file.path(path,"data","pleiotropy","randomSNPs","LinkDeseq","ldSNPs1000Group.RData"))

ldSNPs.subgrps <- data.frame(ldSNPs.subgrps,grp=c("G025","G010","G010","G050","G075","G010","G010","G010","G010")) #####

# create set of SNPs for each gene that meet requirements
# for example, a gene in group 0.50, selected SNPs must not have a mean LD lower than 0.35 or higher than 0.65
SNPfinal <- list()
r2med <- NULL
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
  g <- c("SNPsgene109.bim") ####
  ld <- ldSNPs[[g]] ####
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
      #while (((res < limitmin) | (res > limitmax)) & count < 50000) {
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
  print(c(g,nSNPs+5,res,lubridate::date()))
  r2med <- c(r2med,res,nSNPs+5)
  SNPfinal[[g]] <- SNPran 
}
#save(SNPfinal,file=file.path(path,"data","pleiotropy","randomSNPs","LinkDeseq","SNPSelFinal.RData"))

# numSNPs <- unlist(lapply(SNPfinal,length)) 
# names(numSNPs)[sum(numSNPs < 25)]
# [1] "ENSG00000007080SEL.vcf" Solo uno con menos de 25 SNPs

# # Listado de SNPS para grupo (x5)
# SNPs.grp <- list()
# for (set in c("G010","G025","G050","G075","G090")) {
#   SNPs.grp[[set]] <-NULL
#   for (g in ldSNPs.subgrps %>% filter(grp==set) %>% pull(gene)) {
#     SNPs.grp[[set]] <- c(SNPs.grp[[set]],SNPfinal[[g]])
#   }
#   SNPs.grp[[set]] <- unique(SNPs.grp[[set]])
# }
# unlist(lapply(SNPs.grp,length))
# # G010  G025  G050  G075  G090 
# #19428 16759 16749  9433 11686 
# save(SNPs.grp,file=file.path(path,"data","pleiotropy","randomSNPs","LinkDeseq","SNPsGroups.RData"))


# AQUI SE CONSIGUE SNPs.grp.gene
pathaux <- c("/home/mapardo/procure/pleiotStudy")
load(file=file.path(pathaux,"data","randSNPs1000Genes","LinkDeseq","SNPSelFinal.RData"))
load(file=file.path(pathaux,"data","randSNPs1000Genes","LinkDeseq","ldSNPs1000Group.RData"))
SNPs.grp.gene <- data.frame(gene=names(unlist(SNPfinal)),snpid=unlist(SNPfinal)) %>% mutate(gene=substr(gene,1,15)) %>% 
  left_join(ldSNPs.subgrps %>% select(gene,grp) %>% mutate(gene=substr(gene,1,15)),by=c("gene"))

