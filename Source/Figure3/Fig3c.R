library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’
library(purrr) # ‘1.0.2’
# PLINK v2.00a3 64-bit (17 Feb 2020)             www.cog-genomics.org/plink/2.0/

# unique pleiotropic SNP (4093)
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
SNP.pleio <- leadSNP005 %>% distinct(snpid,chrnum,chrpos)
nSNP.pleio <- dim(SNP.pleio)[1]
# reference SNPs from plink without pleiotropic SNP
# 1. Download 1000G data 
# https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# 2. Prepare 1000G data for PLINK
#     1000GEURplink.sh
# 3. Build list SNP to generate random SNPs sets (pleiotropic SNPs were removed)
#     plink2 --bfile EUR_phase3_autosomes --recode --out EUR_phase3_autosomes_snp
#     cut -f 2 EUR_phase3_autosomes_snp.map > snps.map
#     grep -v -w -f SNPpleio.txt snps.map > snpsNOpleio.map 
SNP.non <- fread("Data/snpsNOpleio.map")
colnames(SNP.non) <- c("chrnum","snpid","chrpos")
nSNP.non <- dim(SNP.non)[1]
# chromosome length
chrlen <- c(249250621,243199373,198022430,191154276,180915260,171115067,
            159138663,146364022,141213431,135534747,135006516,133851895,	
            115169878,107349540,102531392,90354753,81195210,78077248,
            59128983,63025520,48129895,51304566)
names(chrlen) <- seq(1:22)

# number SNP in bins (1000kb, 3000kb and 5000kb)
SNP.count.chr <- NULL
for (b in c(1,3,5)) {
  bin <- b * 1000000
  cbin <- 0
  for (i in seq(1,22)) {
    timestamp()
    print(paste0(b,",",i))
    chr.bin <- data.frame(chrnum=i,nbin=seq(1,round(chrlen[i]/bin)), 
                          start=seq(0,round(chrlen[i]/bin)-1) * bin + 1, 
                          end=seq(1,round(chrlen[i]/bin)) * bin)
    SNP.non.count <- SNP.non %>% inner_join(chr.bin,by=c("chrnum"),relationship="many-to-many") %>% 
      filter(chrpos>=start & chrpos <=end) %>% count(chrnum,nbin) %>% 
      right_join(chr.bin,by=c("chrnum","nbin")) %>% mutate(n=ifelse(is.na(n),0,n)) %>% arrange(chrnum,nbin)
    SNP.pleio.count <- SNP.pleio %>% inner_join(chr.bin,by=c("chrnum"),relationship="many-to-many") %>% filter(chrpos>=start & chrpos <=end) %>% 
      count(chrnum,chrnum,nbin) %>% right_join(chr.bin,by=c("chrnum","nbin")) %>% mutate(n=ifelse(is.na(n),0,n)) %>% arrange(chrnum,nbin)
    SNP.count <- SNP.non.count %>% rename(nNon=n) %>% inner_join(SNP.pleio.count %>% rename(nPleio=n),by=c("chrnum","nbin","start","end")) %>%
      transmute(chrnum,nbin,chrbin=nbin+cbin,sbin=b,start,end,nNon,nPleio)
    SNP.count.chr <- rbind(SNP.count.chr,SNP.count)
    cbin <- cbin + dim(chr.bin)[1]
  }
}
# glimpse(SNP.count.chr)
# SNP.count.chr %>% count(sbin)
# save(SNP.count.chr,file=file.path(pathR,"testPropHotspotMarc.RData"))
#save(SNP.count.chr,file=file.path(pathR,"testPropHotspotPLINK.RData"))

# chi-squared test
chisq.res <- SNP.count.chr %>% mutate(totalNon=nSNP.non-nNon,totalPleio=nSNP.pleio-nPleio) %>% 
  nest(mat=c(nNon,nPleio,totalNon,totalPleio)) %>% mutate(table=purrr$map(mat,~matrix(unlist(.x[1,]),nrow=2))) %>%
  mutate(chiqt=purrr$map(table,~chisq.test(.x))) %>% mutate(pval=purrr$map(chiqt,~ .x$p.value)) %>% 
  select(chrnum,nbin,chrbin,sbin,start,end,pval,mat) %>% unnest(mat) %>% unnest(pval) %>% 
  mutate(OR=(nPleio*totalNon)/(nNon/totalPleio)) %>% mutate(logOR=log10(OR+1))
chisq.res <- data.frame(chisq.res,padj=p.adjust(chisq.res %>% pull(pval),method="fdr"))
# number of cancer traits with pleiotropic SNP within each bin
nCancer.chr <- NULL
for (b in c(1,3,5)) {
  bin <- b * 1000000
  cbin <- 0
  for (i in seq(1,22)) {
    timestamp()
    print(paste0(b,",",i))
    chr.bin <- data.frame(chrnum=i,nbin=seq(1,round(chrlen[i]/bin)), 
                          start=seq(0,round(chrlen[i]/bin)-1) * bin + 1, 
                          end=seq(1,round(chrlen[i]/bin)) * bin)
    nCancer <- leadSNP005 %>% inner_join(chr.bin,by=c("chrnum"),relationship="many-to-many") %>% 
      filter(chrpos>start & chrpos<end) %>% 
      distinct(chrnum,nbin,cancer_trait) %>% count(chrnum,nbin) %>%
      right_join(chr.bin,by=c("chrnum","nbin")) %>% mutate(n=ifelse(is.na(n),0,n)) %>% arrange(chrnum,nbin) %>%
      mutate(chrbin=nbin+cbin) %>% rename(nCancer=n) 
    nCancer.chr <- rbind(nCancer.chr,data.frame(nCancer,sbin=b)) %>% select(chrnum,nbin,sbin,chrbin,start,end,nCancer)
    cbin <- cbin + dim(chr.bin)[1]
  }
}
# plot selected region FDR < 0.05 and nCancer >= 3
# calculate cumulative bp to sum each chromosome and medium point for chromosome
df <- chisq.res %>% inner_join(nCancer.chr,by=c("chrnum","nbin","chrbin","sbin","start","end"))
sumchrlen <- 0
poschr <- NULL
for (i in seq(1,length(chrlen))) {
  sumchrlen <- c(sumchrlen,chrlen[i]+sumchrlen[i])
  poschr <- c(poschr,(sumchrlen[i]+sumchrlen[i+1])/2)
}
chrsum <- data.frame(chrnum=seq(1,23),sumchr=sumchrlen)
cairo_ps(file="Output/Figure3c.ps",height = 6, width = 12)
df %>% inner_join(chrsum,by=c("chrnum")) %>% filter(padj < 0.05 & nCancer >=3) %>% mutate(pos=(start+end)/2+sumchr,sbin=factor(sbin,levels=c(1,3,5))) %>%
  ggplot(aes(x=pos,y=nCancer,color=sbin)) +
  geom_point(size=1) + 
  theme_classic() +
  scale_x_continuous(label=seq(1,22),breaks=poschr,expand=c(0,0)) +
  geom_vline(xintercept = chrsum %>% pull(sumchr), color = "lightblue", linewidth=0.1) +
  xlab("")
dev.off()