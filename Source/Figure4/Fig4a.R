# Graph showing the proportion of pleiotropic variants (all cancers included) mapped in enhancersfrom immune cell types and blood (X-axis)

library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’
library(ggplot2) # ‘3.4.3’
library(data.table) # ‘1.14.8’

load(file="RData/leadSNP005.RData")

leadSNP.uniq <- leadSNP005 %>% distinct(snpid,chrnum,chrpos) %>% 
  transmute(snpid,chr=as.character(paste0("chr",chrnum)),pos=chrpos)

enh.files <- list.files("Data/HTEA")
nSNP <- list()
for (f in enh.files) {
  cell <- gsub("_differentially_expressed_enhancers.bed","",f)
# https://enhancer.binf.ku.dk/human_enhancers/presets
  cell <- gsub("CL_[0-9]*_","",cell)
  cell <- gsub("UBERON_[0-9]*_","",cell)
  enhancer <- fread(file=file.path("Data/HTEA",f),
                                skip=1,header = FALSE) %>% transmute(chr=V1,start=V2,end=V3)
  nSNP[[cell]] <- length(unique(leadSNP.uniq %>% left_join(enhancer,by=c("chr"),relationship="many-to-many") %>% filter(pos>start & pos<end) %>% pull(snpid)))
}

numSNP <- unlist(nSNP)
stats.cell <- data.frame(cell=names(numSNP),n=numSNP,total=4093)

data.test <- expand.grid(obs=stats.cell %>% filter(!(cell %in% c("brain","adipose_tissue"))) %>% pull(cell),ref=c("brain","adipose_tissue")) %>% 
  left_join(stats.cell,by=c("obs"="cell")) %>% rename(nA=n,nTA=total) %>% left_join(stats.cell,by=c("ref"="cell")) %>% rename(nB=n,nTB=total)

PropTest <- data.test %>% mutate(pHatA=nA/nTA) %>% mutate(pHatB=nB/nTB) %>% 
  mutate(pHat=(nA+nB)/(nTA+nTB)) %>% mutate(pHatnum=pHat*(1-pHat)) %>%
  mutate(SE=sqrt((pHatnum/nTA)+(pHatnum/nTB))) %>% mutate(zscore=(pHatA-pHatB)/SE) %>%
  mutate(pval=pnorm(zscore,lower.tail=FALSE)) %>% mutate(sign = case_when(
    pval < 0.05 ~ "*",
    TRUE ~ ""
  ))

PropTest <- data.frame(PropTest,padj=p.adjust(PropTest %>% pull(pval),method="BH")) %>% mutate(signadj = case_when(
  padj < 0.05 ~ "*",
  TRUE ~ ""
))

data.plot <- PropTest %>% distinct(obs,nA,nTA) %>% mutate(prop=nA/nTA)
obs.levels <- data.plot %>% arrange(desc(prop)) %>% pull(obs)
data.plot.ref <- PropTest %>% distinct(ref,nB,nTB) %>% mutate(prop=nB/nTB)
postscript(file="Output/Figure4a.ps")
data.plot %>% mutate(obs=factor(obs,levels=obs.levels)) %>% 
  ggplot(aes(x=obs,y=prop)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust=1), text=element_text(size=14)) +
  xlab("") +
  ylab("Proportion SNPs") +
  geom_hline(yintercept=data.plot.ref %>% filter(ref=="brain") %>% pull(prop),linetype="dashed",color="red") +
  geom_hline(yintercept=data.plot.ref %>% filter(ref=="adipose_tissue") %>% pull(prop),linetype="dashed",color="blue")
dev.off()
