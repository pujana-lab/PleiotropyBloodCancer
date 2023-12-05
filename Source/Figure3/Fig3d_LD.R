library(readxl) # ‘1.4.3’
library(data.table)
library(dplyr) # ‘1.1.3’
library(tidyr) # ‘1.3.0’


# linkage disequilibrium matrix
SNPs.files <- list.files(path="Output/Fig3d",pattern="*.bim")
ldSNPs <- list()
count <- 0
# file <- SNPs.files[seq(1,10)]
#for (file in SNPs.files) {
for (file in SNPs.files[seq(11,200)]) {
  count <- count+1
  print(count)
  print(lubridate::date())
  if (file.info(file.path("Output/Fig3d",file))$size != 0) {
    SNPs.sel <- fread(file.path("Output/Fig3d",file),
                      col.names=c("chrNum","rsID","x","chrPos","A1","A2")) %>% 
      mutate(SNPid=paste("chr",chrNum,":",chrPos,sep="")) %>% pull(SNPid)
    if (length(SNPs.sel)>1000) {
      SNPs.sel <- sample(SNPs.sel,1000)
    }
    ldSNPs[[file]] <- LDlinkR::LDmatrix(SNPs.sel,pop="EUR",r2d="r2",token="e0b6745d1a68")
  }
}
save(ldSNPs,file="Output/Fig3d/ldSNPs11_200.RData")
