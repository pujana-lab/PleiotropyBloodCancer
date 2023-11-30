# read data and save in .RData files

library(data.table) # ‘1.14.8’
library(readxl) # ‘1.4.3’

# read pleiotropic SNP list (pleio_loci.tsv). Filter to conjfdr < 0.05
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
save(leadSNP005, file="RData/leadSNP005.RData")

# metadata info
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
save(metadata, file="RData/metadata.RData")

# label for cancer traits 
cancer.label <- data.frame(LabelInit=c("Onco iCOGs GWAS meta analysis","Onco iCOGs meta analysis luminal A","Breast cancer",
                                       "Cervix cancer","Prostate cancer","Melanoma",
                                       "CIMBA BRCA1/BCAC Triple-Negative meta analysis","Non-Hodgkin's Lymphoma","BRCA1 carrier ovary cancer",
                                       "Onco iCOGs meta analysis Triple N","Onco iCOGs meta analysis luminal B HER2 Neg","Onco iCOGs meta analysis luminal B",
                                       "Colon cancer","Lymphocytic Leukemia","Oral Cavity/Pharynx cancer",
                                       "Lung cancer","BRCA1 carrier breast cancer","Endometrium cancer",
                                       "Thyroid cancer","Pancreas cancer","Onco iCOGs meta analysis HER2 enriched",
                                       "Rectum cancer","Bladder cancer","BRCA2 carrier breast cancer",
                                       "BRCA2 carrier ovary cancer","Kidney cancer","Ovary cancer",
                                       "Esophagus/Stomach cancer"),
                           LabelNew=c("BC#1","BC#1 LumA","BC#2",
                                      "Cervix","Prostate","Melanoma",
                                      "BRCA1 TNBC","NHL","BRCA1 OC",
                                      "BC#1 TNBC","BC#1 LumB/HER2-","BC#1 LumB",
                                      "Colon","Leukemia","Oropharyngeal",
                                      "Lung","BRCA1 BC","Endometrium",
                                      "Thyroid","Pancreas","BC#1 HER2+",
                                      "Rectum","Bladder","BRCA2 BC",
                                      "BRCA2 OC","Kidney","Ovary",
                                      "Gastroesophageal"))
save(cancer.label,file="RData/CancerLabel.RData")

# label for blood traits
blood.label <- data.frame(LabelInit=c("Lymphocyte count","Monocyte count","Eosinophill count","Reticulocyte count",
                                      "White blood cell (leukocyte) count","Reticulocyte percentage","High light scatter reticulocyte count","Haemoglobin concentration",
                                      "Mean corpuscular volume","High light scatter reticulocyte percentage","Platelet count","Neutrophill count",
                                      "Mean platelet (thrombocyte) volume","Neutrophill percentage","Lymphocyte percentage","Eosinophill percentage",
                                      "Monocyte percentage","Red blood cell (erythrocyte) distribution width","Mean corpuscular haemoglobin concentration","Platelet distribution width",
                                      "Haematocrit percentage","Red blood cell (erythrocyte) count","Platelet crit","Immature reticulocyte fraction",
                                      "Basophill count","Basophill percentage","Nucleated red blood cell count"),
                          LabelNew=c("LYMPH#","MONO#","EO#","RET#",
                                     "WBC#","RET%","HLSR#","HGB",
                                     "MCV","HLSR%","PLT#","NEUT#",
                                     "MPV","NEUT%","LYMPH%","EO%",
                                     "MONO%","RDW","MCHC","PDW",
                                     "HCT","RBC#","PCT","IRF",
                                     "BASO#","BASO%","NRBCD"))
save(blood.label,file="RData/BloodLabel.RData")

# Cancer levels
cancer.levels <- c("Thyroid","Rectum","Prostate","Pancreas","Oropharyngeal","NHL","Melanoma","Lung","Leukemia","Kidney",
                   "Gastroesophageal","Endometrium","Colon","Cervix","Bladder","Ovary","BRCA2 OC","BRCA1 OC","BC#2","BRCA2 BC",
                   "BRCA1 BC","BRCA1 TNBC","BC#1 TNBC","BC#1 HER2+","BC#1 LumB/HER2-","BC#1 LumB","BC#1 LumA","BC#1")
save(cancer.levels,file="RData/CancerLevels.RData")

# Blood levels
blood.levels <- c("HCT","HGB","MCHC","NRBCD","RBC#","RDW","MCV",
                  "HLSR#","HLSR%","IRF","RET#","RET%",
                  "MPV","PCT","PDW","PLT#",
                  "BASO#","BASO%",
                  "EO#","EO%", 
                  "LYMPH#","LYMPH%",  
                  "MONO#","MONO%", 
                  "NEUT#","NEUT%",  
                  "WBC#")
save(blood.levels,file="RData/BloodLevels.RData")

# Genes telomere 
telomereSNP <- read_excel("Data/13_genes_telomere-snps.xlsx", 
                          sheet="Sheet1",col_names=paste0("V",seq(1,18)))
telomereSNP005 <- telomereSNP %>% filter(V9 < 0.05)
save(telomereSNP005,file="RData/telomereGene005.RData")

# Number bases at each chr. Start y end for each chr in an acumulative way
chrlen <- c(249250621,243199373,198022430,191154276,180915260,171115067,
            159138663,146364022,141213431,135534747,135006516,133851895,	
            115169878,107349540,102531392,90354753,81195210,78077248,
            59128983,63025520,48129895,51304566)
names(chrlen) <- seq(1:22)

chrend <- cumsum(chrlen)
chrstart <- c(0,chrend[1:21]) + 1 

chrcum <- as.data.frame(cbind(names(chrlen),chrstart,chrend))
colnames(chrcum) <- c("chrnum","chrstart","chrend")
chrcum <- chrcum %>% mutate(chrnum = as.integer(chrnum), chrstart = as.double(chrstart), chrend = as.double(chrend))
save(chrcum,file="RData/ChrCumulative.RData")

# Region SNPs
regionSNPs <- fread("Data/regionSNPs.tsv",sep="\t") %>% as.data.frame()
colnames(regionSNPs) <- c("region","chr","start","end","label")
save(regionSNPs,file="RData/regionSNP.RData")
