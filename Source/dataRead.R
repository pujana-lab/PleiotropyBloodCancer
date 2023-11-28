# read data and save in .RData files

library(data.table) # [1] ‘1.14.8’

# read pleiotropic SNP list (pleio_loci.tsv). Filter to conjfdr < 0.05
leadSNP <- fread("Data/pleio_loci.tsv") %>% as.data.frame() %>%
  rename(cancer_trait=trait1,
         blood_trait=trait2,
         pval_neoplasm=pval_trait1,
         pval_immunological=pval_trait2,
         zscore_neoplasm=zscore_trait1,
         zscore_immunological=zscore_trait2)
leadSNP005 <- leadSNP %>% filter(conjfdr < 0.05)
save(leadSNP005, file="Data/leadSNP005.RData")

# metadata info
metadata <- fread("Data/metadata_pleio.tsv") %>% as.data.frame()
save(metadata, file="Data/metadata.RData")

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
save(cancer.label,file="Data/CancerLabel.RData")

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
save(blood.label,file="Data/BloodLabel.RData")

