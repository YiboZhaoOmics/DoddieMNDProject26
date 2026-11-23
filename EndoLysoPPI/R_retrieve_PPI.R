# Extract PPIs for EndoLysoListv1,2 from HIPPIE, MIST and PINOT

library(readr)
setwd("~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie")
EL <- read.delim("EndoLysoPPI/EndoLysoPPI_genelist_v1.2.txt")


################
#### HIPPIE ####
################
source("Rsource/HIPPIE_network_query.R")
HIPPIE <- read.delim("EndoLysoPPI/PPI/HIPPIE-current.txt") ## query with UniprotID

## HIPPIE queried with uniprot ID
ELPPI_HIPPIE <- HIPPIE_network_query(HIPPIE, EL$Entry)
hist(ELPPI_HIPPIE$score)
rm(HIPPIE)


##############
#### MIST ####
##############

source("Rsource/MIST_network_query.R")
MIST <- read.delim("EndoLysoPPI/PPI/01_28_26_TAX9606_ppi_records.tsv", header = F) ## quety with entrez ID

## MIST queried with entrez ID
# BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
# library(org.Hs.eg.db)
# entrez <- mapIds(org.Hs.eg.db,
#               keys      = EL$Entry,
#               keytype   = "UNIPROT",
#               column    = "ENTREZID",
#               multiVals = "first") # to avoid one:many matching, the first entrez ID is used
# EL$Entrez <- entrez
# write.table(EL, "EndoLysoPPI/EndoLysoPPI_genelist_v1.2.txt", sep = "\t", quote = F, row.names = F)
ELPPI_MIST <- MIST_network_query(MIST, EL$Entrez)
rm(MIST)
## remove comment columns
ELPPI_MIST <- ELPPI_MIST[,1:16]
colnames(ELPPI_MIST) <- c("MasterNetID", "TaxID_A", "GeneA", "TaxID_B", "GeneB", "Rank", "Interaction_type", "Exp_Direct", "Exp_Indirect", "TaxID_interolog", "Source_databases", "Reference", "Reference_type", "Count_direct", "Count_indirect", "Count_paper")
head(ELPPI_MIST)
summary(ELPPI_HIPPIE$score>=0.72)


###############
#### PINOT ####
###############
## need to manually set filter = lenient or stringent
## strigent results are saved in a separate file
# setwd("/Users/yibozhao/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie/PINOT/Dictionary")
# library("httr")
# source("PINOT_unified_REST_fixed_interseed.R")
setwd("EndoLysoPPI/")
ELPPI_PINOT <- read.delim("PINOT_seed_to_seed_network_stringent.tsv")
hist(ELPPI_PINOT$Final.Score)
summary(ELPPI_PINOT$Final.Score>=2)

ELPPI <- list(HIPPIE = ELPPI_HIPPIE,
              MIST = ELPPI_MIST,
              PINOT_Stringent = ELPPI_PINOT
)
saveRDS(ELPPI, "~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie/EndoLysoPPI/ELPPI.rds")
