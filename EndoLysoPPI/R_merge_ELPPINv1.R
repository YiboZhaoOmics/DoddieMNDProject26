# merge PPIs retrieved from 3 datasets
## HIPPIE and MIST needs QC -> remove duplicates including A->B B->A
ELPPI <- readRDS("EndoLysoPPI/ELPPI.rds")
ELPPI_HIPPIE <- ELPPI$HIPPIE
ELPPI_MIST <- ELPPI$MIST
ELPPI_PINOT <- ELPPI$PINOT_Stringent
rm(ELPPI)
library(dplyr)

head(ELPPI_HIPPIE)
head(ELPPI_MIST)
head(ELPPI_PINOT)

create_pair_list <- function(A,B){ # generate a new id for PPI with entrez IDs of interactor A and B, keep the smaller one on the left and larger one on the right, connecting with "-"
        pair_list <- paste(pmin(A,B), pmax(A,B), sep = "-")
        return(pair_list)
}

pair_list_HIPPIE <- create_pair_list(ELPPI_HIPPIE$entrez_A, ELPPI_HIPPIE$entrez_B)
summary(duplicated(pair_list_HIPPIE)) # N = 32 duplicates
ELPPI_HIPPIE <- ELPPI_HIPPIE[!duplicated(pair_list_HIPPIE),]
pair_list_HIPPIE <- pair_list_HIPPIE[!duplicated(pair_list_HIPPIE)]

pair_list_MIST <- create_pair_list(ELPPI_MIST$GeneA, ELPPI_MIST$GeneB)
summary(duplicated(pair_list_MIST))
ELPPI_MIST <- ELPPI_MIST[!duplicated(pair_list_MIST),] # N = 50905
pair_list_MIST <- pair_list_MIST[!duplicated(pair_list_MIST)]

pair_list_PINOT <- create_pair_list(ELPPI_PINOT$EntrezA, ELPPI_PINOT$EntrezB)
summary(duplicated(pair_list_PINOT)) # PINOT QC is completed 

# venn graph
library(ggVennDiagram)
library(ggplot2)

# Combine three vectors into a named list
venn_list <- list(
        HIPPIE = pair_list_HIPPIE,
        MIST = pair_list_MIST,
        PINOT = pair_list_PINOT
)

# Plot Venn diagram
p <- ggVennDiagram(
        venn_list,
        label = "count",
        label_alpha = 0
) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = "Overlap between three platforms") +
        theme(plot.title = element_text(hjust = 0.5))

print(p)

# # overwrite the original raw PPI .rda with QC-ed one
# ELPPI <- list(HIPPIE = ELPPI_HIPPIE,
#               MIST = ELPPI_MIST,
#               PINOT_Stringent = ELPPI_PINOT
# )
# saveRDS(ELPPI, "~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie/EndoLysoPPI/ELPPI.rds")

#############
### Merge ###
#############
pair_list_unique <- unique(c(pair_list_PINOT, pair_list_MIST, pair_list_HIPPIE))
pair_split <- do.call(rbind, strsplit(pair_list_unique, "-"))
ELPPI_merge <- data.frame("EntrezA" = as.numeric(pair_split[,1]),
                          "EntrezB" = as.numeric(pair_split[,2]),
                          "HIPPIE" = ifelse(is.na(match(pair_list_unique,pair_list_HIPPIE)), 0, 1),
                          "MIST" = ifelse(is.na(match(pair_list_unique,pair_list_MIST)), 0, 1),
                          "PINOT" = ifelse(is.na(match(pair_list_unique,pair_list_PINOT)), 0, 1),
                          "Hscore" = ELPPI_HIPPIE$score[match(pair_list_unique, pair_list_HIPPIE)],
                          "Mscore" = ELPPI_MIST$Rank[match(pair_list_unique, pair_list_MIST)],
                          "Pscore" = ELPPI_PINOT$Final.Score[match(pair_list_unique, pair_list_PINOT)],
                          "PairID" = pair_list_unique)

# # save pair_list as searching ID in the 3 individual sets to make later search easier
# ELPPI_HIPPIE$PairID <- pair_list_HIPPIE
# ELPPI_PINOT$PairID <- pair_list_PINOT
# ELPPI_MIST$PairID <- pair_list_MIST
# 
# ELPPI <- list(HIPPIE = ELPPI_HIPPIE,
#               MIST = ELPPI_MIST,
#               PINOT_Stringent = ELPPI_PINOT
# )
# saveRDS(ELPPI, "~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie/EndoLysoPPI/ELPPI.rds")
setwd("EndoLysoPPI")
write.table(ELPPI_merge, "ELPPINv1_edge.txt", sep = "\t", row.names = F, quote = F)
