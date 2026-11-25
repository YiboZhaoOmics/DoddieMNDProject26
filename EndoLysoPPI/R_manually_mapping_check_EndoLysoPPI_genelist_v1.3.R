#########################################################
### Touble Shooting: Mismatched Entrez and Uniprot ID ###
#########################################################
## entrez-entrez edge ID, no duplicates
## uniprot-uniprot ID, 35 duplicates
## while using org.hs.eg.db r package to extract entrez ID, for multiple matches, we kept the first entrez ID
## some proteins are mismatched with uniprot and may thereby affect downstram PPI matching
## update ID list of EndoLysoPPI edge list and node list, using uniprot ID as primary
## generate gene list v1.3 and network edge table v1.1
setwd("EndoLysoPPI")
ELlist <- read.delim("EndoLysoPPI_genelist_v1.2.txt")
ELPPI_merge <- read.delim("ELPPINv1_edge.txt")

library(httr2)
library(dplyr)
library(readr)

base_url <- "https://rest.uniprot.org"

# Unique Entrez IDs from your PPI network
entrez_ids <- unique(c(
        as.character(ELPPI_merge$EntrezA),
        as.character(ELPPI_merge$EntrezB)
))

entrez_ids <- entrez_ids[
        !is.na(entrez_ids) & entrez_ids != ""
]

# 1. Submit mapping job
response <- request(paste0(base_url, "/idmapping/run")) %>%
        req_body_form(
                from = "GeneID",
                to = "UniProtKB",
                ids = paste(entrez_ids, collapse = ",")
        ) %>%
        req_perform()

job_id <- resp_body_json(response)$jobId

# 2. Wait for completion
repeat {
        status <- request(
                paste0(base_url, "/idmapping/status/", job_id)
        ) %>%
                req_perform() %>%
                resp_body_json()
        
        if (!is.null(status$jobStatus)) {
                if (status$jobStatus == "FINISHED") break
                if (status$jobStatus %in% c("ERROR", "FAILED")) {
                        stop("UniProt mapping job failed")
                }
        } else {
                break
        }
        
        Sys.sleep(3)
}

# 3. Download results with annotations
result <- request(
        paste0(base_url, "/idmapping/uniprotkb/results/stream/", job_id)
) %>%
        req_url_query(
                format = "tsv",
                fields = "accession,gene_primary,reviewed,organism_id"
        ) %>%
        req_perform()

mapping_result <- read.delim(text = httr2::resp_body_string(result))

head(mapping_result)
mapping_reviewed <- subset(mapping_result, Reviewed == "reviewed")

# update ID matching in both gene list and network
ELPPI_merge$SwissA <- mapping_reviewed$Entry[match(ELPPI_merge$EntrezA, mapping_reviewed$From)]
ELPPI_merge$SwissB <- mapping_reviewed$Entry[match(ELPPI_merge$EntrezB, mapping_reviewed$From)]

ELPPI_merge$GeneNameA <- mapping_reviewed$Gene.Names..primary.[match(ELPPI_merge$EntrezA, mapping_reviewed$From)]
ELPPI_merge$GeneNameB <- mapping_reviewed$Gene.Names..primary.[match(ELPPI_merge$EntrezB, mapping_reviewed$From)]
missed_ID <- unique(c(ELPPI_merge$EntrezA[is.na(ELPPI_merge$SwissA)], ELPPI_merge$EntrezB[is.na(ELPPI_merge$SwissB)]))

# check missed ID manually
missed_ID[i]
ELPPI_merge$SwissA[ELPPI_merge$EntrezA==missed_ID[i]] <- "A1A4F0"
ELPPI_merge$SwissB[ELPPI_merge$EntrezB==missed_ID[i]] <- "A1A4F0"
ELPPI_merge$GeneNameA[ELPPI_merge$EntrezA==missed_ID[i]] <- "SLC66A1LP"
ELPPI_merge$GeneNameB[ELPPI_merge$EntrezB==missed_ID[i]] <- "SLC66A1LP"
i = i+1
# 8681 with unreviewed uniprot ID, remove
ELPPI_merge <- ELPPI_merge[ELPPI_merge$EntrezA!=8681,]
ELPPI_merge <- ELPPI_merge[ELPPI_merge$EntrezB!=8681,]
summary(is.na(ELPPI_merge))

## for pseudogene, use gene name corresponding uniprot instead.
ELlist_updated <- rbind(
        setNames(
                ELPPI_merge[, c("EntrezA", "SwissA", "GeneNameA")],
                c("Entrez", "Entry", "PrimaryGeneNames")
        ),
        setNames(
                ELPPI_merge[, c("EntrezB", "SwissB", "GeneNameB")],
                c("Entrez", "Entry", "PrimaryGeneNames")
        )
)
ELlist_updated <- unique(ELlist_updated)

# manually check duplicated matches
ELlist_updated$Entrez[duplicated(ELlist_updated$Entry)]
ELlist_updated$Entry[duplicated(ELlist_updated$Entry)]
ELlist_updated$PrimaryGeneNames[duplicated(ELlist_updated$Entry)]

# 721 -> C4B -> P0C0L5
# 728358: Unreviewed, remove

ELlist_updated$Entry[ELlist_updated$Entrez==1667] <- "P59665"
ELlist_updated$Entry[ELlist_updated$Entrez==1667] <- "DEFA1"

ELPPI_merge$SwissA[ELPPI_merge$EntrezA==1667] <- "P59665"
ELPPI_merge$SwissB[ELPPI_merge$EntrezB==1667] <- "P59665"
ELPPI_merge$GeneNameA[ELPPI_merge$EntrezA==1667] <- "DEFA1"
ELPPI_merge$GeneNameB[ELPPI_merge$EntrezB==1667] <- "DEFA1"


ELlist_updated <- ELlist_updated[ELlist_updated$Entrez!=728358,] 

ELPPI_merge <- ELPPI_merge[ELPPI_merge$EntrezA!=728358,] 
ELPPI_merge <- ELPPI_merge[ELPPI_merge$EntrezB!=728358,] 

write.table(ELlist_updated, "EndoLysoPPI_genelist_v1.3.txt", sep = "\t", quote = F, row.names = F)
write.table(ELPPI_merge, "ELPPINv1.1_edge.txt", sep = "\t", quote = F, row.names = F)
