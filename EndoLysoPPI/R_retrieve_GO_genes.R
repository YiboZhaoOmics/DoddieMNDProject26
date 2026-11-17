# extract genes/proteins annotated to GO-BP associated with lysosome and endosome
# note: returned list includes protein_complex
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

source("get_GO_genes.R") # return genes annotated to individual GO terms via QuickGO API
source("quickgo_functions.R") # return log file for GO terms with empty results

retrieve_GO_genes <- function(filename, GOdf){ # colnames: "ID", "Name"
        gene_list <- integer()
        go_terms <- GOdf$ID
        for (i in seq_along(go_terms)) {
                
                GO_results <- tryCatch(
                        map_dfr(
                                go_terms[i],
                                get_GO_genes,
                                taxon = 9606
                        ),
                        error = function(e) {
                                return(data.frame())
                        }
                )
                
                genes <- paste0(unique(GO_results$symbol), collapse = ",")
                gene_list <- c(gene_list, genes)
                #Sys.sleep(2)
        }
        outtable <- data.frame(GOdf,
                               "annotation" = gene_list)
        
        failed_GO <- subset(outtable, annotation == "")
        sapply(failed_GO$ID, is_obsolete)
        GO_check <- check_GO_terms(
                failed_GO$ID,
                taxon = 9606
        )
        write.table(subset(outtable, annotation != ""),paste0(filename, ".txt"), sep = "\t", quote = F, row.names = F)
        write.table(GO_check, paste0(filename,"_log", ".txt"), sep = "\t", quote = F, row.names = F)
        
}


lysosome <- read.delim("EndoLysoPPI/GOBP_Lysosome.txt")
endosome <- read.delim("EndoLysoPPI/GOBP_Endosome.txt")
setwd("EndoLysoPPI/")
retrieve_GO_genes("lysosome_GOBP_gene", lysosome)
retrieve_GO_genes("endosome_GOBP_gene", endosome)
