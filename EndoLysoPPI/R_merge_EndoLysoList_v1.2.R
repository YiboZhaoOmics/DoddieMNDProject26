# Finalise gene list of Endo-Lysosomal PPI network v1.0 
# Merge with 2 previously published interactome v1.1
# UniprotQC - Reviewed only v1.2

setwd("EndoLysoPPI")
L_CC <- read.delim("GO0005764_gene.txt")
E_CC <- read.delim("GO0005768_gene.txt")
L_BP <- read.delim("lysosome_GOBP_gene_proteinonly.txt")
E_BP <- read.delim("endosome_GOBP_gene_proteinonly.txt")

genelist <- unique(c(unlist(strsplit(paste0(L_BP$annotation, collapse = ","), split = ",")),
                     unlist(strsplit(paste0(E_BP$annotation, collapse = ","), split = ",")),
                     L_CC$GeneSymbol,
                     E_CC$GeneSymbol))

write.table(genelist, "EndoLysoPPI_genelist_v1.0.txt", sep = "\t", row.names = F, quote = F)

## manually plotted Venn gram from https://www.biotools.fr/misc/venny
# Checked overlapped and unique genes across the 3 gene/protein list
## v1.1: include the unique parts of Singh and Nanarova PPI list (N = 113)
add_genes <- read.delim("Unique_Singh+Nanarova.txt", header = F)
genelist <- c(genelist, add_genes$V1)
write.table(genelist, "EndoLysoPPI_genelist_v1.1.txt", sep = "\t", row.names = F, quote = F)

## Uniprot QC: Keep only proteins with Reviewed profile
# Note: Uniprot is sensitive to special symbols in the gene name such as HLA-A
library(httr)
library(dplyr)
library(readr)

get_reviewed_genes <- function(symbols, taxon = 9606, batch_size = 50) {
        
        symbols <- unique(trimws(as.character(symbols)))
        symbols <- symbols[!is.na(symbols) & symbols != "" & symbols != "NA"]
        # remove characters that would break the quoted string
        symbols <- gsub('"', "", symbols)
        
        batches <- split(symbols, ceiling(seq_along(symbols) / batch_size))
        
        res <- lapply(batches, function(b) {
                
                gene_query <- paste0('gene_exact:"', b, '"', collapse = " OR ")
                query <- sprintf("(%s) AND (organism_id:%s) AND (reviewed:true)",
                                 gene_query, taxon)
                
                response <- GET(
                        "https://rest.uniprot.org/uniprotkb/search",
                        query = list(
                                query  = query,
                                fields = "accession,gene_primary,gene_names,reviewed",
                                format = "tsv",
                                size   = 500
                        )
                )
                
                if (http_error(response)) {
                        stop("UniProt error ", status_code(response), ": ",
                             content(response, "text", encoding = "UTF-8"),
                             "\nFirst symbols in failing batch: ",
                             paste(head(b, 5), collapse = ", "), call. = FALSE)
                }
                
                read_tsv(content(response, "text", encoding = "UTF-8"),
                         show_col_types = FALSE)
        })
        
        bind_rows(res)
}
reviewed <- get_reviewed_genes(genelist)
UniprotQC <- data.frame("symbol" = genelist,
                        "reviewed" = ifelse(genelist %in% reviewed, "yes", "no"))
write.table(UniprotQC, "EndoLysoPPI_genelist_UniProtQC.txt", sep = "\t", quote = F, row.names = F)
write.table(reviewed, "EndoLysoPPI_genelist_v1.2.txt", sep = "\t", quote = F, row.names = F)
