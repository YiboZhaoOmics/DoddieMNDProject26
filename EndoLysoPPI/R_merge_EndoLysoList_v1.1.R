# Finalise gene list of Endo-Lysosomal PPI network v1.0 & extract interactions among them from HIPPIE & MIST
# PINOT was queried manually via online platform
# STRING score was also checked as a validation evidence
setwd("EndoLysoPPI")
L_CC <- read.delim("GO0005764_gene.txt")
E_CC <- read.delim("GO0005768_gene.txt")
L_BP <- read.delim("lysosome_GOBP_gene.txt")
E_BP <- read.delim("endosome_GOBP_gene.txt")

genelist <- unique(c(unlist(strsplit(paste0(L_BP$annotation, collapse = ","), split = ",")),
                     unlist(strsplit(paste0(E_BP$annotation, collapse = ","), split = ",")),
                     L_CC$GeneSymbol,
                     E_CC$GeneSymbol))

# remove protein complex
complex <- genelist[grepl("^[a-z]", genelist)] # N = 124
genelist <- genelist[!genelist %in% complex]
write.table(genelist, "EndoLysoPPI_genelist_v1.0.txt", sep = "\t", row.names = F, quote = F)

## manually plotted Venn gram from https://www.biotools.fr/misc/venny
# Checked overlapped and unique genes across the 3 gene/protein list
## v1.1: include the unique parts of Singh and Nanarova PPI list (N = 113)
add_genes <- read.delim("Unique_Singh+Nanarova.txt", header = F)
genelist <- c(genelist, add_genes$V1)
write.table(genelist, "EndoLysoPPI_genelist_v1.1.txt", sep = "\t", row.names = F, quote = F)

# save complex in a seperate file juse in case
write.table(complex, "EndoLysoPPI_complex.txt", sep = "\t", row.names = F, quote = F)
