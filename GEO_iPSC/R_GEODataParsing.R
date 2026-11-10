# text parsing from GEO search results (N=954) // 13/Oct/26
library(stringr)
library(dplyr)

parse_geo <- function(path) {
        text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
        
        # Split into records at lines starting with "1. ", "2. ", ...
        chunks <- str_trim(str_split(text, "(?m)^\\s*\\d+\\.\\s+(?=\\S)")[[1]])
        chunks <- chunks[str_detect(chunks, "Accession:")]
        
        grab <- function(chunk, pattern, group = 2) {
                m <- str_match(chunk, pattern)[, group]
                if (is.na(m)) "" else str_trim(m)
        }
        
        rows <- lapply(chunks, function(chunk) {
                lines <- str_split(chunk, "\n")[[1]]
                title <- str_trim(lines[1])
                
                # Abstract = text between the title and the "Organism:" line
                body <- paste0("\n", paste(lines[-1], collapse = "\n"))
                abstract <- str_split(body, "\n\\s*Organism:", n = 2)[[1]][1]
                abstract <- str_squish(abstract)
                abstract <- str_remove(abstract, "^\\(Submitter supplied\\)\\s*")
                abstract <- str_remove(abstract, "\\s*more\\.\\.\\.$")
                
                plat <- str_match(chunk, "Platform:\\s*(.+?)\\s+(\\d+)\\s+Samples?")
                
                data.frame(
                        Title      = title,
                        Abstract   = abstract,
                        Organism   = grab(chunk, "Organism:\\s*(.+)"),
                        Type       = grab(chunk, "Type:\\s*(.+)"),
                        Platform   = if (is.na(plat[1, 2])) "" else str_trim(plat[1, 2]),
                        SampleSize = if (is.na(plat[1, 3])) NA_integer_ else as.integer(plat[1, 3]),
                        FTP        = grab(chunk, "FTP download:\\s*(?:\\S+\\s+\\(\\w+\\)\\s+)?((?:ftp|https?)://\\S+)"),
                        Series     = grab(chunk, "(?m)^\\s*(\\w+)\\s+Accession:"),
                        Accession  = grab(chunk, "Accession:\\s*(\\w+)"),
                        ID         = grab(chunk, "\\bID:\\s*(\\d+)"),
                        stringsAsFactors = FALSE
                )
        })
        
        do.call(rbind, rows)
}

df <- parse_geo("~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCLPDRA2026/Doodie/gds_result.txt")

##### Quality Control #####
df$text <- paste(df$Title, df$Abstract)
# 1.Add MND type columns 
als_pattern <- regex(
        "\\bALS\\b|amyotrophic lateral sclerosis|\\bfALS\\b|\\bsALS\\b|Lou Gehrig",
        ignore_case = FALSE   # keeps "ALS" from matching the ordinary word "als"
)
als_pattern2 <- regex("amyotrophic lateral sclerosis|lou gehrig", ignore_case = TRUE)

ftd_pattern <- regex(
        "\\bFTD\\b|\\bFTLD\\b|frontotemporal (dementia|lobar degeneration)|fronto-temporal dementia",
        ignore_case = TRUE
)

df$ALS   <- str_detect(df$text, als_pattern) | str_detect(df$text, als_pattern2)
df$FTD   <- str_detect(df$text, ftd_pattern)
df$Other <- !df$ALS & !df$FTD
# table(df$Other) N = 488 other MNDs
df_trimmed <- subset(df, !Other)

# 2.Search for iPSC sets and add cell types
## 2.1 Look for iPSC pattern 
ipsc_pattern <- regex(
        "\\b(h?i-?PSCs?|iPS cells?|induced pluripotent stem cells?|induced pluripotent)\\b",
        ignore_case = TRUE
)
df_trimmed$iPSC <- str_detect(df_trimmed$text, ipsc_pattern)
table(df_trimmed$iPSC) # N=68 iPSC

# 2.2: iPSC mention + the next 60 characters
ipsc_ctx <- regex(
        "(?:h?i-?PSCs?|iPS cells?|induced pluripotent stem cells?).{0,60}",
        ignore_case = TRUE
)

## 2.3: Build cell type name dictionary. Abbreviations are case-sensitive, full names are not.
celltype_dict <- c(
        "Microglia"             = "\\bMGs?\\b|(?i:microglia)",
        "Motor neurons"         = "\\bMNs?\\b|(?i:motor neurons?)",
        "Cortical neurons"      = "(?i:cortical (excitatory |glutamatergic )?neurons?)",
        "Dopaminergic neurons"  = "(?i:dopaminergic neurons?)|\\bDANs?\\b",
        "Sensory neurons"       = "(?i:sensory neurons?)",
        "Astrocytes"            = "(?i:astrocytes?)",
        "Oligodendrocytes"      = "(?i:oligodendrocytes?)|\\bOPCs?\\b",
        "Neural progenitors"    = "\\bNPCs?\\b|\\bNSCs?\\b|(?i:neural (progenitor|stem) cells?)",
        "Organoids"             = "(?i:organoids?)",
        "Neurons (unspecified)" = "(?i:neurons?)"
)

extract_celltypes <- function(txt) {
        ctx <- paste(str_extract_all(txt, ipsc_ctx)[[1]], collapse = " | ")
        if (ctx == "") return("")
        hits <- names(celltype_dict)[sapply(celltype_dict, function(p) str_detect(ctx, p))]
        # drop the generic "neurons" label if a specific neuron type was found
        specific <- setdiff(grep("neuron", hits, value = TRUE), "Neurons (unspecified)")
        if (length(specific) > 0) hits <- setdiff(hits, "Neurons (unspecified)")
        paste(hits, collapse = "; ")
}

df_trimmed$iPSC_celltypes <- sapply(df_trimmed$text, extract_celltypes, USE.NAMES = FALSE)
df_trimmed$iPSC_context   <- sapply(df_trimmed$text, function(t)
        paste(str_extract_all(t, ipsc_ctx)[[1]], collapse = " | "), USE.NAMES = FALSE)

df_trimmed[,c("Accession", "iPSC_celltypes", "iPSC_context")]

head(subset(df_trimmed, !iPSC&Type == "Non-coding RNA profiling by high throughput sequencing"))
table(df_trimmed$Type)

# 3. Dataset types (array/seq/non-coding)
type_dict <- c(
        ExpressionArray = "Expression profiling by array",
        ExpressionSeq   = "Expression profiling by high throughput sequencing",
        GenomeSeq       = "Genome binding/occupancy profiling by high throughput sequencing",
        GenomeArray     = "Genome variation profiling by array",
        NonCodeSeq      = "Non-coding RNA profiling by high throughput sequencing",
        NonCodeArray    = "Non-coding RNA profiling by array",
        MethyArray      = "Methylation profiling by (genome tiling )?array",
        MethySeq        = "Methylation profiling by high throughput sequencing"
)

for (nm in names(type_dict)) {
        df_trimmed[[nm]] <- str_detect(df_trimmed$Type, regex(type_dict[[nm]], ignore_case = TRUE))
}

# Other = none of the categories above matched
df_trimmed$Other <- rowSums(df_trimmed[, names(type_dict)]) == 0

# Check
df_trimmed[, c("Accession", "Type", names(type_dict), "Other")]

write.csv(df_trimmed, "geo_structured.csv", row.names = FALSE, fileEncoding = "UTF-8")
