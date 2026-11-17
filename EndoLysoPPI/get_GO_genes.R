# get_GO_genes.R
# Fetch gene annotations for a GO term from the EBI QuickGO REST API.
#
# Usage:
#   source("get_GO_genes.R")
#   res <- get_GO_genes("GO:0006915")                # human (taxon 9606), all pages
#   res <- get_GO_genes("GO:0006915", max_pages = 1) # first page only
#   unique(res$symbol)                               # distinct gene symbols

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tibble)
})

get_GO_genes <- function(go_id, taxon = 9606, limit = 100, max_pages = Inf) {

  go_id <- trimws(as.character(go_id))
  url   <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"

  pages       <- list()
  page        <- 1
  total_pages <- 1

  while (page <= total_pages && page <= max_pages) {

    response <- GET(
      url,
      query = list(
        goId    = go_id,
        taxonId = taxon,
        limit   = limit,
        page    = page
      ),
      accept("application/json")
    )

    stop_for_status(response)

    dat <- fromJSON(
      content(response, "text", encoding = "UTF-8"),
      flatten = TRUE
    )

    if (length(dat$results) == 0) break

    pages[[page]] <- as_tibble(dat$results)

    # total number of pages reported by the API (fall back to current page)
    total_pages <- if (is.null(dat$pageInfo$total)) page else dat$pageInfo$total
    page <- page + 1
  }

  if (length(pages) == 0) {
    return(tibble())
  }

  bind_rows(pages) %>%
    mutate(query_GO = go_id)
}
