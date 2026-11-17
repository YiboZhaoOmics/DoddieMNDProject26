# QuickGO GO-term quality-control utilities
# Usage:
#   source("quickgo_functions.R")
#   GO_check <- check_GO_terms(go_terms, taxon = 9606)
#
# Checks:
#   - GO term existence, name, aspect and obsolete status
#   - Number of direct QuickGO annotation records for the requested taxon
#
# Note: annotation counts are records, not unique genes. Descendant terms
# are not included. QuickGO API errors are reported in the message column.

.check_quickgo_response <- function(response, context) {
  status <- httr::status_code(response)
  if (status != 200) {
    body <- httr::content(response, "text", encoding = "UTF-8")
    stop(
      context, " failed (HTTP ", status, "). Response: ",
      substr(body, 1, 500),
      call. = FALSE
    )
  }
}

check_GO_terms <- function(go_terms, taxon = 9606, pause = 1) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Please install required packages: ",
      "install.packages(c('httr', 'jsonlite', 'dplyr'))",
      call. = FALSE
    )
  }

  go_terms <- unique(trimws(as.character(go_terms)))
  go_terms <- go_terms[!is.na(go_terms) & nzchar(go_terms)]

  if (length(go_terms) == 0) {
    return(data.frame(
      GO_ID = character(),
      GO_name = character(),
      aspect = character(),
      obsolete = logical(),
      taxon = integer(),
      taxon_annotations = integer(),
      message = character(),
      stringsAsFactors = FALSE
    ))
  }

  invalid_format <- !grepl("^GO:[0-9]{7}$", go_terms)
  output <- vector("list", length(go_terms))
  names(output) <- go_terms

  # Check malformed IDs locally; query well-formed IDs in one ontology request.
  valid_format_ids <- go_terms[!invalid_format]
  term_records <- NULL
  ontology_error <- NULL

  if (length(valid_format_ids) > 0) {
    ontology_url <- paste0(
      "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",
      paste(valid_format_ids, collapse = ",")
    )

    ontology_response <- tryCatch(
      httr::GET(
        ontology_url,
        httr::accept("application/json"),
        httr::timeout(30)
      ),
      error = function(e) e
    )

    if (inherits(ontology_response, "error")) {
      ontology_error <- conditionMessage(ontology_response)
    } else if (httr::status_code(ontology_response) != 200) {
      ontology_error <- paste0(
        "Ontology API HTTP ", httr::status_code(ontology_response),
        ": ",
        substr(
          httr::content(ontology_response, "text", encoding = "UTF-8"),
          1, 300
        )
      )
    } else {
      ontology_json <- tryCatch(
        jsonlite::fromJSON(
          httr::content(ontology_response, "text", encoding = "UTF-8"),
          flatten = TRUE
        ),
        error = function(e) e
      )

      if (inherits(ontology_json, "error")) {
        ontology_error <- conditionMessage(ontology_json)
      } else {
        term_records <- ontology_json$results
      }
    }
  }

  for (i in seq_along(go_terms)) {
    id <- go_terms[i]

    if (invalid_format[i]) {
      output[[i]] <- data.frame(
        GO_ID = id, GO_name = NA_character_, aspect = NA_character_,
        obsolete = NA, taxon = as.integer(taxon),
        taxon_annotations = NA_integer_,
        message = "Invalid GO ID format",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!is.null(ontology_error)) {
      output[[i]] <- data.frame(
        GO_ID = id, GO_name = NA_character_, aspect = NA_character_,
        obsolete = NA, taxon = as.integer(taxon),
        taxon_annotations = NA_integer_,
        message = paste("Ontology API error:", ontology_error),
        stringsAsFactors = FALSE
      )
      next
    }

    term <- term_records[term_records$id == id, , drop = FALSE]

    if (nrow(term) == 0) {
      output[[i]] <- data.frame(
        GO_ID = id, GO_name = NA_character_, aspect = NA_character_,
        obsolete = NA, taxon = as.integer(taxon),
        taxon_annotations = NA_integer_,
        message = "GO term not found",
        stringsAsFactors = FALSE
      )
      next
    }

    obsolete <- isTRUE(term$isObsolete[1])
    go_name <- if ("name" %in% names(term)) as.character(term$name[1]) else NA_character_
    aspect <- if ("aspect" %in% names(term)) as.character(term$aspect[1]) else NA_character_

    if (obsolete) {
      output[[i]] <- data.frame(
        GO_ID = id, GO_name = go_name, aspect = aspect,
        obsolete = TRUE, taxon = as.integer(taxon),
        taxon_annotations = NA_integer_,
        message = "Obsolete GO term",
        stringsAsFactors = FALSE
      )
      next
    }

    # Check direct annotations for this GO ID and taxon.
    annotation_response <- tryCatch(
      httr::GET(
        "https://www.ebi.ac.uk/QuickGO/services/annotation/search",
        query = list(
          goId = id,
          taxonId = taxon,
          limit = 1,
          page = 1
        ),
        httr::accept("application/json"),
        httr::timeout(30)
      ),
      error = function(e) e
    )

    if (inherits(annotation_response, "error")) {
      annotation_count <- NA_integer_
      status_message <- paste(
        "Annotation API request error:",
        conditionMessage(annotation_response)
      )
    } else if (httr::status_code(annotation_response) != 200) {
      annotation_count <- NA_integer_
      status_message <- paste(
        "Annotation API HTTP",
        httr::status_code(annotation_response)
      )
    } else {
      annotation_json <- tryCatch(
        jsonlite::fromJSON(
          httr::content(annotation_response, "text", encoding = "UTF-8")
        ),
        error = function(e) e
      )

      if (inherits(annotation_json, "error") ||
          is.null(annotation_json$numberOfHits)) {
        annotation_count <- NA_integer_
        status_message <- "Could not read annotation count from API response"
      } else {
        annotation_count <- as.integer(annotation_json$numberOfHits)
        status_message <- if (annotation_count > 0) {
          "Valid - direct annotations found"
        } else {
          "Valid - no direct annotations for this taxon"
        }
      }
    }

    output[[i]] <- data.frame(
      GO_ID = id,
      GO_name = go_name,
      aspect = aspect,
      obsolete = FALSE,
      taxon = as.integer(taxon),
      taxon_annotations = annotation_count,
      message = status_message,
      stringsAsFactors = FALSE
    )

    if (i < length(go_terms) && pause > 0) Sys.sleep(pause)
  }

  dplyr::bind_rows(output)
}
