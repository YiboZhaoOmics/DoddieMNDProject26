# GEO iPSC Sets Assembly

## 1. Project Overview

This repository contains a manually curated collection of ALS/FTD-related induced pluripotent stem cell (iPSC) datasets obtained from the Gene Expression Omnibus (GEO).

## 2. Systematic Data Search (14 Oct 2026)

Dataset identification and collection follow a systematic search and screening process:

1. **Search Term**: Search term: ((("motor neuron disease"[MeSH Terms] OR motor neuron disease[All Fields]) OR ("amyotrophic lateral sclerosis"[MeSH Terms] OR Amyotrophic Lateral Sclerosis[All Fields])) OR (frontal[All Fields] AND temporal[All Fields] AND ("dementia"[MeSH Terms] OR dementia[All Fields]))) AND "Homo sapiens"[Organism] AND ("gds"[Filter] AND "human"[Organism] OR "Homo sapiens"[Organism])
2. **Text Parsing**: Transform the query results to a structured dataframe. (R_GEODataParsing.R)
3. **Dataset Screening and Selection**: From the 954 results retrieved by the query term, we retained 466 ALS/FTD entries. Among these sets, 68 involved iPSC-derived cell lines.
4. **Manually Curate iPSC sets**: Datasets were manually checked for 3 levels of data: Differential Expression Analysis (DEA) Results, Normalised Expression Levels, Raw Expression Data. For datasets with no clear documentation, we tried to contact the author for more information.
5. **Data Download**: Download eligible datasets and their associated metadata, document relevant dataset information, and organise the files using a consistent naming convention for subsequent analyses.

## 3. Repository File Structure

The repository contains documentation and supporting files, alongside datasets collected from GEO and other publicly available sources.

### 3.1 Documentation and Reference Files

| File name | Description |
|---|---|
| `iPSC_CheckList.xlsx` | Checklist for tracking dataset identification, screening, eligibility, and collection status. |
| `iPSC_FileDictionary.xlsx` | File-name dictionary documenting dataset identifiers, file names, descriptions, and relevant file information. |
| `R_GEODataParsing.R` | R scripts for text parsing and mining of raw GEO search results.|

### 3.2 Data Files

The remaining files in this repository comprise datasets and associated metadata derived from GEO and other publicly available sources.
The goal is to keep the highest level of processed data as possible. In particular: DEA > Normalised > Raw

## 4. Notes
- Original data files should be preserved wherever possible.
- Original dataset identifiers and source information should be retained wherever possible to ensure traceability and reproducibility.
