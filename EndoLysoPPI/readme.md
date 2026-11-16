# Endo-Lysosomal Protein Interactome Construction

## Overview

This repository contains resources used for the construction of an endo-lysosomal protein interactome. It integrates previously published endo-lysosomal protein interactomes with gene/protein lists obtained from Gene Ontology (GO) annotations, covering both lysosomal and endosomal compartments and their associated biological processes.

The resources are organized into two main sections:

1. Previously published endo-lysosomal protein interactomes.
2. Gene/protein lists collected from Gene Ontology Cellular Component (GO:CC) and Biological Process (GO:BP) annotations.

---

## 1. Previously Published Interactomes

Two previously published resources were collected and used as reference datasets for endo-lysosomal protein interactome construction.

### 1.1. First published interactome

- **Number of proteins:** 709
- **Source:** [Nature Communications](https://www.nature.com/articles/s41467-022-33951-0)

### 1.2. Second published interactome

- **Number of proteins:** 768
- **Source:** [Journal of General Physiology](https://rupress.org/jgp/article-pdf/158/4/e202413745/2036727/jgp_202413745.pdf)

**Note:** The bibliographic details of the second reference should be verified, as the supplied publication link does not appear to correspond to an interactome dataset.

### 1.3. Shared proteins between the two interactomes

Comparison of the two published protein lists identified **594 shared proteins**.

| Dataset | Number of proteins |
|:---|---:|
| First published interactome | 709 |
| Second published interactome | 768 |
| Shared proteins | **594** |

---

## 2. Resources of Endo-Lysosomal Gene/Protein Lists

Gene and protein lists associated with the endo-lysosomal system were collected using Gene Ontology (GO) annotations.

### 2.1. Gene Ontology Cellular Component (GO:CC)

Protein and protein-complex annotations were retrieved using the following GO:CC terms. The organism filter was set to *Homo sapiens*.

#### Lysosome

- **GO term:** `GO:0005764`
- **GO term name:** Lysosome
- **Organism:** *Homo sapiens*
- **Number of proteins:** 816
- **Number of protein complexes:** 32

#### Endosome

- **GO term:** `GO:0005768`
- **GO term name:** Endosome
- **Organism:** *Homo sapiens*
- **Number of proteins:** 2,148
- **Number of protein complexes:** 52

The reported counts are after removal of duplicates.

| GO:CC term | GO ID | Organism | Proteins | Protein complexes |
|:---|:---|:---|---:|---:|
| Lysosome | GO:0005764 | *Homo sapiens* | 816 | 32 |
| Endosome | GO:0005768 | *Homo sapiens* | 2,148 | 52 |

### 2.2. Gene Ontology Biological Process (GO:BP)

GO:BP terms were retrieved using keyword searches.

#### Lysosome-associated biological processes

- **Search keyword:** `lysosome`
- **Number of GO:BP terms returned:** 81
- **Output file:** `GOBP_Lysosome.txt`

#### Endosome-associated biological processes

- **Search keyword:** `endosome`
- **Number of GO:BP terms returned:** 96
- **Output file:** `GOBP_Endosome.txt`

Genes/proteins annotated to the retrieved GO:BP terms were extracted using the **GO.db** R package.

The extraction was performed using the following R script:

- `R_GOdb.R`

### 2.3. Key files

| File | Description |
|:---|:---|
| `GOBP_Lysosome.txt` | GO:BP terms retrieved using the keyword search "lysosome" (81 terms). |
| `GOBP_Endosome.txt` | GO:BP terms retrieved using the keyword search "endosome" (96 terms). |
| `R_GOdb.R` | R script using the GO.db package to extract genes/proteins annotated to the selected GO:BP terms. |

---
