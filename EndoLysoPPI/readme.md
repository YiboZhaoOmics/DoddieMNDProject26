# Endo-Lysosomal Protein Interactome Construction

## 1. Analysis Overview

The endo-lysosomal system plays an essential role in intracellular trafficking, protein degradation, and cellular homeostasis. Constructing a comprehensive protein–protein interaction (PPI) network of this system will help characterize its molecular organization and identify potential interactions among proteins involved in endosomal and lysosomal functions.

Two previously published endo-lysosomal PPI datasets were collected as reference resources:

- **Published interactome 1:** 709 proteins ([Nature Communications](https://www.nature.com/articles/s41467-022-33951-0))
- **Published interactome 2:** 768 proteins ([Journal of General Physiology](https://rupress.org/jgp/article-pdf/158/4/e202413745/2036727/jgp_202413745.pdf))

These published interactomes primarily focused on proteins associated with endo-lysosomal cellular compartments (Cellular Component, CC). To establish a more comprehensive network, this project expands the protein list by incorporating genes/proteins associated with both **Cellular Component (CC)** and **Biological Process (BP)** annotations from Gene Ontology (GO).

The resulting protein list is subsequently curated and used to construct an integrated endo-lysosomal PPI network using interaction data from multiple major repositories and platforms.

## 2. Methodology Overview

The construction of the endo-lysosomal protein interactome consists of four main steps:

### 1. Protein List Retrieval (EndoLyso Gene List v1.0)

Genes/proteins associated with the endo-lysosomal system were retrieved from Gene Ontology (GO) annotations, covering both Cellular Component (CC) and Biological Process (BP).

- **GO:CC:** Retrieved proteins and protein complexes annotated to:
  - Lysosome (`GO:0005764`): 816 proteins and 32 protein complexes.
  - Endosome (`GO:0005768`): 2,148 proteins and 52 protein complexes.
- **GO:BP:** Retrieved genes/proteins associated with lysosome- and endosome-related biological processes using the respective keywords in GO:BP. Annotations were retrieved via the QuickGO API.
  - Lysosome: 81 GO:BP terms retrieved.
  - Endosome: 96 GO:BP terms retrieved.

The CC- and BP-associated gene/protein lists were combined, resulting in **EndoLyso Gene List v1.0 (2482 entries)** after removing duplicates and protein complexes.

### 2. Integration with Published Interactomes (EndoLyso Gene List v1.1)

The EndoLyso Gene List v1.0 was compared and integrated with the two previously published endo-lysosomal PPI protein lists.

An additional 113 genes/proteins were identified, generating **EndoLyso Gene List v1.1 (2595 genes/proteins)**.

### 3. Protein Quality Control Using UniProt (EndoLyso Gene List v1.2)

The protein list was further curated using UniProt to retain only proteins with a **Reviewed** status, ensuring that the final list consists of proteins with manually reviewed UniProt entries.

This generated **EndoLyso Gene List v1.2 (2428 entries)**, which was used as the input for PPI extraction.

### 4. PPI Data Extraction and Network Construction (EndoLyso PPIN v1)

Protein–protein interactions among the curated protein list were retrieved from three major PPI repositories and platforms (num. of entries after QC):

- **PINOTv1.1 (Stringent)**: 51177
- **HIPPIEv3.0**: 52515
- **MISTv6.0**: 51924

The retrieved interaction datasets were merged to construct the integrated **EndoLyso Protein–Protein Interaction Network (EndoLyso PPIN v1.0)** (N = 56473 entries).
After remove duplicates (N = 82) due to multi-matched IDs, 56391 edges remained in the **EndoLyso PPIN v1.1**
Out of these edges, 2308 nodes remained in the **EndoLyso Gene List v1.3**

## 3. File Dictionary

| File | Description |
|:---|:---|
| `R_retrieve_GO_genes.R` | R script for extracting genes/proteins annotated to the selected GO:BP terms via the QuickGO API. |
| `quickgo_functions.R` | Source R functions for GO term quality control. |
| `get_GO_genes.R` | Source R functions for connecting to the QuickGO API and retrieving annotated genes/proteins to generate EndoLyso Gene List v1.0. |
| `R_merge_EndoLysoList_v1.2.R` | R script for merging the GO-derived gene/protein lists with the two published interactome lists to generate EndoLyso Gene List v1.2. |
| `R_retrieve_PPI` | R script for retrieve PPIs from HIPPIEv3.0, MIST_humanV6.0, PINOTv1.1. |
| `R_merge_ELPPINv1.0.R` | R script for merge PPI lists from the 3 platforms. |
| `R_manually_mapping_check_EndoLysoPPI_genelist_v1.3.R` | R script for the final QC to remove duplicates due to multi-matched IDs to generate Gene List v1.3 & PPINv1.1. |
