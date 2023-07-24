# This has been generated on:
#> Sys.time()
# [1] "2023-07-19 13:02:07 CEST"
#> packageVersion('biomaRt')
# [1] ‘2.56.1’

# Specify paths to input GTF files for human and mouse
homo_sapiens_gtf <- "data-raw/ensembl_110/Homo_sapiens.GRCh38.110.chr.gtf"
mus_musculus_gtf <- "data-raw/ensembl_110/Mus_musculus.GRCm39.110.chr.gtf"

#' Generate biomaRt Annotation
#'
#' This function imports a GTF file, retrieves annotations from the corresponding BioMart,
#' and renames columns for clarity.
#'
#' @param gtf_path A character string indicating the path to a GTF file.
#' @param organism A character string indicating the organism dataset in BioMart.
#'
#' @return A data frame representing the annotated table from BioMart.
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom rtracklayer import
#' @importFrom dplyr rename
generate_biomaRt_annotation <- function(gtf_path, organism) {
  # Import GTF file and convert to DataFrame
  gtf <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf)

  # Use the correct BioMart for the specified organism
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = organism)

  # Get annotations from BioMart
  biomaRt_annot_table <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "external_gene_source",
      "description",
      "go_id",
      "name_1006",
      "namespace_1003"
    ),
    filters = "ensembl_gene_id",
    values = gtf_df$gene_id,
    mart = ensembl
  )

  # Rename columns for clarity and return the annotated table
  dplyr::rename(biomaRt_annot_table, go_name = "name_1006", go_ontology = "namespace_1003")
}

# Generate annotation for human
biomaRt_annotation_homo_sapiens <- generate_biomaRt_annotation(homo_sapiens_gtf, "hsapiens_gene_ensembl")

# Generate annotation for mouse
biomaRt_annotation_mus_musculus <- generate_biomaRt_annotation(mus_musculus_gtf, "mmusculus_gene_ensembl")

# Save the data for internal use in the package
usethis::use_data(biomaRt_annotation_homo_sapiens, biomaRt_annotation_mus_musculus, internal = TRUE)
