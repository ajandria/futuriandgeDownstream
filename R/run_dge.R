#' Run Differential Gene Expression Analysis
#'
#' This function runs a DESeq2-based differential gene expression analysis
#' on a raw count matrix and associated metadata.
#'
#' @param count_matrix_raw A data frame representing the raw count matrix.
#' It must contain a column named "Geneid" that includes gene identifiers.
#' @param metadata_raw A data frame representing the metadata.
#' It must contain a column named "sample" with sample identifiers
#' and a column named "condition" with condition identifiers.
#' @param comparisons A list of vectors. Each vector includes two elements:
#' the identifiers of two conditions to be compared.
#'
#' @return A data frame representing the results of the differential gene expression analysis.
#' The data frame is sorted by the "Geneid" and "Comparison" columns in ascending order.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr select arrange filter bind_rows mutate
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map_dfr
#'
#' @export
#'
run_dge <- function(count_matrix_raw, metadata_raw, comparisons) {
  # Remove _T.* from sample_id and convert condition to factor
  metadata <- metadata_raw %>%
    dplyr::mutate(sample_id = gsub("_T.*", "", .data$sample_id)) %>%
    dplyr::mutate(condition = factor(.data$condition))

  # Select Geneid and all columns with names as in metadata[["sample_id"]]
  count_matrix <- count_matrix_raw %>%
    dplyr::select(.data$Geneid, dplyr::all_of(gsub("_T.*", "", metadata[["sample_id"]]))) %>%
    tibble::column_to_rownames(var = "Geneid")

  # Check if samples are in the same order in count_matrix and metadata
  if (all(colnames(count_matrix) == metadata$sample_id) == TRUE) {
    print("Sample order as in metadata - processing...")
  } else {
    stop("Samples not in order as in metadata!!")
  }

  # Create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~condition
  )

  # Prefilter for low count genes
  keep <- rowSums(DESeq2::counts(dds)) >= 10

  # Perform filtering
  dds_f <- dds[keep, ]

  # Perform DE
  dds_deseq <- DESeq2::DESeq(dds_f)

  # Normalized counts using DESeq2's median of ratios method allowing for sample-sample comparisons
  dds_norm <- DESeq2::estimateSizeFactors(dds)
  counts_norm <- DESeq2::counts(dds_norm, normalized = TRUE)

  # Prepare normalized matrix for export
  counts_norm_sorted <- counts_norm %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Geneid") %>%
    dplyr::arrange(.data$Geneid)

  # remove X from first colnames starting with numbers...
  colnames(counts_norm_sorted) <- gsub("^X", "", colnames(counts_norm_sorted))

  deseq2_results_list <- list()

  # Loop over comparisons
  for (i in 1:nrow(comparisons)) {
    comparison_i <- comparisons[i, ]
    studied_group <- comparison_i[[1]]
    control_group <- comparison_i[[2]]
    comparison_name <- paste0(studied_group, "_vs_", control_group)

    # Perform DESeq2 analysis and store in results list
    deseq2_results_list[[comparison_name]] <- DESeq2::results(dds_deseq,
                                                              contrast = c(
                                                                "condition",
                                                                strsplit(comparison_name, "_vs_")[[1]][1],
                                                                strsplit(comparison_name, "_vs_")[[1]][2]
                                                              )
    )

    # Prepare for metadata addition
    compared_groups <- gsub(".*condition ", "", deseq2_results_list[[comparison_name]]@elementMetadata[2, 2]) %>%
      stringr::str_replace_all(" ", "_") %>%
      stringr::str_split("_vs_") %>%
      unlist()

    compared_samples <- metadata %>%
      dplyr::filter(.data$condition %in% compared_groups)

    # Add metadata to results
    deseq2_results_list[[comparison_name]]$compared_samples <- c(compared_samples$sample_id, rep(NA, length(deseq2_results_list[[comparison_name]]$log2FoldChange) - length(compared_samples$sample_id)))
    deseq2_results_list[[comparison_name]]$compared_conditions <- c(as.character(compared_samples$condition), rep(NA, length(deseq2_results_list[[comparison_name]]$log2FoldChange) - length(compared_samples$sample_id)))
    deseq2_results_list[[comparison_name]]$analysis_description <- c(deseq2_results_list[[comparison_name]]@elementMetadata$description, rep(NA, length(deseq2_results_list[[comparison_name]]$log2FoldChange) - length(deseq2_results_list[[comparison_name]]@elementMetadata$description)))
  }

  # Create a final object to be returned
  final_object <- list(
    dge_results = deseq2_results_list,
    normalised_counts = dplyr::tibble(counts_norm_sorted),
    raw_counts_no_gtf_miRNA = count_matrix_raw,
    raw_metadata = metadata_raw,
    comparisons = comparisons
  )

  return(final_object)
}
