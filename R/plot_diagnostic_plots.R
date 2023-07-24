plot_diagnostic_plots <- function(count_matrix_raw, metadata_raw) {
  metadata <- metadata_raw %>%
    dplyr::mutate(sample_id = gsub("_T.*", "", sample_id)) %>%
    dplyr::mutate(condition = factor(condition))

  count_matrix <- count_matrix_raw %>%
    dplyr::select(Geneid, dplyr::all_of(gsub("_T.*", "", metadata[["sample_id"]]))) %>%
    tibble::column_to_rownames(var = "Geneid")

  # Create DESeq2 object
  if (all(colnames(count_matrix) == metadata$sample_id) == TRUE) {
    print("Sample order as in metadata - processing...")
  } else {
    stop("Samples not in order as in metadata!!")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~1
  )
  # Prefilter for low count genes
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  # Perform filtering
  dds_f <- dds[keep, ]

  # https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
  # Normalized counts using DESeq2's median of ratios method allowing for sample-sample comparisions
  # (only other tool capable of this kind of normalization is EdgeR)
  dds_norm <- DESeq2::estimateSizeFactors(dds)
  counts_norm <- DESeq2::counts(dds_norm, normalized = TRUE)

  # PCA
  # The running times are shorter when using blind=FALSE and if the function
  # DESeq has already been run, because then it is not necessary to re-estimate
  # the dispersion values.
  vsd <- DESeq2::vst(dds_f, blind = FALSE)

  pcaData <- DESeq2::plotPCA(vsd,
    intgroup = c("sample_id", "condition"),
    ntop = 500,
    returnData = TRUE
  )

  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # by condition
  plot_groups <- ggplot2::ggplot(pcaData, ggplot2::aes(
    x = PC1, y = PC2,
    col = condition,
    sample_id = sample_id
  )) +
    ggplot2::geom_point(size = 2) +
    # ggsci::scale_color_lancet() +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()

  plot_groups %>%
    plotly::ggplotly(tooltip = c("sample_id", "Condition")) %>%
    plotly::layout(margin = list(t = 75)) %>%
    htmlwidgets::saveWidget("PCA_interactive.html")

  # Boxplots ------------------------------------------------------
  df <- reshape2::melt(counts_norm,
    variable.name = "Samples",
    value.name = "count"
  ) %>% # reshape the matrix
    dplyr::rename("sample_id" = "Var2") %>%
    dplyr::left_join(metadata %>% dplyr::select(sample_id, condition)) %>%
    dplyr::arrange(sample_id)

  q <- c(.25, .5, .75)
  df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(
      quant25 = quantile(count, probs = q[1]),
      quant50 = quantile(count, probs = q[2]),
      quant75 = quantile(count, probs = q[3])
    )

  df$sample_id <- factor(df$sample_id, levels = metadata$sample_id)

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = sample_id, y = log2(count + 1), fill = condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::ggtitle("Boxplot of normalized counts") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()

  p2
  ggplot2::ggsave("boxplots.png",
    width = 8,
    height = 10,
    units = "in"
  )

  # Density -----------------------------------------------------------------

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = log(count + 1), colour = sample_id, fill = sample_id)) +
    # ylim(c(0, 0.17)) +
    ggplot2::geom_density(alpha = 0.2, linewidth = 1.25) +
    ggplot2::facet_wrap(~condition) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::xlab(expression(log[2](count + 1))) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 7, face = "bold"),
      legend.text = ggplot2::element_text(size = 7, face = "bold")
    ) +
    ggplot2::theme(legend.key.size = ggplot2::unit(1, "line")) +
    ggplot2::ggtitle("Density of normalized counts") +
    ggplot2::theme_bw()

  p3
  ggplot2::ggsave("density.png",
    width = 15,
    height = 12,
    units = "in"
  )

  p3 %>%
    plotly::ggplotly() %>%
    plotly::layout(margin = list(t = 75)) %>%
    htmlwidgets::saveWidget("density_interactive.html")
}
