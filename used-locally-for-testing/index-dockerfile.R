
# Setup -------------------------------------------------------------------

downstream_ch_file <- "/futuriandgeDownstream/used-locally-for-testing/counts/extended_source-pipeline-out_meta_file-dockerfile.txt"
organism <- 'mus_musculus'

user_file <- "/futuriandgeDownstream/used-locally-for-testing/template-user-request.xlsx"

sample_metadata <- readxl::read_excel(user_file, sheet = 'metadata')
contrasts <- readxl::read_excel(user_file, sheet = 'comparisons')

# Load --------------------------------------------------------------------
library(futuriandgeDownstream)

# F1
counts <- return_count_matrix(downstream_ch_file)

# F2
counts_no_mirnas <- remove_mirna(counts, organism)

# F3
dge_results <- run_dge(count_matrix_raw = counts_no_mirnas,
                       metadata_raw = sample_metadata,
                       comparisons = contrasts)

# save dge results
save(dge_results, file = 'dge_results.RData')

# F4
plot_diagnostic_plots(count_matrix_raw = counts_no_mirnas,
                      metadata_raw = sample_metadata)

# F5
contrasts_for_report <- contrasts %>%
  dplyr::mutate(compared_gorups = paste0(
    studied_effect,
    '_vs_',
    baseline
  ))

for (groups in contrasts_for_report$compared_gorups) {

  sample_metadata$sample <- gsub("_T1",
                                 "",
                                 sample_metadata$sample_id)

  render_dge_html_report(
    comparison_name = groups,
    dge_results_in = data.frame(dge_results$dge_results[[groups]]),
    metadata = sample_metadata,
    norm_counts = dge_results$normalised_counts,
    raw_counts = dge_results$raw_counts_no_gtf_miRNA,
    organism = organism,
    output_name = paste0(groups, ".html")
  )
}

