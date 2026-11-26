# ADORA2B expression versus a TGF-beta signaling expression score across bulk cohorts.
# Inputs: data/bulk_expression/<dataset>/<dataset>_matrix.txt; MSigDB metadata resource
# resources/gene_sets/HALLMARK_TGF_BETA_SIGNALING.v2026.1.Hs.txt (54 source symbols).
# Output: results/07_bulk_tgfb_score_adora2b/ individual and combined correlations.
# Method: intersect available genes; per-gene Z-score within a cohort; sample mean
# across these Z-scores; zero-SD genes contribute zero; require at least 2 genes.
# Correlation: Spearman; lm line and confidence band are visual guides.
# Packages: ggplot2, ggpubr, cowplot; original versions were not saved.
# Source audit: R07 root and nested byte versions share this analysis body.
# The source filename identified ENTPD1, but both saved versions use ADORA2B.
# The earlier ENTPD1/TGF-beta figures have no saved generation version in this project.

source("R/_paths.R")
gene <- "ADORA2B"
score_label <- "TGF-beta signaling score (Mean Z-score)"
input_dir <- file.path(DATA_DIR, "bulk_expression")
gene_set_file <- file.path(RESOURCE_DIR, "gene_sets", "HALLMARK_TGF_BETA_SIGNALING.v2026.1.Hs.txt")
metagene <- read_msigdb_symbols(gene_set_file)
if (length(metagene) != 54) stop("Expected 54 source symbols; verify the specified gene-set version.")
file_list <- list.files(input_dir, pattern = "_matrix[.]txt$", full.names = TRUE,
                        recursive = TRUE)
# Raw GEO series matrices are separate inputs, not processed gene-symbol matrices.
file_list <- file_list[!grepl("_series_matrix[.]txt$", file_list)]
if (!length(file_list)) {
  stop("No processed *_matrix.txt files found in data/bulk_expression. See R/README.md.")
}
dataset_names <- basename(dirname(file_list))
if (anyDuplicated(dataset_names)) stop("Each bulk matrix needs a unique dataset directory name.")
require_packages(c("ggplot2", "ggpubr", "cowplot"))
output_dir <- workflow_results("07_bulk_tgfb_score_adora2b")
plot_list <- list()
cor_values <- c()
summary_rows <- list()

for (file_path in file_list) {
  dataset <- basename(dirname(file_path))
  expression <- read.csv(file_path, sep = "\t", check.names = FALSE, row.names = 1)
  if (!gene %in% rownames(expression)) {
    message("Skipping ", dataset, ": missing ", gene)
    next
  }
  require_gene(expression, gene)
  present_genes <- intersect(metagene, rownames(expression))
  if (length(present_genes) < 2) {
    message("Skipping ", dataset, ": fewer than 2 score genes.")
    next
  }
  sub_matrix <- expression[present_genes, , drop = FALSE]
  if (anyNA(sub_matrix)) stop(dataset, ": source scoring formula requires non-missing score-gene values.")
  scaled_matrix <- t(apply(sub_matrix, 1, function(x) {
    if (sd(x) == 0) return(rep(0, length(x)))
    scale(x)
  }))
  score <- colMeans(scaled_matrix)
  plot_data <- data.frame(expression = as.numeric(expression[gene, ]), score = score)
  correlation <- cor.test(plot_data$expression, plot_data$score, method = "spearman")
  cor_values[dataset] <- correlation$estimate
  summary_rows[[dataset]] <- data.frame(
    dataset = dataset, samples = nrow(plot_data), genes_in_score = length(present_genes),
    spearman_rho = unname(correlation$estimate), p_value = correlation$p.value
  )

  p <- ggplot(plot_data, aes(x = expression, y = score)) +
    geom_point(color = "#34495e", alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", color = "red", fill = "lightgray", linewidth = 0.8) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(title = dataset, subtitle = paste("Genes in Score:", length(present_genes)),
         x = "ADORA2B Expression", y = score_label) + theme_bw() +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(output_dir, paste0(dataset, "_correlation.pdf")), p, width = 4, height = 4)

  plot_list[[dataset]] <- ggplot(plot_data, aes(x = expression, y = score)) +
    geom_point(color = "#34495e", alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", color = "red", fill = "lightgray", linewidth = 0.6) +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) +
    labs(title = dataset, x = "ADORA2B", y = "TGF-beta score") + theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.title = element_text(size = 8), axis.text = element_text(size = 7))
}

sorted_names <- names(sort(cor_values, decreasing = TRUE))
if (!length(sorted_names)) stop("No cohorts produced a finite Spearman correlation.")
combined <- plot_grid(plotlist = plot_list[sorted_names], ncol = 5, align = "hv")
total_rows <- ceiling(length(sorted_names) / 5)
ggsave(file.path(output_dir, "adora2b_tgfb_score_correlations.pdf"), combined,
       width = 20, height = total_rows * 4, limitsize = FALSE)
ggsave(file.path(output_dir, "adora2b_tgfb_score_correlations.png"), combined,
       width = 20, height = total_rows * 4, limitsize = FALSE)
write.csv(do.call(rbind, summary_rows), file.path(output_dir, "correlation_summary.csv"), row.names = FALSE)
record_session(output_dir)
