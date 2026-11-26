# ENTPD1 expression and lung function in GSE32537.
# Inputs: data/lung_function/GSE32537/GSE32537_geneSymbol_mean.csv and
#         GSE32537_lung_function.tsv (clean numeric fvc, dlco, total_score columns).
# Output: results/05_lung_function_entpd1_gse32537/ boxplot and correlation panels.
# Method: exclude dlco == 0; median split in matched clinical samples;
# Wilcoxon rank-sum (exact = FALSE, correct = FALSE); Spearman; lm visual guides.
# Packages: ggpubr, ggplot2, cowplot; original versions were not saved.
# Source audit: R05. The required cleaned inputs and their preprocessing code
# were absent from the project folder. Raw GEO matrix/clinical extracts differ.
# Boxplot labels now identify the measured endpoint; statistical formulas are retained.

source("R/_paths.R")
gene <- "ENTPD1"
correlation_method <- "spearman"
endpoints <- c("fvc", "dlco", "total_score")
input_dir <- file.path(DATA_DIR, "lung_function", "GSE32537")
expression_file <- file.path(input_dir, "GSE32537_geneSymbol_mean.csv")
clinical_file <- file.path(input_dir, "GSE32537_lung_function.tsv")
require_inputs(c(expression_file, clinical_file))
require_packages(c("ggpubr", "ggplot2", "cowplot"))
expression <- read.csv(expression_file, row.names = 1)
clinical <- read.table(clinical_file, sep = "\t", header = TRUE, row.names = 1)
require_columns(clinical, endpoints, "Clean clinical table")
require_numeric(clinical, endpoints, "Clean clinical table")
require_gene(expression, gene)
clinical <- clinical[clinical$dlco != 0, ]
if (!nrow(clinical)) stop("No clinical samples remain after dlco != 0.")
if (!all(rownames(clinical) %in% colnames(expression))) {
  stop("All retained clinical row names must match expression sample columns.")
}
score <- colMeans(expression[gene, rownames(clinical), drop = FALSE], na.rm = TRUE)
score_table <- data.frame(geo = names(score), expression = score)
threshold <- median(score, na.rm = TRUE)
score_table$GROUP <- ifelse(score >= threshold, "ENTPD1_high", "ENTPD1_low")
clinical$geo <- rownames(clinical)
clinical <- merge(clinical, score_table, by = "geo")
groups <- levels(factor(clinical$GROUP))
if (length(groups) != 2) stop("Both expression groups are required for the comparison.")
comparisons <- list(groups)

set.seed(1) # Reproducible jitter positions; no effect on statistical tests.
boxplots <- lapply(endpoints, function(endpoint) {
  ggboxplot(clinical, x = "GROUP", y = endpoint, color = "GROUP",
            xlab = "ENTPD1 expression group", ylab = endpoint,
            palette = c("#2E9FDF", "#FF0F00"), add = "jitter") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       method.args = list(exact = FALSE, correct = FALSE),
                       label = "p.format", tip.length = 0.02, label.size = 3)
})
legend <- get_legend(boxplots[[1]] + theme(legend.position = "right"))
box_panels <- lapply(boxplots, function(p) p + theme(legend.position = "none"))
combined_boxplots <- plot_grid(
  plot_grid(plotlist = box_panels, ncol = 3, align = "hv"),
  legend, rel_widths = c(3, 0.4)
)
correlation_plots <- lapply(endpoints, function(endpoint) {
  ggscatter(clinical, x = "expression", y = endpoint,
            color = "#2E9FDF", shape = 16, size = 1.5, alpha = 0.6,
            add = "reg.line", add.params = list(color = "#FF6666", fill = "lightgray"),
            conf.int = TRUE, fullrange = TRUE) +
    stat_cor(method = correlation_method, label.x.npc = "left",
             label.y.npc = "top", size = 4) +
    labs(x = "ENTPD1 Expression", y = endpoint) + theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(),
          axis.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5))
})
combined_correlations <- plot_grid(plotlist = correlation_plots, ncol = 3,
                                   align = "hv", labels = "AUTO")
output_dir <- workflow_results("05_lung_function_entpd1_gse32537")
ggsave(file.path(output_dir, "entpd1_gse32537_boxplots.pdf"),
       combined_boxplots, width = 10, height = 4)
ggsave(file.path(output_dir, "entpd1_gse32537_correlations.pdf"),
       combined_correlations, width = 12, height = 4)
ggsave(file.path(output_dir, "entpd1_gse32537_correlations.png"),
       combined_correlations, width = 12, height = 4, dpi = 300)
record_session(output_dir)
