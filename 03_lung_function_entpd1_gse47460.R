# ENTPD1 expression and lung function in GSE47460, platform GPL6480.
# Inputs: data/lung_function/GSE47460_GPL6480/{GSE47460_GPL6480_matrix.txt,
#         GSE47460_GPL6480_sampleInfo.csv}.
# Output: results/03_lung_function_entpd1_gse47460/ boxplot and correlation panels.
# Method: median split across expression samples; Wilcoxon rank-sum tests
# (exact = FALSE, correct = FALSE); Spearman correlations; lm visual guides.
# Packages: ggpubr, ggplot2, cowplot; original versions were not saved.
# Source audit: R03, SHA256 261dbed999889cd8eff86329a45a3ed8992cbf434021618ea35e32b6cc7a44b0.
# The original PF-only filter was inactive: retain all clinical groups.
# Boxplot labels now identify the measured endpoint; statistical formulas are retained.

source("R/_paths.R")
gene <- "ENTPD1"
correlation_method <- "spearman"
endpoints <- c("predicted_fev1_pre_bd", "predicted_fvc_pre_bd",
               "predicted_fev1_post_bd", "predicted_fvc_post_bd", "predicted_dlco")
input_dir <- file.path(DATA_DIR, "lung_function", "GSE47460_GPL6480")
expression_file <- file.path(input_dir, "GSE47460_GPL6480_matrix.txt")
clinical_file <- file.path(input_dir, "GSE47460_GPL6480_sampleInfo.csv")
require_inputs(c(expression_file, clinical_file))
require_packages(c("ggpubr", "ggplot2", "cowplot"))
expression <- read.table(expression_file, sep = "\t", header = TRUE, row.names = 1)
clinical <- read.csv(clinical_file, row.names = 1)
require_columns(clinical, c("geoID", endpoints), "Clinical table")
require_numeric(clinical, endpoints, "Clinical table")
require_gene(expression, gene)

score <- colMeans(expression[gene, , drop = FALSE], na.rm = TRUE)
score_table <- data.frame(geo = names(score), expression = score)
threshold <- median(score, na.rm = TRUE)
score_table$GROUP <- ifelse(score >= threshold, "ENTPD1_high", "ENTPD1_low")
clinical$geo <- clinical$geoID
clinical <- merge(clinical, score_table, by = "geo")
if (!nrow(clinical)) stop("Clinical geoID values do not match expression sample columns.")
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
output_dir <- workflow_results("03_lung_function_entpd1_gse47460")
ggsave(file.path(output_dir, "entpd1_gse47460_gpl6480_boxplots.pdf"),
       combined_boxplots, width = 12, height = 8)
ggsave(file.path(output_dir, "entpd1_gse47460_gpl6480_correlations.pdf"),
       combined_correlations, width = 12, height = 8)
ggsave(file.path(output_dir, "entpd1_gse47460_gpl6480_correlations.png"),
       combined_correlations, width = 12, height = 8, dpi = 300)
record_session(output_dir)
