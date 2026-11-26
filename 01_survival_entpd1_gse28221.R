# ENTPD1 survival analysis in GSE28221 (Affymetrix samples).
# Inputs: data/survival/GSE28221/{GSE28221_matrix.txt,GSE28221_sampleInfo.csv}.
# Output: results/01_survival_entpd1_gse28221/ Kaplan-Meier PDFs and session info.
# Method: single-gene expression, median split (>= median is High), log-rank test.
# Packages: survival, survminer, ggplot2; original package versions were not saved.
# Source audit: R01, SHA256 f7baa516e54df33ca750d6565d50c480988b29bec3c5f0311f2ce257a8ba290f.
# Audit note: the original filename identified the wrong dataset; inputs were GSE28221.
# Time units were not established in the source. Set time_label after verification.

source("R/_paths.R")
gene <- "ENTPD1"
time_label <- "Time (source units)"
input_dir <- file.path(DATA_DIR, "survival", "GSE28221")
expression_file <- file.path(input_dir, "GSE28221_matrix.txt")
clinical_file <- file.path(input_dir, "GSE28221_sampleInfo.csv")
require_inputs(c(expression_file, clinical_file))
require_packages(c("survival", "survminer", "ggplot2"))

expression <- read.table(expression_file, sep = "\t", header = TRUE,
                         row.names = 1, check.names = FALSE)
clinical <- read.csv(clinical_file, row.names = 1)
require_columns(clinical, c("platform", "surTime", "surStatus"), "Clinical table")
require_numeric(clinical, c("surTime", "surStatus"), "Clinical table")
require_gene(expression, gene)
clinical$geo <- rownames(clinical)
clinical <- clinical[clinical$platform == "Affymetrix", ]
if (!nrow(clinical)) stop("No clinical samples have platform == 'Affymetrix'.")
if (!all(rownames(clinical) %in% colnames(expression))) {
  stop("All retained clinical row names must match expression sample columns.")
}
gene_expression <- expression[gene, rownames(clinical), drop = FALSE]
score <- colMeans(gene_expression, na.rm = TRUE)
score_table <- data.frame(geo = names(score), expression = score)
threshold <- median(score, na.rm = TRUE)
score_table$GROUP <- ifelse(score_table$expression >= threshold,
                            "ENTPD1_High", "ENTPD1_Low")
clinical <- merge(clinical, score_table, by = "geo")

fit <- survfit(Surv(surTime, surStatus) ~ GROUP, data = clinical)
p <- ggsurvplot(
  fit, data = clinical, pval = TRUE, pval.method = TRUE,
  conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",
  size = 1, palette = c("#E41A1C", "#377EB8"), surv.median.line = "hv",
  legend.title = "Expression", legend.labs = c("ENTPD1 High", "ENTPD1 Low"),
  xlab = time_label, ylab = "Survival Probability",
  ggtheme = theme_bw() + theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10)
  )
)
p$plot <- p$plot + theme(aspect.ratio = 1)
output_dir <- workflow_results("01_survival_entpd1_gse28221")
pdf(file.path(output_dir, "entpd1_gse28221_survival_4x6.pdf"), width = 4, height = 6)
print(p)
dev.off()
pdf(file.path(output_dir, "entpd1_gse28221_survival_5.5x5.5.pdf"), width = 5.5, height = 5.5)
print(p)
dev.off()
record_session(output_dir)
