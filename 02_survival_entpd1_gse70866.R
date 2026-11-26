# ENTPD1 survival analysis in GSE70866, platform GPL14550.
# Inputs: data/survival/GSE70866_GPL14550/{GSE70866_GPL14550_matrix.txt,
#         GSE70866_GPL14550_sampleInfo.csv}.
# Output: results/02_survival_entpd1_gse70866/ Kaplan-Meier PDF and session info.
# Method: above/below mean expression split; equal-to-mean samples excluded; log-rank.
# Packages: survival, survminer, ggplot2; original package versions were not saved.
# Source audit: R02. The original filename incorrectly identified SRSF1/GSE28221.
# DOCUMENTED LEGACY EVENT RECODING: surStatus is set to 0 when surTime > 1000,
# while surTime is unchanged. This preserves the source analysis and requires
# methodological confirmation before interpreting the resulting survival curve.
# No disease-group restriction was active in the source.

source("R/_paths.R")
gene <- "ENTPD1"
legacy_event_threshold <- 1000
time_label <- "Time (Days)"
input_dir <- file.path(DATA_DIR, "survival", "GSE70866_GPL14550")
expression_file <- file.path(input_dir, "GSE70866_GPL14550_matrix.txt")
clinical_file <- file.path(input_dir, "GSE70866_GPL14550_sampleInfo.csv")
require_inputs(c(expression_file, clinical_file))
require_packages(c("survival", "survminer", "ggplot2"))

expression <- read.table(expression_file, sep = "\t", header = TRUE, row.names = 1)
clinical <- read.csv(clinical_file, row.names = 1)
require_columns(clinical, c("geoID", "surTime", "surStatus"), "Clinical table")
require_numeric(clinical, c("surTime", "surStatus"), "Clinical table")
require_gene(expression, gene)
gene_expression <- as.data.frame(t(expression[gene, , drop = FALSE]))
gene_expression$geoID <- rownames(gene_expression)
clinical <- merge(clinical, gene_expression, by = "geoID")
if (!nrow(clinical)) stop("Clinical geoID values do not match expression sample columns.")
threshold <- mean(clinical$ENTPD1)
clinical$surStatus <- ifelse(clinical$surTime > legacy_event_threshold, 0,
                            clinical$surStatus)
clinical$GROUP <- ifelse(clinical$ENTPD1 > threshold, "ENTPD1_High",
                         ifelse(clinical$ENTPD1 < threshold, "ENTPD1_Low", " "))
clinical <- clinical[clinical$GROUP %in% c("ENTPD1_High", "ENTPD1_Low"), ]

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
output_dir <- workflow_results("02_survival_entpd1_gse70866")
pdf(file.path(output_dir, "entpd1_gse70866_gpl14550_survival_4x6.pdf"), width = 4, height = 6)
print(p)
dev.off()
record_session(output_dir)
