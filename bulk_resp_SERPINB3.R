## ─────────────────────────────────────────────────────────────
## RUN ON: local Mac
##   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
##   Rscript bulk_resp_SERPINB3.R
## SERPINB3 in TCGA-HNSC bulk RNA-seq (STAR TPM, log2(TPM+1), primary tumor -01),
## cetuximab Non-responder (PD/SD) vs Responder (CR/PR). Wilcoxon.
## Input : Source_Data/Bulk/TCGA-HNSC.star_tpm.tsv (only the SERPINB3 row is read)
##         ../tcga_hnsc_cetuximab_patients.csv
## Output: ../Result/bulk_resp_SERPINB3_boxplot.tiff
## ─────────────────────────────────────────────────────────────
library(ggplot2)
base <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics"
tpm_file <- file.path(base, "Source_Data/Bulk/TCGA-HNSC.star_tpm.tsv")
GENE_ID <- "ENSG00000057149"   # SERPINB3

hdr <- strsplit(readLines(tpm_file, n = 1), "\t")[[1]]
row <- system(sprintf("grep '^%s' '%s'", GENE_ID, tpm_file), intern = TRUE)
v <- as.numeric(strsplit(row, "\t")[[1]][-1]); names(v) <- hdr[-1]
v <- v[substr(names(v), 14, 15) == "01"]
v <- tapply(v, substr(names(v), 1, 12), mean)

grp <- read.csv(file.path(base, "cetuximab_improve/DeepspotM/tcga_hnsc_cetuximab_patients.csv"))
grp$group <- ifelse(grepl("Complete Response|Partial Response", grp$cetuximab_response), "R",
             ifelse(grepl("Progressive Disease|Stable Disease", grp$cetuximab_response), "NR", NA))
grp <- grp[!is.na(grp$group) & grp$patient_id %in% names(v), ]
grp$SERPINB3 <- v[grp$patient_id]
grp$group <- factor(grp$group, levels = c("NR", "R"))
print(table(grp$group))
p <- wilcox.test(SERPINB3 ~ group, data = grp)$p.value
cat(sprintf("SERPINB3 median NR=%.2f, R=%.2f, Wilcoxon p=%.3f\n",
            median(grp$SERPINB3[grp$group == "NR"]), median(grp$SERPINB3[grp$group == "R"]), p))

g <- ggplot(grp, aes(group, SERPINB3, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.15, size = 0.8) +
  scale_fill_manual(values = c(NR = "#E41A1C", R = "#377EB8")) +
  labs(x = "Cetuximab response", y = "SERPINB3 log2(TPM+1)", title = sprintf("p = %.3f", p)) +
  theme_classic(base_size = 7) + theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5))
ggsave(file.path(base, "cetuximab_improve/DeepspotM/Result/bulk_resp_SERPINB3_boxplot.tiff"),
       g, width = 1.6, height = 2, units = "in", dpi = 600, device = "tiff", compression = "lzw")
