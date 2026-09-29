# ============================================================================
# RUN ON: local Mac
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
#   Rscript tcga25_topgenes_barplot.R
# Top 10 genes, tumor-enriched pseudobulk, High risk vs Low risk (25 TCGA
# cetuximab patients). Layout from Win4221 topgenes_barplot.R, but UNPAIRED:
# group means across patients, ranked by |mean_High - mean_Low|.
#   LEFT  = up in High risk, RIGHT = up in Low risk
# 4.1 x 1.5 inch, 600 dpi TIFF, Arial 7pt
# Input:  ../Result/tcga25_tumor_pseudobulk.csv  (tcga25_tumor_pseudobulk.py)
#         ../Result/cetuximab_25_OS_groups.csv   (KM_cetuximab25_OS_risk.R)
# Output: ../Result/tcga25_topgenes_tumor_High_vs_Low.tiff
#         ../Result/tcga25_DE_tumor_High_vs_Low.csv  (all genes, with Wilcoxon p)
# ============================================================================

library(dplyr)

base_dir   <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
result_dir <- file.path(base_dir, "Result")
N <- 10

pb  <- read.csv(file.path(result_dir, "tcga25_tumor_pseudobulk.csv"),
                row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(result_dir, "cetuximab_25_OS_groups.csv"))
grp <- grp[grp$patient_id %in% rownames(pb), ]
cat(sprintf("High=%d, Low=%d\n", sum(grp$group == "High"), sum(grp$group == "Low")))

gene_cols <- setdiff(names(pb), c("n_tumor_spots", "n_total_spots"))
matH <- as.matrix(pb[grp$patient_id[grp$group == "High"], gene_cols])
matL <- as.matrix(pb[grp$patient_id[grp$group == "Low"],  gene_cols])

# ---- Unpaired DE: A = Low, B = High ----
mean_A <- colMeans(matL, na.rm = TRUE)
mean_B <- colMeans(matH, na.rm = TRUE)
de <- data.frame(gene = gene_cols, mean_A = mean_A, mean_B = mean_B,
                 log2FC = log2((mean_B + 1e-6) / (mean_A + 1e-6)),
                 stringsAsFactors = FALSE, row.names = NULL)
de$wilcox_p <- sapply(gene_cols, function(g)
  suppressWarnings(wilcox.test(matH[, g], matL[, g])$p.value))
de$padj <- p.adjust(de$wilcox_p, method = "BH")
de <- de[is.finite(de$log2FC), ]

out_de <- de %>% rename(mean_Low = mean_A, mean_High = mean_B) %>%
  mutate(expr_diff = abs(mean_High - mean_Low)) %>% arrange(-expr_diff)
write.csv(out_de, file.path(result_dir, "tcga25_DE_tumor_High_vs_Low.csv"), row.names = FALSE)

# ---- Plot top N barplot (2-panel), same style as template ----
plot_topgenes <- function(de_df, labelA, labelB, colA, colB, out_file) {
  de_df$expr_diff <- abs(de_df$mean_B - de_df$mean_A)
  up <- de_df %>% filter(log2FC > 0) %>% arrange(-expr_diff) %>% head(N)
  dn <- de_df %>% filter(log2FC < 0) %>% arrange(-expr_diff) %>% head(N)

  global_max <- max(c(up$expr_diff, dn$expr_diff), na.rm = TRUE)
  up_norm <- up$expr_diff / global_max
  dn_norm <- dn$expr_diff / global_max

  tmp_png <- sub("\\.tiff$", "_tmp.png", out_file)
  png(tmp_png, width = 4.1, height = 1.5, units = "in", res = 600, type = "quartz",
      family = "Arial", pointsize = 7)
  par(mfrow = c(1, 2), oma = c(0, 0, 0, 0),
      mar = c(3.2, 2.9, 1.1, 0.3), mgp = c(1.8, 0.4, 0))

  # Panel 1: Up in B (High risk)
  bp <- barplot(up_norm, col = colB, border = NA,
                ylab = "Expression score",
                main = paste0("Up in ", labelB), col.main = colB,
                cex.main = 0.9, font.main = 2,
                ylim = c(0, 1.15), axes = FALSE)
  axis(2, las = 1, tcl = -0.2, lwd = 0.3, lwd.ticks = 0.3)
  box(bty = "l", lwd = 0.3)
  text(bp, -0.04, labels = up$gene, srt = 45, adj = 1, xpd = TRUE, cex = 0.75)

  # Panel 2: Up in A (Low risk)
  bp <- barplot(dn_norm, col = colA, border = NA,
                ylab = "Expression score",
                main = paste0("Up in ", labelA), col.main = colA,
                cex.main = 0.9, font.main = 2,
                ylim = c(0, 1.15), axes = FALSE)
  axis(2, las = 1, tcl = -0.2, lwd = 0.3, lwd.ticks = 0.3)
  box(bty = "l", lwd = 0.3)
  text(bp, -0.04, labels = dn$gene, srt = 45, adj = 1, xpd = TRUE, cex = 0.75)

  dev.off()
  system(sprintf("sips -s format tiff -s dpiWidth 600 -s dpiHeight 600 '%s' --out '%s'", tmp_png, out_file))
  file.remove(tmp_png)
  cat("Saved:", out_file, "\n")

  cat("\nTop up in", labelB, ":\n"); print(up[, c("gene", "mean_A", "mean_B", "log2FC", "wilcox_p")])
  cat("\nTop up in", labelA, ":\n"); print(dn[, c("gene", "mean_A", "mean_B", "log2FC", "wilcox_p")])
}

plot_topgenes(de, "Low risk", "High risk", "#377EB8", "#E41A1C",
              file.path(result_dir, "tcga25_topgenes_tumor_High_vs_Low.tiff"))
cat("Done!\n")
