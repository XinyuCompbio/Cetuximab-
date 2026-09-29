# ============================================================================
# RUN ON: local Mac
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
#   python3 phoenix_resp_tumor_pseudobulk.py   # makes the input CSV
#   Rscript phoenix_resp_topgenes_barplot.R
# Top 10 genes, Phoenix 112um tumor-region pseudobulk (377-gene panel),
# Non-responder vs Responder (TCGA cetuximab patients). Layout from Win4221
# topgenes_barplot.R, UNPAIRED: group means, ranked by |mean_NR - mean_R|.
# Values are Phoenix log1p, so no per-patient total scaling.
#   LEFT = up in Non-responder, RIGHT = up in Responder
# 4.1 x 1.5 inch, 600 dpi TIFF, Arial 7pt
# Input:  ../Result/phoenix_resp_tumor_pseudobulk.csv
# Output: ../Result/phoenix_resp_topgenes_tumor_NR_vs_R.tiff
#         ../Result/phoenix_resp_DE_tumor_NR_vs_R.csv  (all genes, Welch t p)
# ============================================================================

library(dplyr)

base_dir   <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
result_dir <- file.path(base_dir, "Result")
N <- 10

pb  <- read.csv(file.path(result_dir, "phoenix_resp_tumor_pseudobulk.csv"),
                row.names = 1, check.names = FALSE)
grp <- data.frame(patient_id = rownames(pb), group = pb$group)
cat(sprintf("NR=%d, R=%d\n", sum(grp$group == "NR"), sum(grp$group == "R")))
gene_cols <- setdiff(names(pb), "group")
matH <- as.matrix(pb[grp$patient_id[grp$group == "NR"], gene_cols])
matL <- as.matrix(pb[grp$patient_id[grp$group == "R"],  gene_cols])

# ---- Vectorized Welch t (NR vs R), all genes at once ----
welch_t <- function(A, B) {             # A = NR, B = R (patients x genes)
  mA <- colMeans(A); mB <- colMeans(B)
  vA <- apply(A, 2, var) / nrow(A); vB <- apply(B, 2, var) / nrow(B)
  se <- sqrt(vA + vB)
  t  <- (mA - mB) / se
  df <- (vA + vB)^2 / (vA^2 / (nrow(A) - 1) + vB^2 / (nrow(B) - 1))
  t[!is.finite(t)] <- 0; df[!is.finite(df)] <- 1
  list(t = t, p = 2 * pt(-abs(t), df))
}

# ---- Unpaired DE: A = R, B = NR ----
mean_A <- colMeans(matL, na.rm = TRUE)
mean_B <- colMeans(matH, na.rm = TRUE)
de <- data.frame(gene = gene_cols, mean_A = mean_A, mean_B = mean_B,
                 log2FC = mean_B - mean_A,          # log1p space, so difference = log fold change
                 stringsAsFactors = FALSE, row.names = NULL)
de$welch_p <- welch_t(matH, matL)$p
de$padj <- p.adjust(de$welch_p, method = "BH")
de <- de[is.finite(de$log2FC), ]

out_de <- de %>% rename(mean_R = mean_A, mean_NR = mean_B) %>%
  mutate(expr_diff = abs(mean_NR - mean_R)) %>% arrange(-expr_diff)
write.csv(out_de, file.path(result_dir, "phoenix_resp_DE_tumor_NR_vs_R.csv"), row.names = FALSE)

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

  # Panel 1: Up in B (Non-responder)
  bp <- barplot(up_norm, col = colB, border = NA,
                ylab = "Expression score",
                main = paste0("Up in ", labelB), col.main = colB,
                cex.main = 0.9, font.main = 2,
                ylim = c(0, 1.15), axes = FALSE)
  axis(2, las = 1, tcl = -0.2, lwd = 0.3, lwd.ticks = 0.3)
  box(bty = "l", lwd = 0.3)
  text(bp, -0.04, labels = up$gene, srt = 45, adj = 1, xpd = TRUE, cex = 0.75)

  # Panel 2: Up in A (Responder)
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

  cat("\nTop up in", labelB, ":\n"); print(up[, c("gene", "mean_A", "mean_B", "log2FC", "welch_p")])
  cat("\nTop up in", labelA, ":\n"); print(dn[, c("gene", "mean_A", "mean_B", "log2FC", "welch_p")])
}

plot_topgenes(de, "Responder", "Non-responder", "#377EB8", "#E41A1C",
              file.path(result_dir, "phoenix_resp_topgenes_tumor_NR_vs_R.tiff"))
cat("Done!\n")
