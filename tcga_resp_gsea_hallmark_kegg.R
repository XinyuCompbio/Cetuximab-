############################################################
# RUN ON: local Mac
#   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
#   Rscript tcga_resp_gsea_hallmark_kegg.R
# GSEA Hallmark + KEGG, tumor-enriched pseudobulk, Non-responder vs Responder
# (TCGA cetuximab patients with response data). From Win4221 gsea_hallmark_v2.R,
# gsea_kegg_v2.R, gsea_NES_barplot.R, kegg_NES_barplot.R, but UNPAIRED:
# ranking = Welch t statistic (NR vs R) after scaling each patient to the
# same total expression, ties broken by effect size.
#   NES > 0 = up in Non-responder (red, right)
#   NES < 0 = up in Responder  (blue, left)
# Input : ../Result/tcga25_tumor_pseudobulk.csv
#         ../tcga_hnsc_cetuximab_patients.csv  (response labels)
# KEGG = current KEGG pathways via clusterProfiler::gseKEGG (needs internet),
#        standard pathway names (e.g. "Focal adhesion"). Human disease
#        pathways (hsa05xxx) are hidden from the plot except cancer (hsa052xx).
# Output: ../Result/tcga_resp_GSEA_hallmark_NR_vs_R.csv  + _NES_barplot.png
#         ../Result/tcga_resp_GSEA_kegg_NR_vs_R.csv      + _NES_barplot.png
############################################################
rm(list = ls()); graphics.off()
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggplot2)

base_dir   <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
result_dir <- file.path(base_dir, "Result")

PCUT   <- 0.05          # raw pvalue filter for the NES plot
N_SIDE <- 10            # top N per side
MAX_CHAR <- 40          # max label length before trimming with "..."
LOWEXP_Q <- 0.25        # drop genes in the bottom 25% of mean expression
DROP_DISEASE <- TRUE    # hide KEGG human disease pathways (keep cancer) in the plot
COL_HIGH <- "#E41A1C"   # same as KM / topgenes
COL_LOW  <- "#377EB8"

# ---- Data ----
pb  <- read.csv(file.path(result_dir, "tcga25_tumor_pseudobulk.csv"),
                row.names = 1, check.names = FALSE)
# ---- Groups straight from the response CSV (no extra step) ----
# R = Complete/Partial Response, NR = Progressive/Stable Disease, others dropped
grp <- read.csv(file.path(base_dir, "tcga_hnsc_cetuximab_patients.csv"))
grp$group <- ifelse(grepl("Complete Response|Partial Response", grp$cetuximab_response), "R",
             ifelse(grepl("Progressive Disease|Stable Disease", grp$cetuximab_response), "NR", NA))
grp <- grp[!is.na(grp$group), ]
grp <- grp[grp$patient_id %in% rownames(pb), ]
gene_cols <- setdiff(names(pb), c("n_tumor_spots", "n_total_spots"))
# ---- Per-patient scaling: every patient to the same total expression ----
# removes per-patient total-expression differences before ranking
tot <- rowSums(pb[, gene_cols])
pb[, gene_cols] <- pb[, gene_cols] / tot * mean(tot)
cat("total per patient before scaling, NR:", round(mean(tot[grp$patient_id[grp$group == "NR"]])),
    " R:", round(mean(tot[grp$patient_id[grp$group == "R"]])), "\n")

# ---- Drop near-zero genes (bottom LOWEXP_Q by mean) ----
# e.g. olfactory receptors sit at ~0.01, pure noise that GSEA picks up
gm <- colMeans(pb[, gene_cols])
gene_cols <- gene_cols[gm > quantile(gm, LOWEXP_Q)]
cat("genes kept after low-expression filter:", length(gene_cols), "\n")

matH <- as.matrix(pb[grp$patient_id[grp$group == "NR"], gene_cols])
matL <- as.matrix(pb[grp$patient_id[grp$group == "R"],  gene_cols])
cat(sprintf("NR=%d, R=%d, genes=%d\n", nrow(matH), nrow(matL), length(gene_cols)))

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
wt <- welch_t(matH, matL); tstat <- wt$t; pval <- wt$p
res <- data.frame(symbol = gene_cols,
                  mean_R = colMeans(matL, na.rm = TRUE),
                  mean_NR = colMeans(matH, na.rm = TRUE),
                  t = tstat, P.Value = pval, stringsAsFactors = FALSE)
res$diff <- res$mean_NR - res$mean_R

gl_df <- res %>% filter(is.finite(t), symbol != "") %>% arrange(desc(t), desc(abs(diff)))
gl <- gl_df$t; names(gl) <- gl_df$symbol
gl <- gl[!duplicated(names(gl))]
cat(sprintf("ranked genes: %d, t range %.2f to %.2f\n", length(gl), min(gl), max(gl)))

# ---- Gene sets, version-safe ----
to_t2g <- function(h) {
  sym <- if ("gene_symbol" %in% names(h)) "gene_symbol" else "human_gene_symbol"
  out <- unique(as.data.frame(h[, c("gs_name", sym)]))
  names(out) <- c("gs_name", "gene_symbol"); out
}
get_hallmark <- function() {
  h <- try(msigdbr(species = "Homo sapiens", collection = "H"), silent = TRUE)
  if (inherits(h, "try-error") || nrow(h) == 0)
    h <- msigdbr(species = "Homo sapiens", category = "H")
  to_t2g(h)
}
# ---- GSEA ----
run_gsea <- function(sets, tag) {
  cat(sprintf("\n==== %s: %d sets, overlap genes %d ====\n", tag,
              length(unique(sets$gs_name)),
              length(intersect(names(gl), unique(sets$gene_symbol)))))
  set.seed(42)
  g <- GSEA(geneList = gl, TERM2GENE = sets,
            minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 1, pAdjustMethod = "BH",
            eps = 0, nPermSimple = 10000, seed = TRUE)
  df <- as.data.frame(g)
  csv <- file.path(result_dir, paste0("tcga_resp_GSEA_", tag, "_NR_vs_R.csv"))
  write.csv(df, csv, row.names = FALSE)
  cat("Saved:", csv, "\n")
  cat("\n-- top 10 up in Non-responder --\n")
  print(head(df[order(-df$NES), c("ID","setSize","NES","pvalue","p.adjust")], 10), row.names = FALSE)
  cat("\n-- top 10 up in Responder --\n")
  print(head(df[order(df$NES), c("ID","setSize","NES","pvalue","p.adjust")], 10), row.names = FALSE)
  df
}

run_kegg <- function() {
  map <- suppressWarnings(bitr(names(gl), fromType = "SYMBOL", toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db))
  map <- map[!duplicated(map$SYMBOL) & !duplicated(map$ENTREZID), ]
  gk <- gl[map$SYMBOL]; names(gk) <- map$ENTREZID
  gk <- sort(gk, decreasing = TRUE)
  cat(sprintf("\n==== kegg: %d genes mapped to Entrez ====\n", length(gk)))
  set.seed(42)
  g <- gseKEGG(geneList = gk, organism = "hsa",
               minGSSize = 10, maxGSSize = 500,
               pvalueCutoff = 1, pAdjustMethod = "BH",
               eps = 0, nPermSimple = 10000, seed = TRUE)
  g <- setReadable(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  df <- as.data.frame(g)
  df$Description <- sub(" - Homo sapiens \\(human\\)$", "", df$Description)
  csv <- file.path(result_dir, "tcga_resp_GSEA_kegg_NR_vs_R.csv")
  write.csv(df, csv, row.names = FALSE)
  cat("Saved:", csv, "\n")
  if (DROP_DISEASE) df <- df[!(grepl("^hsa05", df$ID) & !grepl("^hsa052", df$ID)), ]
  cat("\n-- top 10 up in Non-responder --\n")
  print(head(df[order(-df$NES), c("Description","setSize","NES","pvalue","p.adjust")], 10), row.names = FALSE)
  cat("\n-- top 10 up in Responder --\n")
  print(head(df[order(df$NES), c("Description","setSize","NES","pvalue","p.adjust")], 10), row.names = FALSE)
  df
}

# ---- NES diverging barplot (template style) ----
ACRONYMS <- c("DNA","RNA","G2M","E2F","MYC","KRAS","TNFA","NFKB","TGF","MTORC1",
              "MTOR","P53","UV","DN","IL2","STAT5","IL6","JAK","STAT3","PI3K",
              "AKT","WNT","ROS","ATP","EMT","V1","V2","IL","TNF",
              "ECM","TCA","EGF","EGFR","ERK","RAS","PKC","PLC","NF","MHC","TCR","BCR","CoA","NADPH","HCM","ARVC","MAPK","STAT","VEGF","GNRH","PPAR","ABC","SNARE","MRNA")

plot_nes <- function(df, title, lab_size) {
  d0 <- df %>% filter(!is.na(NES), !is.na(pvalue), pvalue < PCUT)
  cat(sprintf("\n%s significant sets (raw pvalue < %.2f): %d\n", title, PCUT, nrow(d0)))
  if (nrow(d0) == 0) { cat("nothing to plot\n"); return(invisible(NULL)) }

  up   <- d0 %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(N_SIDE)
  down <- d0 %>% filter(NES < 0) %>% arrange(NES)       %>% head(N_SIDE)
  cat("up in Non-responder:", nrow(up), " up in Responder:", nrow(down), "\n")

  d <- bind_rows(up, down)
  if (grepl("^HALLMARK_", d$ID[1])) {
    d$label <- gsub("_", " ", sub("^HALLMARK_", "", d$ID))
    d$label <- paste0(substr(d$label, 1, 1), tolower(substr(d$label, 2, nchar(d$label))))
    for (a in ACRONYMS)
      d$label <- gsub(paste0("\\b", a, "\\b"), a, d$label, ignore.case = TRUE)
  } else {
    d$label <- d$Description          # KEGG: official pathway name
  }
  # MEDICUS names are long, trim so labels stay inside the panel
  d$label <- ifelse(nchar(d$label) > MAX_CHAR,
                    paste0(substr(d$label, 1, MAX_CHAR - 3), "..."), d$label)
  d$side <- ifelse(d$NES > 0, "high", "low")
  d <- d %>% arrange(NES)
  d$ID <- factor(d$ID, levels = d$ID)
  xmax <- max(abs(d$NES)) * 1.15

  g <- ggplot(d, aes(x = NES, y = ID, fill = side)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
    geom_text(aes(x = ifelse(NES > 0, -0.06, 0.06), label = label,
                  hjust = ifelse(NES > 0, 1, 0)),
              size = lab_size / .pt, family = "Arial", colour = "black") +
    scale_fill_manual(values = c(high = COL_HIGH, low = COL_LOW), guide = "none") +
    scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0)) +
    labs(x = "NES", y = NULL, title = title) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          axis.text.x  = element_text(colour = "black", size = 7),
          axis.title.x = element_text(size = 7),
          axis.line.x  = element_line(linewidth = 0.3),
          axis.ticks.x = element_line(linewidth = 0.3),
          plot.title   = element_text(size = 8, family = "Arial", face = "bold",
                                      hjust = 0.5, margin = margin(b = 1)),
          plot.margin  = margin(4, 2, 2, 2)) +
    annotate("text", x = -xmax * 0.98, y = nrow(d) + 0.9, label = "Up in Responder",
             hjust = 0, vjust = 0, colour = COL_LOW, fontface = "bold",
             size = 7 / .pt, family = "Arial") +
    annotate("text", x =  xmax * 0.98, y = nrow(d) + 0.9, label = "Up in Non-responder",
             hjust = 1, vjust = 0, colour = COL_HIGH, fontface = "bold",
             size = 7 / .pt, family = "Arial") +
    coord_cartesian(ylim = c(1, nrow(d)), clip = "off")

  out <- file.path(result_dir, paste0("tcga_resp_GSEA_", tolower(title), "_NES_barplot.png"))
  ggsave(out, g, width = 3, height = 3, units = "in", dpi = 600)
  cat("Saved:", out, "\n")
}

hm <- run_gsea(get_hallmark(), "hallmark")
plot_nes(hm, "Hallmark", 7)

kg <- run_kegg()
plot_nes(kg, "KEGG", 6)

cat("\nDone.\n")
