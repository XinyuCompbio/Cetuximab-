## ─────────────────────────────────────────────────────────────
## RUN ON: local Mac
##   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
##   Rscript tcga_resp_gigatime_NR_vs_R.R
## GigaTIME virtual mIF (21 channels), TCGA cetuximab patients,
## Non-responder (NR: Progressive/Stable Disease) vs Responder (R: Complete/Partial Response). Per-slide values from Source_Data/TCGA_gigatime.
## Stats: Wilcoxon rank-sum (NR vs R), BH.
## Input : Source_Data/TCGA_gigatime/gigatime_activation_density.csv
##         Source_Data/TCGA_gigatime/gigatime_mean_intensity.csv
##         ../tcga_hnsc_cetuximab_patients.csv  (response labels)
## Output: ../Result/tcga_resp_gigatime_NR_vs_R.csv  (both metrics, all markers)
##         ../Result/tcga_resp_gigatime_density_boxplot.tiff (activation density)
## ─────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)

base_dir <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics"
gt_dir   <- file.path(base_dir, "Source_Data/TCGA_gigatime")
res_dir  <- file.path(base_dir, "cetuximab_improve/DeepspotM/Result")
COL_HIGH <- "#E41A1C"; COL_LOW <- "#377EB8"   # NR red, R blue

# R = Complete/Partial Response, NR = Progressive/Stable Disease, others dropped
grp <- read.csv(file.path(base_dir, "cetuximab_improve/DeepspotM/tcga_hnsc_cetuximab_patients.csv"))
grp$group <- ifelse(grepl("Complete Response|Partial Response", grp$cetuximab_response), "R",
             ifelse(grepl("Progressive Disease|Stable Disease", grp$cetuximab_response), "NR", NA))
grp <- grp[!is.na(grp$group), c("patient_id", "group")]
grp$group <- factor(grp$group, levels = c("NR", "R"))

compare <- function(file, metric) {
  d <- read.csv(file, check.names = FALSE)
  d <- merge(grp, d, by = "patient_id")
  markers <- setdiff(names(d), names(grp))
  out <- do.call(rbind, lapply(markers, function(m) {
    h <- d[d$group == "NR", m]; l <- d[d$group == "R", m]
    data.frame(metric = metric, marker = m,
               mean_NR = mean(h), mean_R = mean(l),
               diff_NR_minus_R = mean(h) - mean(l),
               wilcox_p = wilcox.test(h, l)$p.value)
  }))
  out$padj <- p.adjust(out$wilcox_p, method = "BH")
  list(stats = out[order(out$wilcox_p), ], data = d, markers = markers)
}

dens <- compare(file.path(gt_dir, "gigatime_activation_density.csv"), "activation_density")
intn <- compare(file.path(gt_dir, "gigatime_mean_intensity.csv"),     "mean_intensity")
cat(sprintf("patients matched: %d\n", nrow(dens$data)))

all_stats <- rbind(dens$stats, intn$stats)
write.csv(all_stats, file.path(res_dir, "tcga_resp_gigatime_NR_vs_R.csv"), row.names = FALSE)
cat("\n-- activation density --\n");  print(dens$stats, row.names = FALSE, digits = 3)
cat("\n-- mean intensity --\n");      print(intn$stats, row.names = FALSE, digits = 3)

# ---- Boxplot, activation density, markers ordered by p ----
long <- dens$data %>%
  pivot_longer(all_of(dens$markers), names_to = "marker", values_to = "value")
lab <- setNames(sprintf("%s\np=%.3f", dens$stats$marker, dens$stats$wilcox_p), dens$stats$marker)
long$marker <- factor(long$marker, levels = dens$stats$marker, labels = lab[dens$stats$marker])

g <- ggplot(long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.6) +
  facet_wrap(~ marker, scales = "free_y", ncol = 7) +
  scale_fill_manual(values = c(NR = COL_HIGH, R = COL_LOW)) +
  labs(x = "Cetuximab response", y = "GigaTIME activation density") +
  theme_classic(base_size = 7) +
  theme(legend.position = "none", strip.background = element_blank())

out <- file.path(res_dir, "tcga_resp_gigatime_density_boxplot.tiff")
ggsave(out, g, width = 7, height = 4.5, units = "in", dpi = 600,
       device = "tiff", compression = "lzw")
cat("\nSaved:", out, "\n")
