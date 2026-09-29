## ─────────────────────────────────────────────────────────────
## RUN ON: local Mac
##   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
##   Rscript tcga25_gigatime_High_vs_Low.R
## GigaTIME virtual mIF (21 channels), 25 TCGA cetuximab patients,
## High risk vs Low risk. Per-slide values from Source_Data/TCGA_gigatime.
## Stats: Wilcoxon rank-sum (High vs Low), BH, Spearman rho with OS days
##        (continuous, no split).
## Input : Source_Data/TCGA_gigatime/gigatime_activation_density.csv
##         Source_Data/TCGA_gigatime/gigatime_mean_intensity.csv
##         ../Result/cetuximab_25_OS_groups.csv
## Output: ../Result/tcga25_gigatime_High_vs_Low.csv  (both metrics, all markers)
##         ../Result/tcga25_gigatime_density_boxplot.tiff (activation density)
## ─────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)

base_dir <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics"
gt_dir   <- file.path(base_dir, "Source_Data/TCGA_gigatime")
res_dir  <- file.path(base_dir, "cetuximab_improve/DeepspotM/Result")
COL_HIGH <- "#E41A1C"; COL_LOW <- "#377EB8"

grp <- read.csv(file.path(res_dir, "cetuximab_25_OS_groups.csv"))
grp <- grp[, c("patient_id", "group", "OS_time", "OS_event")]
grp$group <- factor(grp$group, levels = c("High", "Low"))

compare <- function(file, metric) {
  d <- read.csv(file, check.names = FALSE)
  d <- merge(grp, d, by = "patient_id")
  markers <- setdiff(names(d), names(grp))
  out <- do.call(rbind, lapply(markers, function(m) {
    h <- d[d$group == "High", m]; l <- d[d$group == "Low", m]
    data.frame(metric = metric, marker = m,
               mean_High = mean(h), mean_Low = mean(l),
               diff_High_minus_Low = mean(h) - mean(l),
               wilcox_p = wilcox.test(h, l)$p.value,
               spearman_rho_OS = cor(d[[m]], d$OS_time, method = "spearman"),
               spearman_p_OS = suppressWarnings(cor.test(d[[m]], d$OS_time,
                                                         method = "spearman")$p.value))
  }))
  out$padj <- p.adjust(out$wilcox_p, method = "BH")
  list(stats = out[order(out$wilcox_p), ], data = d, markers = markers)
}

dens <- compare(file.path(gt_dir, "gigatime_activation_density.csv"), "activation_density")
intn <- compare(file.path(gt_dir, "gigatime_mean_intensity.csv"),     "mean_intensity")
cat(sprintf("patients matched: %d\n", nrow(dens$data)))

all_stats <- rbind(dens$stats, intn$stats)
write.csv(all_stats, file.path(res_dir, "tcga25_gigatime_High_vs_Low.csv"), row.names = FALSE)
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
  scale_fill_manual(values = c(High = COL_HIGH, Low = COL_LOW)) +
  labs(x = "Risk", y = "GigaTIME activation density") +
  theme_classic(base_size = 7) +
  theme(legend.position = "none", strip.background = element_blank())

out <- file.path(res_dir, "tcga25_gigatime_density_boxplot.tiff")
ggsave(out, g, width = 7, height = 4.5, units = "in", dpi = 600,
       device = "tiff", compression = "lzw")
cat("\nSaved:", out, "\n")
