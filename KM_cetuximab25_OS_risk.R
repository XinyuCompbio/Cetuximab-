## ─────────────────────────────────────────────────────────────
## RUN ON: local Mac (Elements drive mounted)
##   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
##   Rscript split_cetuximab_OS_median.R   # makes the 35-patient input CSV first
##   Rscript KM_cetuximab25_OS_risk.R
## KM for the 25 cetuximab patients with a DeepSpotM slide (.h5ad.gz) in
## /Volumes/Elements/DeepspotM/TCGA, OS, re-split at the median OS of these 25
## Short (<= median) = High risk, Long = Low risk
## Panel 2.5 x 2 inch, 600 dpi, TIFF, Arial
## Input:  ../Result/cetuximab_35_OS_groups.csv
## Output: ../Result/cetuximab_25_OS_groups.csv
##         ../Result/KM_cetuximab25_OS_risk.tiff
## ─────────────────────────────────────────────────────────────
library(survival)
library(survminer)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)

# ----- Paths -----
base_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                     error = function(e) getwd())
proj_dir <- dirname(base_dir)
out_dir  <- file.path(proj_dir, "Result")

df <- read.csv(file.path(out_dir, "cetuximab_35_OS_groups.csv"))

# keep patients with a slide in the DeepSpotM TCGA folder
slide_dir <- "/Volumes/Elements/DeepspotM/TCGA"
have <- unique(substr(list.files(slide_dir, pattern = "\\.h5ad\\.gz$"), 1, 12))
df <- df[df$patient_id %in% have, ]
cat(sprintf("Patients with slides: %d\n", nrow(df)))

# re-split at median OS of these patients
med <- median(df$OS_time)
df$OS_group <- ifelse(df$OS_time > med, "Long_OS", "Short_OS")
df$OS <- df$OS_time / 30.44                     # days -> months
df$group <- factor(ifelse(df$OS_group == "Short_OS", "High", "Low"),
                   levels = c("High", "Low"))
cat(sprintf("Median OS = %.0f days; High=%d, Low=%d\n",
            med, sum(df$group == "High"), sum(df$group == "Low")))
write.csv(df[order(df$OS_time), ], file.path(out_dir, "cetuximab_25_OS_groups.csv"),
          row.names = FALSE)

# ===== POSITION CONTROLS (months axis) =====
ylab_size   <- 7
x_max  <- 80
x_by   <- 20
pval_x <- x_max * 0.05
pval_y <- 0.12
plot_margin_t <- 2; plot_margin_r <- 2; plot_margin_b <- 3; plot_margin_l <- 12
leg_x <- x_max * 0.62; leg_y <- 1.00; leg_dy <- 0.10
leg_seg <- x_max * 0.09; leg_gap <- x_max * 0.03; leg_size <- 7
nrisk_y <- 1.35; nrisk_x <- -x_max * 0.025
high_y <- 0.72; low_y <- 0.28
high_label_x <- -x_max * 0.12; low_label_x <- -x_max * 0.12   # left-aligned start; more negative = further left
table_xnum_mt <- 1
time_y <- -1.15
small_size <- 6
table_margin_t <- 8; table_margin_r <- 2; table_margin_b <- 9; table_margin_l <- 12
table_panel_h <- 0.40
COL_HIGH <- "#E41A1C"; COL_LOW <- "#377EB8"
# ====================================================

n_at_risk <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$n.risk
  if (length(s) == 0) 0 else s
}

save_km_panel <- function(dat, tcol, ecol, ylab_text, out_path) {
  dat$time  <- dat[[tcol]]
  dat$event <- dat[[ecol]]
  fit <- survfit(Surv(time, event) ~ group, data = dat)

  p_km <- ggsurvplot(fit, data = dat,
                     pval = TRUE, pval.size = 7 / .pt,
                     pval.coord = c(pval_x, pval_y),
                     risk.table = FALSE,
                     palette = c(COL_HIGH, COL_LOW),
                     legend = "none",
                     legend.labs = c("High", "Low"),
                     xlab = "", ylab = ylab_text,
                     break.x.by = x_by,
                     xlim = c(0, x_max),
                     censor.size = 1,
                     size = 0.4,
                     ggtheme = theme_classic() +
                       theme(text = element_text(family = "Arial", size = 6),
                             axis.title = element_text(size = ylab_size),
                             axis.text = element_text(size = small_size),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.line = element_line(linewidth = 0.3),
                             axis.ticks = element_line(linewidth = 0.3),
                             axis.ticks.length = unit(0.04, "cm"),
                             legend.background = element_blank(),
                             legend.margin = margin(0, 0, 0, 0),
                             plot.title = element_blank(),
                             plot.margin = margin(plot_margin_t, plot_margin_r,
                                                  plot_margin_b, plot_margin_l)))

  time_breaks <- seq(0, x_max, by = x_by)

  p_km$plot <- p_km$plot +
    scale_x_continuous(breaks = time_breaks, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = c(0, x_max), clip = "off") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank())

  # ----- Manual legend -----
  leg_txt <- leg_x + leg_seg + leg_gap
  y_hi    <- leg_y - leg_dy
  y_lo    <- leg_y - 2 * leg_dy

  p_km$plot <- p_km$plot +
    annotate("text", x = leg_x, y = leg_y, label = "Risk", hjust = 0,
             size = leg_size / .pt, family = "Arial") +
    annotate("segment", x = leg_x, xend = leg_x + leg_seg, y = y_hi, yend = y_hi,
             colour = COL_HIGH, linewidth = 0.4) +
    annotate("text", x = leg_txt, y = y_hi, label = "High", hjust = 0,
             size = leg_size / .pt, family = "Arial") +
    annotate("segment", x = leg_x, xend = leg_x + leg_seg, y = y_lo, yend = y_lo,
             colour = COL_LOW, linewidth = 0.4) +
    annotate("text", x = leg_txt, y = y_lo, label = "Low", hjust = 0,
             size = leg_size / .pt, family = "Arial") +
    theme(legend.position = "none")

  # ----- Risk table -----
  risk_data <- data.frame()
  for (grp in c("High", "Low")) {
    sf <- survfit(Surv(time, event) ~ 1, data = dat[dat$group == grp, ])
    for (tb in time_breaks) {
      risk_data <- rbind(risk_data, data.frame(time = tb, group = grp, n = n_at_risk(sf, tb)))
    }
  }
  risk_data$y <- ifelse(risk_data$group == "High", high_y, low_y)

  p_table <- ggplot(risk_data, aes(x = time, y = y, label = n)) +
    geom_text(hjust = 0.5, size = small_size / .pt, family = "Arial", color = "black") +
    annotate("text", x = high_label_x, y = high_y, label = "High",
             color = COL_HIGH, size = small_size / .pt, family = "Arial", hjust = 0) +
    annotate("text", x = low_label_x, y = low_y, label = "Low",
             color = COL_LOW, size = small_size / .pt, family = "Arial", hjust = 0) +
    scale_x_continuous(breaks = time_breaks, expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, 1), clip = "off") +
    annotate("text", x = nrisk_x, y = nrisk_y, label = "Number at risk",
             hjust = 0, size = small_size / .pt, family = "Arial") +
    annotate("text", x = x_max / 2, y = time_y, label = "Time (months)",
             hjust = 0.5, size = 7 / .pt, family = "Arial") +
    labs(x = NULL, y = "Risk", title = NULL) +
    theme_void() +
    theme(text = element_text(family = "Arial", size = small_size),
          plot.title = element_blank(),
          axis.text.x = element_text(size = small_size, margin = margin(t = table_xnum_mt)),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 6, angle = 90),
          axis.line.x = element_line(linewidth = 0.3),
          axis.line.y = element_line(linewidth = 0.3),
          axis.ticks.x = element_line(linewidth = 0.3),
          axis.ticks.length = unit(0.04, "cm"),
          legend.position = "none",
          plot.margin = margin(table_margin_t, table_margin_r,
                               table_margin_b, table_margin_l))

  # ----- Combine -----
  g1 <- ggplotGrob(p_km$plot)
  g2 <- ggplotGrob(p_table)
  g1$widths <- unit.pmax(g1$widths, g2$widths)
  g2$widths <- unit.pmax(g1$widths, g2$widths)
  panel_row <- g2$layout$t[g2$layout$name == "panel"]
  g2$heights[panel_row] <- unit(table_panel_h, "cm")
  combined <- arrangeGrob(g1, g2, ncol = 1,
                          heights = unit.c(unit(1, "null"), sum(g2$heights)))
  ggsave(out_path, combined, width = 2.5, height = 2, units = "in",
         dpi = 600, device = "tiff", compression = "lzw")
  cat(sprintf("Saved: %s\n", basename(out_path)))
}

cat("OS:\n")
save_km_panel(df, "OS", "OS_event", "Overall Survival",
              file.path(out_dir, "KM_cetuximab25_OS_risk.tiff"))
cat("\nDone.\n")
