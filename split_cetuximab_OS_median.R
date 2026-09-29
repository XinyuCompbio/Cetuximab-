# Split 35 TCGA-HNSC cetuximab patients into 2 groups by OS time (median split)
library(survival)

base <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics"
pts  <- read.csv(file.path(base, "cetuximab_improve/DeepspotM/tcga_hnsc_cetuximab_patients.csv"))
clin <- read.delim(file.path(base, "Source_Data/Metadata/TCGA-HNSC.clinical.tsv"),
                   check.names = FALSE, na.strings = c("", "'--"))

clin <- clin[clin$submitter_id %in% pts$patient_id,
             c("submitter_id", "vital_status.demographic",
               "days_to_death.demographic", "days_to_last_follow_up.diagnoses")]
colnames(clin) <- c("patient_id", "vital_status", "days_to_death", "days_to_last_fu")

# one row per patient (take max follow-up across duplicate rows)
clin <- do.call(rbind, lapply(split(clin, clin$patient_id), function(x) {
  data.frame(patient_id = x$patient_id[1],
             vital_status = x$vital_status[1],
             days_to_death = suppressWarnings(max(x$days_to_death, na.rm = TRUE)),
             days_to_last_fu = suppressWarnings(max(x$days_to_last_fu, na.rm = TRUE)))
}))
clin[sapply(clin, is.infinite)] <- NA

clin$OS_event <- as.integer(clin$vital_status == "Dead")
clin$OS_time  <- ifelse(clin$OS_event == 1 & !is.na(clin$days_to_death),
                        clin$days_to_death, clin$days_to_last_fu)

df <- merge(pts, clin, by = "patient_id")
med <- median(df$OS_time)
df$OS_group <- ifelse(df$OS_time > med, "Long_OS", "Short_OS")
df <- df[order(df$OS_time), ]

cat("Median OS time (days):", med, "\n")
print(table(df$OS_group, Event = df$OS_event))
print(survdiff(Surv(OS_time, OS_event) ~ OS_group, data = df))

write.csv(df, file.path(base, "cetuximab_improve/DeepspotM/Result/cetuximab_35_OS_groups.csv"),
          row.names = FALSE)
