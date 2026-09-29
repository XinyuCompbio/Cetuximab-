## ─────────────────────────────────────────────────────────────
## RUN ON: local Mac
##   cd /Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM/Code
##   Rscript make_response_groups.R
## Cetuximab response groups for patients with a DeepSpotM slide.
##   R  (Responder)     = Complete Response or Partial Response
##   NR (Non-responder) = Clinical Progressive Disease or Stable Disease
##   Not available / Unknown / Not applicable only -> dropped
## Input : ../tcga_hnsc_cetuximab_patients.csv
##         ../Result/cetuximab_25_OS_groups.csv (the 25 with slides)
## Output: ../Result/cetuximab_response_groups.csv
## ─────────────────────────────────────────────────────────────
base_dir <- "/Users/xinyucui/Desktop/Cancer_Bioinformatics/cetuximab_improve/DeepspotM"
resp <- read.csv(file.path(base_dir, "tcga_hnsc_cetuximab_patients.csv"))
os   <- read.csv(file.path(base_dir, "Result/cetuximab_25_OS_groups.csv"))
resp <- resp[resp$patient_id %in% os$patient_id, ]

resp$group <- ifelse(grepl("Complete Response|Partial Response", resp$cetuximab_response), "R",
              ifelse(grepl("Progressive Disease|Stable Disease", resp$cetuximab_response), "NR", NA))
resp <- resp[!is.na(resp$group), ]
resp <- merge(resp, os[, c("patient_id", "OS_time", "OS_event")], by = "patient_id")

print(resp[order(resp$group), ], row.names = FALSE)
print(table(resp$group))
write.csv(resp, file.path(base_dir, "Result/cetuximab_response_groups.csv"), row.names = FALSE)
cat("Saved: Result/cetuximab_response_groups.csv\n")
