# Spatial Analysis of Cetuximab Resistance in HNSCC

## Background

Cetuximab is a monoclonal antibody targeting the epidermal growth factor receptor (EGFR) and is used to treat head and neck squamous cell carcinoma (HNSCC). However, treatment outcomes vary substantially, and many patients exhibit either primary or acquired resistance. The biological mechanisms underlying cetuximab resistance remain incompletely understood, limiting the development of effective biomarkers and combination therapies.

Available cetuximab-response datasets primarily contain bulk molecular measurements and provide limited information about the spatial organization of tumor, stromal, and immune cells. To address this limitation, we used deep-learning models to reconstruct spatial molecular features from routinely collected hematoxylin and eosin (H&E) pathology images.

## Objective

This project investigated the molecular, cellular, and spatial mechanisms associated with cetuximab resistance in HNSCC. Resistance-associated features were then used to identify potential therapeutic targets and prioritize rational combination therapies that may improve response to cetuximab.

## Study Design

Patients were classified as cetuximab responders or non-responders using available treatment-response information. Deep-learning models were applied to pretreatment H&E pathology images to generate virtual spatial molecular profiles.

The models included:

- **Phoenix**, which predicted spatially resolved gene expression at single-cell resolution.
- **DeepSpot-M**, which generated transcriptome-wide virtual spatial gene-expression profiles.
- **GigaTIME**, which predicted virtual multiplex immunofluorescence (mIF) profiles of tumor and immune markers.

## Spatial Molecular Profiling

Virtual spatial molecular profiles were generated from H&E whole-slide images. These profiles characterized the distribution and organization of tumor, immune, and stromal compartments within the tumor microenvironment.

## Responder–Non-Responder Comparison

Inferred molecular and spatial features were compared between cetuximab responders and non-responders. The analysis examined:

- Differential gene and protein-marker expression
- EGFR and downstream signaling activity
- Alternative and bypass signaling pathways
- Immune-cell composition and functional states
- Tumor–stroma and tumor–immune interactions
- Spatial co-localization and exclusion patterns
- Cellular neighborhoods associated with treatment resistance
- Intratumoral heterogeneity

## Resistance-Mechanism Analysis

Resistance-associated molecular and spatial features were integrated to identify biological programs linked to reduced cetuximab sensitivity. The analysis focused on:

- Reactivation of EGFR downstream signaling
- Activation of alternative receptor tyrosine kinases
- Epithelial–mesenchymal transition
- Altered cell-proliferation and survival pathways
- Immune suppression and immune exclusion
- Stromal remodeling
- Spatially localized resistant tumor-cell populations

## Survival Analysis

TCGA-HNSC cases with available pathology images, clinical outcomes, and treatment information were used for a complementary clinical analysis.

Resistance-associated features were combined into a risk score, and patients were divided into high- and low-risk groups. Kaplan–Meier analysis and Cox proportional-hazards models were used to examine associations between the identified resistance programs and patient survival.

## Combination-Therapy Discovery

Resistance-associated pathways, cellular states, and spatial interactions were mapped to potentially actionable therapeutic targets. Candidate drugs and combination strategies were prioritized according to:

- Their ability to inhibit resistance-associated pathways
- Their biological complementarity with EGFR inhibition
- Their relevance to the identified tumor and immune-cell states
- Existing evidence in HNSCC or related cancers
- Their potential to reverse immune suppression or spatial immune exclusion

This analysis produced a set of mechanism-based hypotheses for combining cetuximab with other targeted or immune-modulating therapies.

## Project Outputs

This project generated:

1. Virtual spatial molecular profiles of cetuximab-responsive and non-responsive HNSCC.
2. Candidate molecular, cellular, and spatial mechanisms of cetuximab resistance.
3. Resistance-associated features linked to patient survival.
4. Potentially actionable targets associated with resistant tumor states.
5. A prioritized set of rational cetuximab combination-therapy hypotheses.

## Long-Term Significance

This work demonstrated a computational framework for using routine H&E pathology images to investigate cetuximab resistance. The identified mechanisms and therapeutic targets provide a foundation for developing more effective, mechanism-based combination therapies for patients with HNSCC.
