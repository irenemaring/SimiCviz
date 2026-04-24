# SimiCviz NEWS

## SimiCviz 0.99.0 (2026-04-24)

### New package release!
- Initial public release of SimiCviz.
- Includes the S4 containers:
  - SimiCvizExperiment for coordinated GRN weights, activity scores, labels, and metadata.
  - AUCProcessor for activity-score computation workflows.
- Includes import utilities for:
  - SimiCPipeline outputs.
  - Generic CSV-based GRN and AUC inputs.
  - Optional Python-backed pickle/h5ad reading via reticulate interoperability.
- Includes activity-score functionalities:
  - TF activity score calculation from expression and weights.
  - Dissimilarity score calculation.
- Includes visualization functions for:
  - TF-target weight barplots.
  - TF network heatmaps.
  - AUC distributions, cumulative curves, and summary statistics.
  - Dissimilarity heatmaps across labels/cell groups.
  - Adjusted R2 distribution plots.

### Interoperability
- Python-based pickle/h5ad reading is optional.
- Interoperability support is included for workflows that commonly use SimiCPipeline python outputs.
- Outputs from alternative algorithms for gene regulatory network are also compativle in standard csv format input.