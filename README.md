# SimiCviz

SimiCviz — Visualization tools for SimiC and SimiCPipeline outputs.

A lightweight R/Bioconductor-oriented package to import, summarize and plot single cell gene regulatory network outputs (weights, TF activity/AUC matrices, networks and related figures) produced by SimiC or other GRN inference tools.

## Status
- Current: Install from GitHub (development).
- Planned: Submit to Bioconductor for formal distribution.

## Installation

Install from GitHub:
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("irenemaring/SimiCviz")
```

(When available on Bioconductor, installation instructions will be updated.)

## Quick start

```r
library(SimiCviz)

# Load example data included with the package
simic <- load_SimiCPipeline(project_dir = system.file("extdata", package = "SimiCviz"),
                            run_name = "example1",
                            lambda1 = 0.01, lambda2 = 0.001)

# Set labels and plot an example
simic <- setLabelNames(simic, label_names = c("control","PD-L1","DAC","Combo"),
                       colors = c("#e0e0e0","#a8c8ff","#ffb6b6","#c1a9e0"))
# Alternativ data formats

# Weight file
weight_path <-  system.file("extdata","example_weights.csv", 
                            package = "SimiCviz")
weight_df <- read_weights_csv(weight_path)
head(weight_df)

# Activity Score File
auc_path <-  system.file("extdata","example_auc.csv", package = "SimiCviz")
auc_df <- read_auc_csv(auc_path)
head(auc_df)

simic3 <- SimiCvizExperiment(weights = weight_df,
                              auc = auc_df,
                              tf_ids = unique(auc_df$tf),
                              target_ids = unique(weight_df$target),
                              meta = list(lambda1 = 0.01, lambda2= 0.001),
                              label_names = c("control","PD-L1","DAC","Combo"),
                              colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))

plot_tf_weights(simic, tf_names = simic@tf_ids[1:4], top_n    = 30)

plot_auc_distributions(simic3, tf_names = simic3@tf_ids[1:4], 
                       labels = c("control","Combo")
                       fill = TRUE, 
                       alpha = .6,
                       bw_adjust = 1/8,
                       rug = TRUE,
                       grid = c(2,2))
```

## Main features
- Read SimiCPipeline outputs (pickle/CSV) and standardize into a SimiCvizExperiment object.
- Visualize TF weights and target weights (barplots), regulatory networks, and AUC/activity distributions.
- Integration with Seurat for UMAP visualization.
- Export organized data and publication-ready plots.
- Support for importing results from other methods via CSV.

## Graphical abstract / Logo
(placeholder — add an image here when available)

![Graphical abstract](docs/figures/graphical_abstract.png)
<!-- Replace the path above with your generated image -->

## Documentation
Vignettes and function documentation are included in the package. After installing, open the main vignette:
```r
browseVignettes("SimiCviz")
```
## Contributing
Contributions, issues and feature requests are welcome. Please open issues or pull requests on the GitHub repository.

## Contact
Irene Marín-Goñi — imarin.4@alumni.unav.es

## License
MIT
