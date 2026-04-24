## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  comment = "##"
)
options(width=60)
library(SimiCviz)


## ----echo = FALSE-------------------------------------------------------------
simic_full <- readRDS(system.file("extdata", file.path("simic_full.rds"), 
                                  package = "SimiCviz"))


## ----eval = FALSE, results='hold'---------------------------------------------
# library(SimiCviz)
# 
# # Load entire SimiCPipeline run automatically
# simic_full <- load_SimiCPipeline(
#   project_dir = "path/to/simic_run",
#   run_name = "example1",
#   lambda1 = "0.01",
#   lambda2 = "0.001"
# )
# 
# # Set display names and colors for visualization (Part 3)
# simic_full <- setLabelNames(
#   simic_full,
#   label_names = c("control", "PD-L1", "DAC", "Combination"),
#   colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0")
# )


## -----------------------------------------------------------------------------
simic_full


## -----------------------------------------------------------------------------

# Load weights from pickle
weights_file <- system.file("extdata", 
  file.path("outputSimic/example1_simic_weights.pickle"), 
  package = "SimiCviz")

simic_weights <- read_weights_pickle(weights_file)
simic_weights[[1]][, 1:6]



## -----------------------------------------------------------------------------
# Load GRN weights (method agnostic)

weight_path <- system.file("extdata", "example_weights.csv", 
                           package = "SimiCviz")

# Read as data.frame
weights_df <- read_weights_csv(weight_path)
head(weights_df)

# If your method uses different column names you need to rename them


## -----------------------------------------------------------------------------
# Load from CSV (recommended format: columns 'cell', 'label')
cell_labels_path <- system.file("extdata",
  file.path("inputFiles", "treatment_annotation.csv"), 
  package = "SimiCviz")

cell_labels <- load_cell_labels(cell_labels_path, header = TRUE, sep = ",")
head(cell_labels)


## ----eval=FALSE---------------------------------------------------------------
# # From vector (will generate cell_1, cell_2, ... names)
# cell_labels <- c(0, 0, 1, 1, 2, 2)  # Must match order of AUC rows
# 
# # From named vector
# cell_labels <- c(cell_A = 0, cell_B = 0, cell_C = 1)
# 
# # From data.frame
# cell_labels <- data.frame(cell = c("cell_A", "cell_B"), label = c(0, 1))


## -----------------------------------------------------------------------------
# Load from multiple formats
expression_mat_path <- system.file("extdata",
  file.path("inputFiles", "example1_expression.pickle"),
  package = "SimiCviz")

# Will auto-detect format and load
expression_mat <- load_expression_matrix(expression_mat_path)
print(class(expression_mat))
print(dim(expression_mat))


## -----------------------------------------------------------------------------
# From .pickle SimiC files

# Extract Adjusted R squared from SimiC outputs
out <- read_pickle(weights_file)
adjusted_r_squared <- out$adjusted_r_squared

viz_obj_simic <- SimiCvizExperiment(
  weights = simic_weights,
  auc = NULL,  # Will compute this in the next section later but can be loaded as well
  cell_labels = cell_labels,
  label_names = c("control","PD-L1","DAC","Combination"),
  colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"),
  meta=list(adjusted_r_squared=adjusted_r_squared))

viz_obj_simic



## -----------------------------------------------------------------------------
# From SimiC input files
# adj_r2_threshold: For SimiCPipeline style outputs it will look in metadata of viz_obj_simic, for "simic@meta$adjusted_r_squared" and filter out targets with lower R²)
# n_cores: Number of workers if backend = 'multissession' or 
#          Number of cores if 'multiprocess'

viz_obj_simic <- calculate_activity_scores(
              viz_obj_simic,
              expression = expression_mat_path,
              adj_r2_threshold = 0.7, # For SimiC style outputs 
              sort_by="expression", # Rank targets by expression or weight
              select_top_k = NULL,  # Use all targets (or limit to top K)
              percent_of_target = 1.0,  
              n_cores = 2,
              backend = "multisession",
              verbose = TRUE
            )


# Access computed scores
auc_scores <- viz_obj_simic@auc$collected
head(auc_scores[, 1:5])


## -----------------------------------------------------------------------------
# Initialize processor with weights and expression
AS_processor <- AUCProcessor(
  weights = weights_df,
  expression = expression_mat_path,  # or matrix/data.frame
  cell_labels = cell_labels,
  qc_type = "adj_p_val",
  qc_threshold = 0.05,  # Filter targets above this threshold
  n_cores = 2,
  backend = "multisession"
)

# Compute with custom parameters
AS_processor <- compute_auc(
  AS_processor,
  sort_by = "expression",
  select_top_k = NULL,            # Use all targets (or limit to top K)
  percent_of_target = 1.0,        # Use all targets (or subset %)
  verbose = TRUE
)


## -----------------------------------------------------------------------------
# Extract results in wide format (default)
auc_wide <- get_auc(AS_processor, format = "wide")
head(auc_wide)


## -----------------------------------------------------------------------------
# Extract results long format
auc_long <- get_auc(AS_processor, format = "long") # Long format
head(auc_long)


## ----eval = F-----------------------------------------------------------------
# weights_df <- read.csv("path/to/your/weights.csv")
# weights_df_filtered <- weights_df[weights_df$p_value < 0.01,]
# # From CSV files
# viz_obj <- SimiCvizExperiment(
#   weights = weights_df_filtered,
#   auc = NULL,  # Will compute this inthe next section later but can be loaded as well
#   cell_labels = cell_labels,
#   label_names = c("control","PD-L1","DAC","Combination"),
#   colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"),
#   meta=list() # Anything you want to store in a list format
#   )
# 
# viz_obj
# viz_obj <- calculate_activity_scores(
#               viz_obj,
#               expression = expression_mat_path,
#               adj_r2_threshold = 0.7, # For SimiC style outputs
#               sort_by="expression", # Rank targets by expression or weight
#               select_top_k = NULL,  # Use all targets (or limit to top K)
#               percent_of_target = 1.0,
#               n_cores = 2,
#               backend = "multisession",
#               verbose = TRUE
#             )


## ----eval = FALSE-------------------------------------------------------------
# # SimiC: Filter by adjusted R² (goodness of fit)
# 
# processor_simic <- AUCProcessor(
#   weights = simic_weights,
#   expression = expr_mat,
#   cell_labels = cell_labels,
#   adj_r2_list = adjusted_r_squared, # a list length as simic_weights
#   qc_type = "adj_r2",
#   qc_threshold = 0.7  # Keep targets with R² ≥ 0.7
# )
# 
# # SCENIC / Pando: Filter by adjusted p-value
# processor_scenic <- AUCProcessor(
#   weights = weights_df,
#   expression = expr_mat,
#   cell_labels = cell_labels,
#   qc_type = "p_value",
#   qc_threshold = 0.05,  # Keep targets with adj_p_val ≤ 0.05
#   n_cores = 4,
#   backend = "multisession"
#   )
# 
# # Compute with the same data, different parameters
# processor_scenic <- compute_auc(processor_scenic, sort_by = "weight")


## -----------------------------------------------------------------------------
# Recall above examples
 simic_full # Complete SimiCpipeline output
 viz_obj_simic # SimiCPipeline weights -> `calculate_activity_scores`


## -----------------------------------------------------------------------------
# Create SimiCvizExperiment
simic <- SimiCvizExperiment(weights = simic_weights,
                             auc = auc_wide,
                             cell_labels = cell_labels,
                             label_names = c("control","PD-L1","DAC","Combination"),
                             colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
simic


## ----echo = FALSE-------------------------------------------------------------
plot_dir <- file.path(getwd(),"SimiCviz_output","plots")

## ----eval = FALSE-------------------------------------------------------------
# plot_dir <- file.path(getwd(),"SimiCviz_output")
# dir.create(plot_dir,recursive = TRUE)


## -----------------------------------------------------------------------------
# Plot distribution of adjusted R² values across targets
# Extract Adjusted R squared from SimiC outputs
out <- read_pickle(weights_file)
adjusted_r_squared <- out$adjusted_r_squared
plot_r2_distribution(adjusted_r_squared, simic, grid = c(2, 2), 
                     save = FALSE, out_dir = plot_dir)


## -----------------------------------------------------------------------------
# Select targets by adjusted R2
unselected_targets <- list()
selected_targets <- list()
lab_keys <- names(simic@label_names)
for (lab in lab_keys){
    # Save selected for plotting
    selected_targets[[lab]] <- simic@target_ids[which(adjusted_r_squared[[lab]] >= 0.7)]
    # Save unselected for reporting
    label <- simic@label_names[[lab]]
    unselected_targets[[label]] <- simic@target_ids[which(adjusted_r_squared[[lab]] < 0.7)]
}
print("Number of unselected targets per label:")
print(sapply(unselected_targets, length)) 



## ----fig.height=10, fig.width=10----------------------------------------------
plot_tf_weights(
  simic,
  tf_names = simic@tf_ids[1:4],
  top_n = 25,
  allowed_targets = selected_targets,  # Filter by R² if desired
  grid = c(2, 2),
  save = FALSE,
  out_dir = plot_dir,
  filename = "TF_weights_barplot.pdf"
)


## ----fig.height=5, fig.width=5------------------------------------------------
plot_target_weights(
  simic,
  target_names = simic@target_ids[1:4],
  labels = c("control", "Combination"),
  grid = c(2, 2),
  save = FALSE,
  out_dir = plot_dir,
  filename = "Target_weights_barplot.pdf"
)


## ----include=FALSE------------------------------------------------------------
all_tfs_barplots <- plot_tf_weights(
                          simic,
                          top_n = 25, 
                          grid = NULL,
                          allowed_targets = selected_targets)



## ----eval = F-----------------------------------------------------------------
# all_tfs_barplots <- plot_tf_weights(
#                           simic,
#                           top_n = 25,
#                           grid = NULL,
#                           allowed_targets = selected_targets)

## -----------------------------------------------------------------------------
all_tfs_barplots[[1]]


## -----------------------------------------------------------------------------
network <- get_tf_network(simic_full, "Tet2", r2_threshold = 0.7)
print(head(network))

plot_tf_network_heatmap(simic_full, "Tet2", 
                        save = FALSE, 
                        top_n = 15,
                        r2_threshold = 0.7,
                        show_values = T, 
                        cmap = c("purple","white","yellow"))



## ----collapse=TRUE, results='hold'--------------------------------------------
dis_score <- calculate_dissimilarity(simic)
top_tfs <- rownames(dis_score)


## -----------------------------------------------------------------------------
plot_dissimilarity_heatmap(simic, 
                           top_n = 5, 
                           cmap = "viridis",
                           save = FALSE)


## -----------------------------------------------------------------------------
metadata <- read.csv(system.file("extdata/metadata.csv", 
  package = "SimiCviz"))
 
# Build cell groups from metadata (e.g. Seurat clusters, cell types, etc.)
cell_groups  <- lapply(unique(metadata$final_annotation_functional), 
                       function(celltype) {
  cell_labels$cell[metadata$final_annotation_functional == celltype]
})
names(cell_groups) <- unique(metadata$final_annotation_functional)

dissim_grouped <- calculate_dissimilarity(simic, labels = c(0,2),
                                          cell_groups = cell_groups)

# For all labels
plot_dissimilarity_heatmap(simic,
                            cell_groups = cell_groups, 
                            top_n = 5,
                            cmap=c("magma"),
                            save = FALSE, out_dir = plot_dir,
                            filename = "dissimilarity_heatmap_grouped.pdf")


## -----------------------------------------------------------------------------
# For labels 0,2
plot_dissimilarity_heatmap(simic,
                            labels = c(0,2),
                            cell_groups = cell_groups, 
                            top_n = 5, 
                            sort_by = "Proliferating.cells",
                            cmap=c("red", "white", "blue"),
                            save = FALSE)


## ----fig.height=10, fig.width=10----------------------------------------------
# Plot distributions for top TFs
plot_auc_distributions(
  simic,
  tf_names = top_tfs[1:4],
  fill = TRUE,
  alpha = 0.6,
  bw_adjust = 1/8,
  rug = TRUE,
  save = FALSE,
  out_dir = plot_dir,
  filename = "AUC_distributions.pdf",
  grid = c(2, 2)
)

## ----fig.height=5, fig.width=10-----------------------------------------------
# Plot top 4 TFs density distributions
plot_auc_distributions(simic,
                       labels = c(0,3),
                       tf_names = top_tfs[1:2],
                       fill = FALSE,
                       bw_adjust = 0.5,
                       rug = FALSE,
                       out_dir = plot_dir,
                       filename="AUC_distributions_notfilled_multipage.pdf",
                       save = FALSE,
                       grid = c(1,2))


## ----fig.height=10, fig.width=10----------------------------------------------
plot_auc_cumulative(
  simic,
  tf_names = top_tfs[1:4],
  rug = TRUE,
  grid = c(2, 2),
  include_table = TRUE,
  save = FALSE,
  out_dir = plot_dir
)


## ----eval=TRUE----------------------------------------------------------------
ecdf_metrics <- calculate_ecdf_auc(simic, tf_names = simic@tf_ids[1:4])
head(ecdf_metrics)


## -----------------------------------------------------------------------------
plot_auc_heatmap(simic, top_n = 20)


## ----eval=TRUE----------------------------------------------------------------
summary_plot <- plot_auc_summary_statistics(simic)


## -----------------------------------------------------------------------------
sessionInfo()

