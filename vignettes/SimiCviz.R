## ----setup, include=FALSE---------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  comment = "##"
)
options(width=60)
# install.packages("magick")
# library(magick)
library(SimiCviz)

## ----collapse=TRUE, results='hold'------------------------
library(SimiCviz)
simic <- load_SimiCPipeline(project_dir = "/home/workdir/SimiCviz/inst/extdata",
                               run_name = "example1",
                               lambda1 = 0.01,
                               lambda2 = 0.001)

message("\nSetting lables and colors")
simic <- setLabelNames(simic, label_names  = c("control","PD-L1","DAC","Combo"),
              colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
simic


## ---------------------------------------------------------
# Weights output dictionary from SimiCPipeline
weights_file <- system.file("extdata",file.path("outputSimic/matrices/example1",
                  "example1_L1_0.01_L2_0.001_simic_matrices_filtered_BIC.pickle"), 
                            package = "SimiCviz")
# This function extracts only the weights matrices from the pickle file, 
weights_results <- read_weights_pickle(weights_file)
head(weights_results[[1]][,1:6])

# If you want to access all the content use `read_pickle()` function
out <- read_pickle(weights_file)
str(out,1)

## ---------------------------------------------------------
auc_collect_file <- system.file("extdata",file.path("outputSimic/matrices/example1"
                                ,"example1_L1_0.01_L2_0.001_wAUC_matrices_filtered_BIC_collected.csv"),
                                package = "SimiCviz")
auc_collect <- read.csv(auc_collect_file, header=TRUE, row.names = 1)
head(auc_collect[,1:5])

# In New SimiCpipeline the input cell labels has "cell" and "label" columns.
cell_labels_path <- system.file("extdata",
                                file.path("inputFiles","treatment_annotation.csv"), 
                                package = "SimiCviz")
cell_labels_new = load_cell_labels(cell_labels_path, header = TRUE, sep = ",")
head(cell_labels_new)

## ---------------------------------------------------------
auc_file <- system.file("extdata", file.path("outputSimic/matrices/example1"
                        ,"example1_L1_0.01_L2_0.001_wAUC_matrices_filtered_BIC.pickle"), 
                        package = "SimiCviz")
auc_pkl <- read_auc_pickle(auc_file)

# In old SimiCpipeline the input cell labels is just a vector of labels, 
# we need to load it and convert it to a dataframe mapping "cell" and "label" columns.
cell_labels_path <- system.file("extdata", 
                               file.path("inputFiles","treatment_annotation.txt"), 
                               package = "SimiCviz")
cell_labels_old = read.csv(cell_labels_path, header = FALSE)$V1
head(cell_labels_old)

## ---------------------------------------------------------
simic2 <- SimiCvizExperiment(weights = weights_results,
                             auc = auc_collect,
                             cell_labels = cell_labels_new,
                             label_names = c("control","PD-L1","DAC","Combo"),
                             colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
simic2

## ----eval = FALSE-----------------------------------------
#  # Equivalent (not run)
#  simic2_bis <- SimiCvizExperiment(weights = weights_results,
#                               auc = auc_pkl,
#                               cell_labels = cell_labels_old,
#                               label_names = c("control","PD-L1","DAC","Combo"),
#                               colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
#  simic2_bis

## ----collapse=TRUE,results='hold'-------------------------
weight_path <-  system.file("extdata","example_weights.csv", package = "SimiCviz")
auc_path <-  system.file("extdata","example_auc.csv", package = "SimiCviz")
cell_labels_path <- system.file("extdata", 
                                file.path("inputFiles","treatment_annotation.csv"),
                                package = "SimiCviz")
# Loading using the paths
simic3 <- load_from_csv(weights_file = weight_path, 
                        auc_file = auc_path, 
                        cell_labels_file = cell_labels_path,
                        meta = list(lambda1 = 0.01, lambda2= 0.001))

simic3 <- setLabelNames(simic3, label_names  = c("control","PD-L1","DAC","Combo"),
                        colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
simic3

## ----eval = FALSE-----------------------------------------
#  # Equivalent (not run)
#  weight_df <- read_weights_csv(weight_path)
#  head(weight_df)
#  auc_df <- read_auc_csv(auc_path)
#  head(auc_df)
#  
#  # Loading using the dataframe (equivalent)
#  simic3_bis <- SimiCvizExperiment(weights = weight_df,
#                                   auc = auc_df,
#                                   cell_labels = cell_labels_new, # Need to have "cell" and "label"
#                                   meta = list(lambda1 = 0.01, lambda2= 0.001),
#                                   label_names = c("control","PD-L1","DAC","Combo"),
#                                   colors = c("#e0e0e0", "#a8c8ff", "#ffb6b6", "#c1a9e0"))
#  simic3_bis

## ---------------------------------------------------------
plot_dir <- file.path(getwd(),"SimiCviz_output","plots")
dir.create(plot_dir,recursive = TRUE)

## ---------------------------------------------------------
# load_SimiCPipeline()` function will already load the `adjusted_r_squared` values in the meta slot.
adjusted_r_squared <- simic@meta$adjusted_r_squared
plot_r2_distribution(adjusted_r_squared, simic, grid=c(2,2), width= 16,save =TRUE, out_dir = plot_dir)

## ----eval = FALSE-----------------------------------------
#  # For manual loading you can extract the adjusted R2 values from the weights pickle file.
#  out <- read_pickle(weights_file)
#  adjusted_r_squared <- out$adjusted_r_squared
#  plot_r2_distribution(adjusted_r_squared, simic2_bis, grid=c(2,2), width= 16,save =TRUE, out_dir = plot_dir)

## ----fig.height=8, fig.width=10,echo =FALSE, message = FALSE----
plot_r2_distribution(adjusted_r_squared, simic, grid=c(2,2))

## ---------------------------------------------------------
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

## ---------------------------------------------------------
 plot_tf_weights(simic,
                tf_names        = simic@tf_ids[1:4],
                top_n           = 25,
                allowed_targets = selected_targets,
                grid            = c(4,1),
                save            = TRUE,
                out_dir = plot_dir, 
                filename = "TF_weights_barplot.pdf"
                ) 


## ----fig.height=8, fig.width=10, echo = FALSE, message = FALSE----
 plot_tf_weights(simic,
                tf_names        = simic@tf_ids[1:2],
                top_n           = 25,
                allowed_targets = selected_targets,
                grid            = c(2,1),
                save            = FALSE
                ) 


## ---------------------------------------------------------
plot_target_weights(simic, 
                target_names = simic@target_ids[1:4],
                labels   = c("control","DAC"),
                save = TRUE, 
                width= 18, # Adjust width and height as needed for grid layout
                height = 10, 
                grid     = c(1, 2),
                out_dir = plot_dir, 
                filename = "Target_weights_barplot.pdf")


## ----fig.height=10, fig.width=15, echo = FALSE, message = FALSE----
plot_target_weights(simic, 
                    target_names = simic@target_ids[1:4],
                    labels   = c("control","DAC"),
                    save = FALSE, 
                    width= 18, height = 10,
                    grid     = c(2, 2))

## ----eval = FALSE-----------------------------------------
#  plot_tf_weights(simic,
#                  save = TRUE,
#                  width= 18, height = 10,
#                  grid     = c(3, 2),
#                  out_dir = plot_dir,
#                  filename = "All_tf_weights_barplot.pdf")
#  plot_target_weights(simic,
#                  save = TRUE,
#                  width= 18, height = 10,
#                  grid     = c(3, 2),
#                  out_dir = plot_dir,
#                  filename = "All_target_weights_barplot.pdf")

## ----fig.height=8, fig.width=10, eval = FALSE-------------
#  plot_network(simic, top_n = 4)

## ----message=FALSE----------------------------------------
# Plot top 4 TFs density distributions
plot_auc_distributions(simic,
                       tf_names = simic@tf_ids[1:4],
                       fill = TRUE,
                       alpha = .6,
                       bw_adjust = 1/8,
                       rug = TRUE,
                       save = TRUE,
                       out_dir = plot_dir,
                       filename = "AUC_distributions.pdf",
                       grid = c(2, 2))


## ----echo = FALSE ,fig.height=10, fig.width=15------------
# Plot top 4 TFs density distributions
plot_auc_distributions(simic,
                       tf_names = simic@tf_ids[1:4],
                       fill = TRUE,
                       alpha = 0.6,
                       bw_adjust = 1/8,
                       rug = TRUE,
                       save = FALSE,
                       grid = c(2, 2))

## ----message=FALSE ,fig.height=10, fig.width=15-----------
# Plot top 6 TFs density distributions
plot_auc_distributions(simic,
                       labels = c(0,2,3),
                       tf_names = simic@tf_ids[1:6],
                       fill = TRUE,
                       alpha = 0.9,
                       bw_adjust = 0.8,
                       rug = TRUE,
                       out_dir = plot_dir,
                       filename="AUC_distributions_filled.pdf",
                       save = TRUE,
                       grid = c(2, 3))

## ----echo = FALSE ,fig.height=10, fig.width=15------------
# Plot top 6 TFs density distributions
plot_auc_distributions(simic,
                       labels = c(0,2,3),
                       tf_names = simic@tf_ids[1:6],
                       fill = TRUE,
                       alpha = 0.9,
                       bw_adjust = 0.8,
                       rug = TRUE,
                       save = FALSE,
                       grid = c(2, 3))

## ---------------------------------------------------------
# Plot top 4 TFs density distributions
plot_auc_distributions(simic,
                       labels = c(0,3),
                       tf_names = simic@tf_ids[1:2],
                       fill = FALSE,
                       alpha = 0.9,
                       bw_adjust = 0.5,
                       rug = TRUE,
                       out_dir = plot_dir,
                       filename="AUC_distributions_notfilled_multipage.pdf",
                       save = TRUE,
                       grid = c(1,2))

## ----echo = FALSE ,fig.height=10, fig.width=15------------
# Plot top 4 TFs density distributions
plot_auc_distributions(simic,
                       labels = c(0,3),
                       tf_names = simic@tf_ids[1:2],
                       fill = FALSE,
                       alpha = 0.9,
                       bw_adjust = 0.5,
                       rug = TRUE,
                       save = FALSE,
                       grid = c(2, 2))

## ----message=FALSE, warning = FALSE-----------------------
plot_auc_cumulative(simic,
                    tf_names = simic@tf_ids[1:2],
                    rug = FALSE,
                    grid = c(2, 2),
                    save = TRUE, 
                    include_table = TRUE,
                    width = 12, height = NULL,
                    out_dir = plot_dir)

## ----echo = FALSE, fig.height=8, fig.width=15-------------
plot_auc_cumulative(simic,
                    tf_names = simic@tf_ids[1:2],
                    rug = FALSE,
                    grid = c(1, 2),
                    save = FALSE)

## ----eval = FALSE-----------------------------------------
#  auc_file <- system.file("extdata", file.path("outputSimic/matrices/example1"
#                          ,"example1_L1_0.01_L2_0.001_wAUC_matrices_filtered_BIC.pickle"),
#                              package = "SimiCviz")
#  auc_pkl <- read_auc_pickle(auc_file)
#  
#  cell_labels_path <- system.file("extdata",file.path("inputFiles","treatment_annotation.csv"), package = "SimiCviz")
#  cell_labels = load_cell_labels(cell_labels_path)
#  auc_df <- auc_list_to_df(auc_pkl, cell_labels_df = cell_labels )
#  head(auc_df[,1:5])
#  
#  save_collected_auc(auc_df, file = "./collected_auc.csv", overwrite = FALSE)

## ----eval = FALSE-----------------------------------------
#  ecdf_metrics <- calculate_ecdf_auc(simic, tf_names = simic@tf_ids[1:6])
#  head(ecdf_metrics)

## ----eval = FALSE-----------------------------------------
#  plot_auc_summary_statistics(simic)

## ----fig.height=8, fig.width=10, eval=FALSE---------------
#  plot_auc_heatmap(simic, top_n = 20)

## ----eval=FALSE-------------------------------------------
#  plot_auc_heatmap(simic)

## ----eval=FALSE-------------------------------------------
#  # getwd()
#  # out_dir <- "/home/workdir/SimiCviz/SimiCviz_output"
#  # dir.create(file.path(out_dir))
#  # export_SimiCviz_all(simic, out_dir = out_dir, prefix = "SimiC_example")

