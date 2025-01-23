run_deseq_analysis <- function(counts, metadata, min_count,
                               min_samples, purity_var, comparison,
                               additional_covariates = c("APC_CELLS", "NON_APC_CELLS", 
                                                         "Endothelial", "Fibroblasts"),
                               group_var = "group") {
  
  # Validate inputs
  if(!all(additional_covariates %in% colnames(metadata))) {
    stop("Not all additional covariates found in metadata")
  }
  if(!purity_var %in% colnames(metadata)) {
    stop(sprintf("Purity variable '%s' not found in metadata", purity_var))
  }
  
  # Construct design formula
  covariates_formula <- paste(additional_covariates, collapse = " + ")
  design_formula <- as.formula(
    sprintf("~ %s + %s + %s:%s + %s",
            group_var, purity_var, group_var, purity_var, covariates_formula)
  )
  
  message("Using design formula: ", deparse(design_formula))
  
  # Create DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = design_formula
  )
  
  # Filter low counts
  keep <- rowSums(counts(dds) > min_count) > min_samples
  message("Keeping ", sum(keep), " out of ", length(keep), " genes")
  dds <- dds[keep,]
  
  # Run DESeq
  dds <- DESeq(dds, quiet = TRUE)
  
  # Verify the comparison exists
  available_comparisons <- resultsNames(dds)
  if(!comparison %in% available_comparisons) {
    stop(sprintf("Requested comparison '%s' not found. Available comparisons: %s",
                 comparison, paste(available_comparisons, collapse = ", ")))
  }
  
  # Get results
  res <- data.frame(results(dds, name = comparison)) %>%
    rownames_to_column("gene")
  
  return(list(
    results = res,
    design = design_formula
  ))
}

#' Create and Preprocess Seurat Object
#' @description Creates a Seurat object and performs standard preprocessing steps
#' @param counts Raw count matrix
#' @param metadata Cell metadata (optional)
#' @param project_name Project name for the Seurat object (default: "project")
#' @param min_cells Minimum number of cells expressing a gene (default: 10)
#' @param min_features Minimum number of features per cell (default: 100)
#' @param norm_method Normalization method (default: "LogNormalize")
#' @param scale_factor Scale factor for normalization (default: 10000)
#' @param var_features_method Method for variable feature selection (default: "vst")
#' @param n_var_features Number of variable features to select (default: 2000)
#' @param scale_all_features Logical, whether to scale all features or only variable features (default: TRUE)
#' @return Preprocessed Seurat object
#' @import Seurat
#' @export
#' 
seurat_wrapper <- function(
    counts,
    metadata = NULL,
    project_name = "project",
    min_cells = 10,
    min_features = 100,
    norm_method = "LogNormalize",
    scale_factor = 10000,
    var_features_method = "vst",
    n_var_features = 2000,
    scale_all_features = TRUE
) {
  # Input validation
  tryCatch({
    if (!is.matrix(counts) && !is.data.frame(counts)) {
      stop("counts must be a matrix or data frame")
    }
    
    if (!is.null(metadata) && nrow(metadata) != ncol(counts)) {
      stop("Number of rows in metadata must match number of columns in counts matrix")
    }
    
    # Create Seurat object
    message("Creating Seurat object...")
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      meta.data = metadata,
      project = project_name,
      min.cells = min_cells,
      min.features = min_features
    )
    
    # Print initial dimensions
    message("Initial dimensions: ", 
            paste(dim(seurat_obj), collapse = " x "))
    
    # Normalize data
    message("Normalizing data...")
    seurat_obj <- NormalizeData(
      seurat_obj,
      normalization.method = norm_method,
      scale.factor = scale_factor,
      verbose = FALSE
    )
    
    # Find variable features
    message("Finding variable features...")
    seurat_obj <- FindVariableFeatures(
      seurat_obj,
      selection.method = var_features_method,
      nfeatures = n_var_features,
      verbose = FALSE
    )
    
    # Scale data
    message("Scaling data...")
    features_to_scale <- if (scale_all_features) {
      message("Scaling all features...")
      rownames(seurat_obj)
    } else {
      message("Scaling only variable features...")
      VariableFeatures(seurat_obj)
    }
    
    seurat_obj <- ScaleData(
      seurat_obj,
      features = features_to_scale,
      verbose = FALSE
    )
    
    # Print summary statistics
    message("\nPreprocessing complete:")
    message("Final dimensions: ", 
            paste(dim(seurat_obj), collapse = " x "))
    message("Number of variable features: ", 
            length(VariableFeatures(seurat_obj)))
    
    return(seurat_obj)
    
  }, error = function(e) {
    stop("Error in Seurat preprocessing: ", e$message)
  })
}
