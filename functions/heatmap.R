library(smoother)
library(splus2R)
library(stringr)
library(zoo) # for moving average smoothing
library(viridis)
library(SingleCellExperiment)
library(slingshot)
library(pheatmap)
library(dplyr)
library(scales)

#### heatmap with customized order, adapted from Cachero, Mittleton et al. 2026 biorxiv, https://doi.org/10.1101/2025.07.16.664682  ####
ordered_heat_max <- function(object, 
                             features, 
                             window_size = 0.05, 
                             alpha_val = 0.5,
                             span.value = 21, 
                             byMax = TRUE,
                             outfolder = "",
                             filename = "") {
  
  # Expression matrix (adjust if using Seurat or monocle3 object)
  matrix <- Nb@assays$RNA@data 
  genes_keep <- features
  pseudotimes <- LPS_ori_tbl_umap$pseudotime
  ordered.pseudo <- pseudotimes[order(pseudotimes, decreasing = F)]
  
  # Filter genes and order cells by pseudotime
  matrix.filter <- matrix[genes_keep, , drop = FALSE]
  matrix.ordered <- matrix.filter[, order(pseudotimes, decreasing = F)]
  
  if (byMax) {
    ## --- Order by peak position ---
    raw_peaks <- t(apply(matrix.ordered, 1, FUN = peaks, span = span.value, strict = TRUE, endbehavior = 1))
    
    peak.positions <- sapply(1:nrow(matrix.ordered), function(i) {
      peak.indices <- which(raw_peaks[i, ])
      if (length(peak.indices) > 0) {
        return(peak.indices[which.max(matrix.ordered[i, peak.indices])])
      } else {
        return(which.max(matrix.ordered[i, ]))
      }
    })
    
    ordering <- order(peak.positions)
  } else {
    ## --- Use custom order as given in 'features' ---
    ordering <- match(features, rownames(matrix.ordered))
  }
  
  # Apply ordering
  matrix.ordered.by.peak <- matrix.ordered[ordering, , drop = FALSE]
  
  # Smooth expression
  matrix.smooth <- t(apply(matrix.ordered.by.peak, 1, smth.gaussian, 
                           window = window_size, alpha = alpha_val, tails = TRUE))
  
  # Heatmap color breaks
  ## Heatmap breaks and plotting
  #col.br <- c(seq(0, 2.5, length.out = 50), seq(2.6, 4, length.out = 50)) # virilis is made with that
  col.br <- seq(0, 4, length.out = 100)
  # Save heatmap
  png(filename = file.path(outfolder, paste0(filename, ".png")), units = "cm", width = 25, height = 25, bg = "white", res = 600)
  pheatmap(matrix.smooth, 
           #color = viridis(100, option = "plasma"),  # or "plasma", "inferno", "cividis" or nth so it is the blue/yellow classic
           breaks = col.br,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           scale = "row",
           show_colnames = FALSE,
           treeheight_col = 0)
  dev.off()
}
