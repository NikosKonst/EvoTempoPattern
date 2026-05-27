# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering,
#that leads to the generation of the species whole optic lobe data in the cloeon_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers cloeon_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

larval1.data <- Read10X(data.dir = "/Cdipterum_larva1_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

larval2.data <- Read10X(data.dir = "/Cdipterum_larva2_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

larval3.data <- Read10X(data.dir = "/Cdipterum_larva3_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_197152_correspondence_GE_subsample.tsv", header = F)

larval1.data@Dimnames[[1]] <- mgi$V1
larval2.data@Dimnames[[1]] <- mgi$V1
larval3.data@Dimnames[[1]] <- mgi$V1

#create Seurat object with standard filters
larv1 <-CreateSeuratObject(counts =larval1.data, min.cells = 3, min.features = 200,project="lib1")
larv2 <-CreateSeuratObject(counts =larval2.data, min.cells = 3, min.features = 200,project="lib2")
larv3 <-CreateSeuratObject(counts =larval3.data, min.cells = 3, min.features = 200,project="lib3")
rm(larval1.data,larval2.data,larval3.data)

# decide thresholds by visualizations, first plot vln plots and hist_features, 
# decide thresholds and visualize them 

#### QC Figure S4 ####
larv_list <- list(
  larv1 = larv1,
  larv2 = larv2,
  larv3 = larv3
)

# functions for each plot

# Histogram of nFeature_RNA with vertical reference line
hist_feat <- function(obj) {
  # Histogram for nFeature_RNA with a red dashed line at x = 
  ggplot2::qplot(
    obj$nFeature_RNA,
    geom = "histogram",
    bins = 101,
    fill = I("white"),
    col  = I("black")
  ) +
    ggplot2::geom_vline(
      xintercept = 800,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 800, # or different
               linetype = "dashed",
               color = "red", size = 1) +
    guides(fill = "none")
}


# list with the plots i want to generate
plot_types <- list(
  
  vlnplot        = vlnplot,
  hist_feat      = hist_feat
)

output_dir <- "/plots/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_type_names <- names(plot_types)
larv_names      <- names(larv_list)

for (obj_name in larv_names) {
  obj <- larv_list[[obj_name]]
  
  for (plot_name in plot_type_names) {
    plot_fun <- plot_types[[plot_name]]
    p <- plot_fun(obj)
    
    ggsave(
      filename = file.path(output_dir, paste0("plot_", plot_name, "_", obj_name, ".jpg")),
      plot     = p,
      units    = "cm", width = 40, height = 30, bg = "white"
    )
  }
}

#### subset 
larv1 <- subset(larv1, subset = nFeature_RNA > 800)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 800)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 800)
larv3 <- NormalizeData(larv3, normalization.method = "LogNormalize", scale.factor = 10000)
larv3 <- FindVariableFeatures(larv3, selection.method = "vst", nfeatures = 2000)


reference.list <- c(larv1, larv2, larv3)
opticlobes.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:150)

larvaOL.integrated <- IntegrateData(anchorset = opticlobes.anchors, dims = 1:150)

rm(reference.list, opticlobes.anchors)

larvaOL.integrated <- ScaleData(object = larvaOL.integrated, verbose = FALSE)
larvaOL.integrated <- RunPCA(object = larvaOL.integrated, npcs = 150, verbose = FALSE)
larvaOL.integrated <- RunUMAP(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindNeighbors(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindClusters(object = larvaOL.integrated, resolution = 4)

saveRDS(larvaOL.integrated, "/cloeon_OL.rds")
larvaOL.integrated <- readRDS("/cloeon_OL.rds")
larvaOL.integrated <-FindClusters(object = larvaOL.integrated, resolution = 8)

DimPlot(larvaOL.integrated, label = TRUE, raster=F, shuffle = T) + NoLegend()
DimPlot(larvaOL.integrated, split.by = "orig.ident",label = TRUE) + NoLegend()

#### choοse and save medulla NBs ####
# find medulla NB cluster based on conserved markers
FeaturePlot(larvaOL.integrated, features = "rna_dpn") #NBs

nw_object <- c(2)
NbS_all <- subset(larvaOL.integrated_4, idents = nw_object)
NbS_all <- FindVariableFeatures(NbS_all, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS_all, verbose = FALSE)
NbS_all <- RunPCA(NbS_all, npcs = 50)
NbS_all <- RunUMAP(NbS_all, dims = 1:50)
NbS_all <- FindNeighbors(NbS_all, dims = 1:50)
NbS_all <- FindClusters(NbS_all, resolution = 1.5)

# dpn + elav+
FeaturePlot(NbS_all, features = c("rna_dpn","rna_CLODIP-2-CD15777"))
DimPlot(NbS_all, label = T)

#second round of subsetting
nw_object <- c(0,1,2,4)
NbS <- subset(NbS_all, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)

ScaleData(object = NbS_all, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)

DimPlot(NbS, label = T)

DefaultAssay(NbS) <- "RNA"
rna_assay <- NbS[["RNA"]]
# Get current feature (gene) names
rnames <- rownames(rna_assay)
# Your rename map
rename <- c("CLODIP-2-CD15498" = "B-H1",
            "CLODIP-2-CD15283" = "erm",
            "CLODIP-2-CD02704" = "Oaz",
            "CLODIP-2-CD07852" = "ey",
            "CLODIP-2-CD00329" = "D",
            "CLODIP-2-CD08459" = "hth",
            "CLODIP-2-CD14302" = "dmrt99B"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(NbS[["RNA"]]) <- rnames

saveRDS(NbS, "/cloeon_medulla_NBs.rds")

# remove low quality clusters (that do not integrate properly among the three libraries) from musca_OL.rds for visualization
#no need for cloeon

