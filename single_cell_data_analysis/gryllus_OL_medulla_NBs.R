# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering,
#that leads to the generation of the species whole optic lobe data in the gryllus_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers gryllus_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#### stage 1 ####
larval1.data <- Read10X(data.dir = "/Gbimaculatus_stage11_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval2.data <- Read10X(data.dir = "/Gbimaculatus_stage12_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval3.data <- Read10X(data.dir = "/Gbimaculatus_stage13_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_GRYBI_correspondence_GE_subsample.tsv", header = F)

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
  # Histogram for nFeature_RNA with a red dashed line at x = ..
  ggplot2::qplot(
    obj$nFeature_RNA,
    geom = "histogram",
    bins = 101,
    fill = I("white"),
    col  = I("black")
  ) +
    ggplot2::geom_vline(
      xintercept = 850,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 850, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 850)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 850)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 850)
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

saveRDS(larvaOL.integrated, "/stage1.rds")

#### stage3 ####
larval1.data <- Read10X(data.dir = "/Gbimaculatus_stage31_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval2.data <- Read10X(data.dir = "/Gbimaculatus_stage32_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval3.data <- Read10X(data.dir = "/Gbimaculatus_stage33_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_GRYBI_correspondence_GE_subsample.tsv", header = F)

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
  # Histogram for nFeature_RNA with a red dashed line at x = library1:1400, lbrary2:1300, library3: 1400
  ggplot2::qplot(
    obj$nFeature_RNA,
    geom = "histogram",
    bins = 101,
    fill = I("white"),
    col  = I("black")
  ) +
    ggplot2::geom_vline(
      xintercept = 1400,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 1400, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 1400)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 1300)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 1400)
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

saveRDS(larvaOL.integrated, "/stage3.rds")

#### stage5 ####
larval1.data <- Read10X(data.dir = "/Gbimaculatus_stage51_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval2.data <- Read10X(data.dir = "/Gbimaculatus_stage52_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")
larval3.data <- Read10X(data.dir = "/Gbimaculatus_stage53_GeneExt_merged_subsample100m_force20k/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_GRYBI_correspondence_GE_subsample.tsv", header = F)

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
  # Histogram for nFeature_RNA with a red dashed line at x = library1:550, library2:500, library3:550
  ggplot2::qplot(
    obj$nFeature_RNA,
    geom = "histogram",
    bins = 101,
    fill = I("white"),
    col  = I("black")
  ) +
    ggplot2::geom_vline(
      xintercept = 550,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 550, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 550)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 500)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 550)
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

saveRDS(larvaOL.integrated, "/stage5.rds")

#### all_together to extract medulla NBs ####
stage1 <- readRDS("/stage1.rds")
stage3 <- readRDS("/stage3.rds")
stage5 <- readRDS("/stage5.rds")

stage1@meta.data$Stages <- "Stage1"
stage3@meta.data$Stages <- "Stage3"
stage5@meta.data$Stages <- "Stage5" 

# first merge and then integrate layers with new seurat IntegrateLayers
merged_obj <- merge(x = stage1, y = list(stage3, stage5))
DefaultAssay(merged_obj)<-"RNA"
merged_obj <-NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(object = merged_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged_obj <-ScaleData(object = merged_obj, verbose = FALSE)
merged_obj <- RunPCA(object = merged_obj, npcs = 150, verbose = FALSE)
merged_obj <- IntegrateLayers(object = merged_obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.cca", dims = 1:150, reduction.name = "umap.cca")
merged_obj <- FindNeighbors(merged_obj,reduction = "integrated.cca", dims = 1:150)
merged_obj <- FindClusters(merged_obj, resolution = 4, cluster.name = "cca_clusters")

# merging together stages 1,3,5
saveRDS(merged_obj, "/gryllus_OL.rds")

# find medulla NB cluster based on conserved markers
# NBS markers other than dpn,wor https://journals.biologists.com/view-large/figure/7693046/4297tbl1a.jpeg
FeaturePlot(merged_obj, features = c("rna_GBIM-06054","rna_GBIM-06848","rna_esg","rna_numb","rna_pros","rna_polo","rna_shg")) # dpn:GBIM-06054, mira:GBIM-06848
#also checked cell cycle genes, to compare with what we see in D.melanogaster medulla Nb cluster for the same genes
FeaturePlot(merged_obj, features = c("rna_CycA","rna_stg","rna_CycE","rna_E2f1","rna_CycD","rna_Cdk2","rna_Cdk4"))
#chose cluster 56 as medulla NBs

#### choοse and save medulla NBs####
# find lamina and lobula plate gcm, eya, sim, tll, dac, acj6
FeaturePlot(merged_obj, features = c("rna_eya","rna_acj6","rna_gcm","rna_sim","rna_dac","rna_GBIM-04903"), raster = F) # GBIM-04903: tll

nw_object <- c(56)
NbS <- subset(merged_obj, idents = nw_object)
rm(larvaOL.integrated)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)

NbS <- RunPCA(NbS, npcs = 150)
ElbowPlot(NbS, ndims = 150)
NbS <- RunUMAP(NbS, dims = 1:150)
NbS <- FindNeighbors(NbS, dims = 1:150)
NbS <- FindClusters(NbS, resolution = 0.8)
DimPlot(NbS, label = TRUE) + NoLegend()

rename <- c("GBIM-01738" = "D",
            "GBIM-16025" = "B-H1",
            "GBIM-16201" = "Oaz",
            "GBIM-12506" = "hth"
)

# rename rownames in the matrix
rename_features_in_matrix <- function(mat, rename_map) {
  rn <- rownames(mat)
  rn[rn %in% names(rename_map)] <- rename_map[rn[rn %in% names(rename_map)]]
  rownames(mat) <- rn
  return(mat)
}

# renaming to all assay layers
NbS@assays$RNA@counts     <- rename_features_in_matrix(NbS@assays$RNA@counts, rename)
NbS@assays$RNA@data       <- rename_features_in_matrix(NbS@assays$RNA@data, rename)
NbS@assays$RNA@scale.data <- rename_features_in_matrix(NbS@assays$RNA@scale.data, rename)

saveRDS(NbS, "/gryllus_medulla_NBs.rds")

