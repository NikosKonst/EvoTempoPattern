# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering,
#that leads to the generation of the species whole optic lobe data in the nasonia_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers nasonia_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

larval1.data <- Read10X(data.dir = "/Nvitripennis_larva1_all_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

larval2.data <- Read10X(data.dir = "/Nvitripennis_larva2_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

larval3.data <- Read10X(data.dir = "/Nvitripennis_larva3_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_NASVI_correspondence_GE_subsample.tsv", header = F)

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
  # Histogram for nFeature_RNA with a red dashed line at x = 300 for library1, and 380 for library 2 and 3
  ggplot2::qplot(
    obj$nFeature_RNA,
    geom = "histogram",
    bins = 101,
    fill = I("white"),
    col  = I("black")
  ) +
    ggplot2::geom_vline(
      xintercept = 300,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 300, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 300)
larv1 = NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 = FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 380)
larv2 <-  NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 380)
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

saveRDS(larvaOL.integrated, "/nasonia_OL.rds")
larvaOL.integrated <- readRDS("/nasonia_OL.rds")

DimPlot(larvaOL.integrated, label = TRUE, raster=F, shuffle = T) + NoLegend()
DimPlot(larvaOL.integrated, split.by = "orig.ident",label = TRUE) + NoLegend()

#### chose and save medulla NBs ####
# find medulla NB cluster based on conserved markers
FeaturePlot(larvaOL.integrated, features = "rna_shg") # epithel, NBs
FeaturePlot(larvaOL.integrated, features = "rna_repo") # glia
FeaturePlot(larvaOL.integrated, features = "rna_mira", raster = F) # TcasGA2-TC007824: mira, NBs
FeaturePlot(larvaOL.integrated, features = "rna_dpn") #NBs
FeaturePlot(larvaOL.integrated, features = "rna_LOC100114879") # wor NBs,
FeaturePlot(larvaOL.integrated,features = c("rna_dpn","rna_mira","rna_LOC100114879"))

#cluster 23
nw_object = 23
NbS <- subset(larvaOL.integrated, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)
DimPlot(Nb, label = TRUE, raster=F) + NoLegend()

# remove cluster 1 because it has pics of tll,dac and hth likely coming from another nb pool
nw_object = c(0,2,3)
NbS <- subset(Nb, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)
DimPlot(NbS, label = TRUE, raster=F) + NoLegend()

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")

saveRDS(NbS, "/nasonia_medulla_NBs.rds")

# remove low quality clusters (that do not integrate properly among the three libraries) from musca_OL.rds for visualization
# remove clusters 0,1,2,3,4,5,7,11 because they have low counts based on VlnPlot and check DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE, split.by = "orig.ident", raster=F) + NoLegend()
# clusters that have low counts but they form clear groups of cells in all three objects where kept
rm_clust <- c(0,1,2,3,4,5,7,11)

int_nas_filt <- subset(larvaOL.integrated, idents = setdiff(levels(Idents(int_nas)),as.character(rm_clust)))
DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE, raster=F, shuffle = T) + NoLegend()
FeaturePlot(larvaOL.integrated, features=c("rna_repo","rna_elav"))
saveRDS(larvaOL.integrated,"/nasonia_OL.rds")

#rename known terminal selectors gene names from D. melanogaster inside the nasonia_OL objects, even if not 1:1 orthologs, based on FLybase orthologs, orthoDB and blastp, to plot them easier
nasonia <-readRDS("/nasonia_OL.rds")
#### rename selectors nasonia ####
DefaultAssay(nasonia) <- "RNA"
# Extract RNA assay
rna_assay <- nasonia[["RNA"]]
# Get current feature (gene) names
rnames <- rownames(rna_assay)
# Your rename map
rename <- c("LOC101740863" = "ab",
            "LOC100122175" = "apt",
            "LOC107980969" = "Awh",
            "LOC100124164" = "bab1",
            "LOC100121348" = "caup",
            "LOC100121564" = "CG32105",
            "LOC100118780" = "CG34340",
            "LOC100121564" = "CG4328",
            "LOC116415728" = "CG9650",
            "LOC100123387" = "ct",
            "LOC100678366" = "D",
            "LOC116417139" = "danr",
            "LOC103316520" = "disco-r",
            "LOC100679741" = "dm",
            "LOC100302336" = "dmrt99B",
            "LOC100119062" = "erm",
            "LOC100123890" = "foxo",
            "LOC100122721" = "ham",
            "LOC100117267" = "klu",
            "LOC100116733" = "Lim1",
            "LOC100120126" = "Lim3",
            "LOC100679002" = "Lin29",
            "LOC100118014" = "lov",
            "LOC100115265" = "luna",
            "LOC100678062" = "rn",
            "LOC100116939" = "salm",
            "LOC100117203" = "slou",
            "LOC100117833" = "tj",
            "LOC100116264" = "tsh",
            "LOC100118683" = "vfl",
            "LOC100114448" = "Vsx1", 
            "LOC100120173" = "zfh1",
            "LOC100116002" = "zfh2"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(nasonia[["RNA"]]) <- rnames
saveRDS(nasonia,"/nasonia_OL.rds")

