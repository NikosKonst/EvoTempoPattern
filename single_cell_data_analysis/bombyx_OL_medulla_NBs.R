# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering, also DoubletFinder was used to identify and remove doublets.
#that leads to the generation of the species whole optic lobe data in the bombyx_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers bombyx_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

larval1.data <- Read10X(data.dir = "/Bmori_larva1_all_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

larval2.data <- Read10X(data.dir = "/Bmori_larva2_all_GeneExt_merged_subsample100/outs/filtered_feature_bc_matrix/")

larval3.data <- Read10X(data.dir = "/Bmori_larva3_all_GeneExt_merged_subsample100m/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/all_genes_OTO_ortho_bombyx.csv", header = F)

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
      xintercept = 1700,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 1700, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 1700)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 1700)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 1700)
larv3 <- NormalizeData(larv3, normalization.method = "LogNormalize", scale.factor = 10000)
larv3 <- FindVariableFeatures(larv3, selection.method = "vst", nfeatures = 2000)

# identify and remove doublets for each library
library(DoubletFinder)
lib1 <- larv1
lib1 <- ScaleData(object = lib1)
lib1 <- RunPCA(object = lib1 , npcs = 50)
lib1 <- FindNeighbors(object = lib1, dims = 1:50)
lib1 <- FindClusters(object = lib1)
lib1 <- RunUMAP(object = lib1, dims = 1:50)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib1 <- paramSweep(lib1, PCs = 1:50, sct = FALSE) 
sweep.res.list_lib1 
sweep.stats_lib1 <- summarizeSweep(sweep.res.list_lib1, GT = FALSE)
sweep.stats_lib1   
bcmvn_lib1 <- find.pK(sweep.stats_lib1)  

ggplot(bcmvn_lib1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))  # choose the first value from the list of pk 
pK  #

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib1@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           ## ex: model the % of homotypic doublet based on the user provided annotations 
homotypic.prop
nExp_poi <- round(0.176*nrow(lib1@meta.data))  ## the expected number of doublet for ~22000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib1 <- doubletFinder(lib1,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib1@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plota = DimPlot(lib1, reduction = 'umap', group.by = "DF.classifications_0.25_0.07_2679")
plotb = FeaturePlot(lib1, features = "rna_dpn")
plota + plotb

# number of singlets and doublets
table(lib1@meta.data$DF.classifications_0.25_0.07_2679)
#Doublet Singlet 
#2679    18775

object_NoDoublets <- subset(lib1, subset = DF.classifications_0.25_0.07_2679 == "Singlet")

saveRDS(object_NoDoublets , "/lib1_singlets.rds")
rm(lib1, object_NoDoublets)

##
lib2 <- larv2

lib2 <- ScaleData(object = lib2)
lib2 <- RunPCA(object = lib2 , npcs = 50)
lib2 <- FindNeighbors(object = lib2, dims = 1:50)
lib2 <- FindClusters(object = lib2)
lib2 <- RunUMAP(object = lib2, dims = 1:50)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib2 <- paramSweep(lib2, PCs = 1:50, sct = FALSE) 
sweep.res.list_lib2
sweep.stats_lib2 <- summarizeSweep(sweep.res.list_lib2, GT = FALSE)
sweep.stats_lib2   
bcmvn_lib2 <- find.pK(sweep.stats_lib2)  

ggplot(bcmvn_lib2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))  # choose the first value from the list of pk 
pK  #

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib2@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           ## ex: model the % of homotypic doublet based on the user provided annotations 
homotypic.prop
nExp_poi <- round(0.178*nrow(lib2@meta.data))  ## the expected number of doublet 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib2 <- doubletFinder(lib2,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib2@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plota = DimPlot(lib2, reduction = 'umap', group.by = "DF.classifications_0.25_0.08_3176")
plotb = FeaturePlot(lib2, features = "rna_dpn")
plota + plotb

# number of singlets and doublets
table(lib2@meta.data$DF.classifications_0.25_0.08_3176)
#Doublet Singlet 
#3176    22002

object_NoDoublets <- subset(lib2, subset = DF.classifications_0.25_0.08_3176 == "Singlet")
saveRDS(object_NoDoublets,  "/lib2_singlets.rds")

#rm(lib2, object_NoDoublets)

lib3 <- larv3

lib3 <- ScaleData(object = lib3)
lib3 <- RunPCA(object = lib3 , npcs = 50)
lib3 <- FindNeighbors(object = lib3, dims = 1:50)
lib3 <- FindClusters(object = lib3)
lib3 <- RunUMAP(object = lib3, dims = 1:50)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib3 <- paramSweep(lib3, PCs = 1:50, sct = FALSE) 

sweep.stats_lib3 <- summarizeSweep(sweep.res.list_lib3, GT = FALSE)

bcmvn_lib3 <- find.pK(sweep.stats_lib3)  

ggplot(bcmvn_lib3, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib3 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))  # choose the first value from the list of pk 
pK  #

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib3@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           ## ex: model the % of homotypic doublet based on the user provided annotations 
homotypic.prop
nExp_poi <- round(0.178*nrow(lib3@meta.data))  ## 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib3 <- doubletFinder(lib3,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib3@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plota = DimPlot(lib3, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_3540")
plotb = FeaturePlot(lib3, features = "rna_dpn")
plota + plotb

# number of singlets and doublets
table(lib3@meta.data$DF.classifications_0.25_0.3_3540)
#Doublet Singlet 
#3540    21912

object_NoDoublets <- subset(lib3, subset = DF.classifications_0.25_0.3_3540 == "Singlet")

saveRDS(object_NoDoublets , "/lib3_singlets.rds")
rm(lib3, object_NoDoublets)

#### integration singlets from each library ####
larv1 <- readRDS("/lib1_singlets.rds")

singlets1 <- larv1[["RNA"]]$counts

signlets1 <- CreateSeuratObject(counts =singlets1, project="lib1")
rm(larv1)

signlets1 <- NormalizeData(signlets1)
signlets1 <- FindVariableFeatures(signlets1)
VlnPlot(signlets1, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0)

larv2 <- readRDS("/lib2_singlets.rds")

singlets <- larv2[["RNA"]]$counts

signlets2 <- CreateSeuratObject(counts =singlets, project="lib2")
rm(larv2)

signlets2 <- NormalizeData(signlets2)
signlets2 <- FindVariableFeatures(signlets2)

larv3 <- readRDS("/lib3_singlets.rds")
singlets <- larv3[["RNA"]]$counts

signlets3 <- CreateSeuratObject(counts =singlets, project="lib3")
rm(larv3, singlets)

signlets3 <- NormalizeData(signlets3)
signlets3 <- FindVariableFeatures(signlets3)

reference.list <- c(signlets1, signlets2, signlets3)
opticlobes.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:150)
larvaOL.integrated <- IntegrateData(anchorset = opticlobes.anchors, dims = 1:150)
rm(reference.list, opticlobes.anchors)

larvaOL.integrated <- ScaleData(object = larvaOL.integrated, verbose = FALSE)
larvaOL.integrated <- RunPCA(object = larvaOL.integrated, npcs = 150, verbose = FALSE)

ElbowPlot(larvaOL.integrated, ndims = 150)

larvaOL.integrated <- RunUMAP(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindNeighbors(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindClusters(object = larvaOL.integrated, resolution = 8)

saveRDS(larvaOL.integrated, "/bombyx_OL.rds")
rm(signlets1,signlets2,signlets3,larvaOL.integrated)

larvaOL.integrated <- readRDS("/bomyx_OL.rds")

DimPlot(larvaOL.integrated, label = TRUE, raster=F, shuffle = T) + NoLegend()
DimPlot(larvaOL.integrated, split.by = "orig.ident",label = TRUE) + NoLegend()

#### choοse and save medulla NBs ####
# find medulla NB cluster based on conserved markers
FeaturePlot(larvaOL.integrated, features = "rna_mira", raster = F) # mira, NBs
FeaturePlot(larvaOL.integrated, features = "rna_dpn") #NBs
FeaturePlot(larvaOL.integrated, features = "rna_LOC101744234") #wor, NBs
FeaturePlot(larvaOL.integrated, features = "rna_shg") #epith, NBs
FeaturePlot(NbS, features = "rna_elav", raster = F) #neurons, already higher in gmcs
FeaturePlot(NbS, features = "rna_LOC101744918", raster = F) #nerfin:LOC101744918 higher in gmcs

# find lamina and lobula plate gcm, eya, sim, tll, dac, acj6
FeaturePlot(larvaOL.integrated, features = c("rna_LOC101742696","rna_eya","rna_sim","rna_tll","rna_dac","rna_LOC100750224")) #gcm:LOC101742696, acj6:100750224
FeaturePlot(larvaOL.integrated, features = c("rna_repo"), min.cutoff = 1.5) #glia

nw_object <- c(48)
NbS <- subset(larvaOL.integrated, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)
DimPlot(NbS, label = TRUE) + NoLegend()

DefaultAssay(NbS) <- "RNA"
rna_assay <- NbS[["RNA"]]
# Get current feature (gene) names
rnames <- rownames(rna_assay)
# Your rename map
rename <- c("LOC101747109" = "erm",
            "LOC101739648" = "Oaz"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(NbS[["RNA"]]) <- rnames

saveRDS(NbS, "/bombyx_medulla_NBs.rds")

# remove low quality clusters (that do not integrate properly among the three libraries) from musca_OL.rds for visualization
# remove clusters 0,2,3,4,5,7,9,15,24,28,29,36,39,73,74 because they have low counts based on VlnPlot and check DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE, split.by = "orig.ident", raster=F) + NoLegend()
# clusters that have low counts but they form clear groups of cells in all three objects where kept
rm_clust <- c(0,2,3,4,5,7,9,15,24,28,29,36,39,73,74)
int_bombyx_filt <- subset(larvaOL.integrated, idents = setdiff(levels(Idents(larvaOL.integrated)),as.character(rm_clust)))
DimPlot(int_bombyx_filt, reduction = "umap", label = TRUE, raster=F, shuffle = T) + NoLegend()
saveRDS(int_bombyx_filt,"/bombyx_OL.rds")



