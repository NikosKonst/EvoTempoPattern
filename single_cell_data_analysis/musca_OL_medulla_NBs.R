# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering, also DoubletFinder was used to identify and remove doublets.
#that leads to the generation of the species whole optic lobe data in the musca_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers musca_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

larval1.data <- Read10X(data.dir = "/Mdomestica_larva1_all_GeneExt_merged_subsample100m_Ey_orphan/outs/filtered_feature_bc_matrix/")

larval2.data <- Read10X(data.dir = "/Mdomestica_larva2_all_GeneExt_merged_subsample100m_Ey_orphan/outs/filtered_feature_bc_matrix/")

larval3.data <- Read10X(data.dir = "/Mdomestica_larva3_all_GeneExt_merged_subsample100m_Ey_orphan/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_MUSDO_correspondence_GE_subsample_Ey.csv", header = F)

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
      xintercept = 1300,
      color      = "red",
      linetype   = "dashed",
      linewidth  = 1
    )
}

Vln_thres <- function(obj) {
  VlnPlot(obj, "nFeature_RNA", pt.size = 0) +
    geom_hline(yintercept = 1300, # or different
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
larv1 <- subset(larv1, subset = nFeature_RNA > 1300)
larv1 <- NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 <- FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 1300)
larv2 <- NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 <- FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 1300)
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
nExp_poi <- round(0.228*nrow(lib1@meta.data))  ## the expected number of doublet for ~27589 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib1 <- doubletFinder(lib1,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib1@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plota = DimPlot(lib1, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_5710")
plotb = FeaturePlot(lib1, features = "rna_dpn")
plota + plotb

# number of singlets and doublets
table(lib1@meta.data$DF.classifications_0.25_0.3_5710)
#Doublet Singlet 
#5710   21874 

object_NoDoublets <- subset(lib1, subset =  DF.classifications_0.25_0.3_5710 == "Singlet")

saveRDS(object_NoDoublets , "/lib1_singlets_ey_1300_filt.rds")
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
nExp_poi <- round(0.32*nrow(lib2@meta.data))  ## the expected number of doublet for ~31534 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib2<- doubletFinder(lib2,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib2@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plotc = DimPlot(lib2, reduction = 'umap', group.by = "DF.classifications_0.25_0.2_8775")
plotd = FeaturePlot(lib2, features = "rna_dpn")
plotc + plotd

# number of singlets and doublets
table(lib2@meta.data$DF.classifications_0.25_0.2_8775)
#Doublet Singlet 
#8775   22753 

object_NoDoublets <- subset(lib2, subset =  DF.classifications_0.25_0.2_8775 == "Singlet")

saveRDS(object_NoDoublets , "/lib2_singlets_ey_1300_filt.rds")
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
nExp_poi <- round(0.308*nrow(lib3@meta.data))  ## the expected number of doublet for ~30229 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # to adjust by homotypic doublet 
lib3<- doubletFinder(lib3,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)  
lib3@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
plote = DimPlot(lib3, reduction = 'umap', group.by = "DF.classifications_0.25_0.16_8312")
plotf = FeaturePlot(lib3, features = "rna_dpn")
plote + plotf

# number of singlets and doublets
table(lib3@meta.data$DF.classifications_0.25_0.16_8312)
#Doublet Singlet 
#8312   21975 

object_NoDoublets <- subset(lib3, subset =  DF.classifications_0.25_0.16_8312 == "Singlet")

saveRDS(object_NoDoublets , "/lib3_singlets_ey_1300_filt.rds")
rm(lib3, object_NoDoublets)

#### integration singlets from each library ####
larv1 <- readRDS("/lib1_singlets_ey_1300_filt.rds")

singlets1 <- larv1[["RNA"]]$counts

signlets1 <- CreateSeuratObject(counts =singlets1, project="lib1")
rm(larv1)

signlets1 <- NormalizeData(signlets1)
signlets1 <- FindVariableFeatures(signlets1)
VlnPlot(signlets1, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0)

larv2 <- readRDS("/lib2_singlets_ey_1300_filt.rds")

singlets <- larv2[["RNA"]]$counts

signlets2 <- CreateSeuratObject(counts =singlets, project="lib2")
rm(larv2)

signlets2 <- NormalizeData(signlets2)
signlets2 <- FindVariableFeatures(signlets2)

larv3 <- readRDS("/lib3_singlets_ey_1300_filt.rds")
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
larvaOL.integrated<-FindClusters(object = larvaOL.integrated, resolution = 4)

saveRDS(larvaOL.integrated, "/musca_OL.rds")
rm(signlets1,signlets2,signlets3,larvaOL.integrated)

larvaOL.integrated <- readRDS("/musca_OL.rds")

DimPlot(larvaOL.integrated, label = TRUE, raster=F, shuffle = T) + NoLegend()
DimPlot(larvaOL.integrated, split.by = "orig.ident",label = TRUE) + NoLegend()

#### choοse and save medulla NBs ####
# find medulla NB cluster based on conserved markers
FeaturePlot(larvaOL.integrated, features = "rna_mira", raster = F) # mira, NBs
FeaturePlot(larvaOL.integrated, features = "rna_dpn") #NBs
#could not idenitfy wor homolog

# find lamina and lobula plate gcm, eya, sim, tll, dac, acj6
FeaturePlot(larvaOL.integrated, features = c("rna_gcm","rna_eya","rna_sim","rna_tll","rna_dac","rna_acj6")) 

nw_object <- c(20)
NbS <- subset(larvaOL.integrated, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)
DimPlot(NbS, label = TRUE) + NoLegend()

#2nd round cleaning
nw_object <- c(5,3,2,1) # 0,4,6 more likely gmcs based on higher levels of expression of elav (neuronal marker), nerfin-1(higher at gmcs), and repetition of windows like hth, D,slp1
NbS <- subset(NbS, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 0.8)

DefaultAssay(NbS) <- "RNA"

# Extract RNA assay
rna_assay <- NbS[["RNA"]]

# Get current feature (gene) names
rnames <- rownames(rna_assay)

# Your rename map
rename <- c("LOC101894095" = "opa",
            "LOC131804155" = "scro",
            "LOC101888975" = "slp2"
)

# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]

# Assign back *only* to the RNA assay
rownames(NbS[["RNA"]]) <- rnames
saveRDS(NbS, "/musca_medulla_NBs.rds")

# remove low quality clusters (that do not integrate properly among the three libraries) from musca_OL.rds for visualization
# remove clusters 0,1,2,3,5,7,10,12,14,19,31,33,35,36,76 because they have low counts based on VlnPlot and check DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE, split.by = "orig.ident", raster=F) + NoLegend()
# clusters that have low counts but they form clear groups of cells in all three objects where kept
rm_clust <- c(0,1,2,3,5,7,10,12,14,19,31,33,35,36,76)
int_mus_filt <- subset(larvaOL.integrated, idents = setdiff(levels(Idents(larvaOL.integrated)),as.character(rm_clust)))
DimPlot(int_mus_filt, reduction = "umap", label = TRUE, raster=F, shuffle = T) + NoLegend()
saveRDS(int_mus_filt,"/musca_OL.rds")

#rename known terminal selectors gene names from D. melanogaster inside the musca_OL objects, even if not 1:1 orthologs, based on FLybase orthologs, orthoDB and blastp, to plot them easier
musca <-readRDS("/musca_OL.rds")
#### rename selectors musca ####
DefaultAssay(musca) <- "RNA"
# Extract RNA assay
rna_assay <- musca[["RNA"]]
# Get current feature (gene) names
rnames <- rownames(rna_assay)
# Your rename map
rename <- c("LOC101894095" = "opa",
            "LOC131804155" = "scro",
            "LOC101888975" = "slp2",
            "LOC101895495" = "ct",
            "LOC101892502" = "vvl",
            "LOC105262121" = "pros",
            "LOC101889122" = "Ets65A",
            "LOC101891159" = "bi",
            "LOC101889596" = "CG9650",
            "LOC101887740" = "ab",
            "LOC101891515" = "apt",
            "LOC101896023" = "bab1",
            "LOC101899928" = "bab2",
            "LOC101892866" = "br",
            "LOC101895484" = "Camta",
            "LOC101890442" = "CG32105",
            "LOC101890442" = "CG32105",
            "LOC109613748" = "CG34340",
            "LOC101891125" = "CG4328", 
            "LOC101901331" = "CG43689", 
            "LOC101894034" = "CG9932", 
            "LOC101894733" = "dac", 
            "LOC101899406" = "disco-r", 
            "LOC101887976" = "dm", 
            "LOC101899338" = "dve", 
            "LOC101889153" = "foxo", 
            "LOC101895855" = "grn", 
            "LOC101890702" = "ham",
            "LOC101899402" = "klu", 
            "LOC101895069" = "mirr",
            "LOC101887472" = "oc", 
            "LOC101900955" = "pdm3",
            "LOC105262410" = "rn", 
            "LOC101895997" = "salm",  
            "LOC101895647" = "salr", 
            "LOC101895624" = "slou", 
            "LOC101894180" = "Sox102F", #1:1 from gabaldon as well
            "LOC131800511" = "sr", 
            "LOC105262321" = "Stat92E",  
            "LOC101894165" = "vfl",
            "LOC101898105" = "zfh1",
            "LOC101893110" = "zfh2"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(musca[["RNA"]]) <- rnames
saveRDS(musca,"/musca_OL.rds")


