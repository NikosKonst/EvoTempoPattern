# This R file contains all the steps for the analysis of single-cel mRNA seq data, from preprocessing to quality control and filtering,
#that leads to the generation of the species whole optic lobe data in the tribolium_OL.rds 
#It also contains the choice of the medulla neuroblasts cluster based on markers tribolium_medulla_OL.rds.
#Some genes that do not have an 1:1 ortholog have been identified via Flybase orthologs, or orthoDB and blast and have been reanmed and stored inside the objects.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

larval1.data <- Read10X(data.dir = "/Tcastaneum_larva1_all_GeneExt_merged_subsample100m_GCA_annot/outs/filtered_feature_bc_matrix/")

larval2.data <- Read10X(data.dir = "/Tcastaneum_larva2_all_GeneExt_merged_subsample100m_GCA_annot/outs/filtered_feature_bc_matrix/")

larval3.data <- Read10X(data.dir = "/Tcastaneum_larva3_all_GeneExt_merged_subsample100m_GCA_annot/outs/filtered_feature_bc_matrix/")

#rename one to one ortholog genes in each library, based on TABLE S1
mgi <- read.csv(data.dir = "/LOC_combined_OtO_TRICA_correspondence_GE_subsample.tsv", header = F)

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
larv1 = NormalizeData(larv1, normalization.method = "LogNormalize", scale.factor = 10000)
larv1 = FindVariableFeatures(larv1, selection.method = "vst", nfeatures = 2000)

larv2 <- subset(larv2, subset = nFeature_RNA > 1700)
larv2 = NormalizeData(larv2, normalization.method = "LogNormalize", scale.factor = 10000)
larv2 = FindVariableFeatures(larv2, selection.method = "vst", nfeatures = 2000)

larv3 <- subset(larv3, subset = nFeature_RNA > 1700)
larv3 = NormalizeData(larv3, normalization.method = "LogNormalize", scale.factor = 10000)
larv3 = FindVariableFeatures(larv3, selection.method = "vst", nfeatures = 2000)


reference.list <- c(larv1, larv2, larv3)
opticlobes.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:150)

larvaOL.integrated <- IntegrateData(anchorset = opticlobes.anchors, dims = 1:150)

rm(reference.list, opticlobes.anchors)

larvaOL.integrated <- ScaleData(object = larvaOL.integrated, verbose = FALSE)
larvaOL.integrated <- RunPCA(object = larvaOL.integrated, npcs = 150, verbose = FALSE)
larvaOL.integrated <- RunUMAP(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindNeighbors(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated<-FindClusters(object = larvaOL.integrated, resolution = 6)

saveRDS(larvaOL.integrated, "/tribolium_OL.rds")
larvaOL.integrated <- readRDS("/tribolium_OL.rds")

DimPlot(larvaOL.integrated, label = TRUE, raster=F, shuffle = T) + NoLegend()
DimPlot(larvaOL.integrated, split.by = "orig.ident",label = TRUE) + NoLegend()

#### choοse and save medulla NBs ####
# find medulla NB cluster based on conserved markers
FeaturePlot(larvaOL.integrated, features = "rna_shg") # epithel, NBs
FeaturePlot(larvaOL.integrated, features = "rna_repo") # glia
FeaturePlot(larvaOL.integrated, features = "rna_TcasGA2-TC008437", raster = F) # TcasGA2-TC008437: ase gmcs
FeaturePlot(larvaOL.integrated, features = "rna_TcasGA2-TC007824", raster = F) # TcasGA2-TC007824: mira, NBs
FeaturePlot(larvaOL.integrated, features = "rna_dpn") #NBs
FeaturePlot(larvaOL.integrated, features = "rna_TcasGA2-TC014474") # TcasGA2-TC014474:wor NBs, this gene hits both in uniprot protein blast to wor+esg, both are Snail TF family members 

# find lamina and lobula plate gcm, eya, sim, tll, dac, acj6
FeaturePlot(larvaOL.integrated, features = c("rna_TcasGA2-TC014730","rna_eya","rna_TcasGA2-TC016205","rna_tll","rna_dac","rna_TcasGA2-TC003196")) #gcm:TcasGA2_TC014730 , sim:TcasGA2_TC016205(found in uniprot, and ibeetle and ncbi protein), acj6:TcasGA2_TC003196

nw_object <- c(25,61)
NbS <- subset(larvaOL.integrated, idents = nw_object)
NbS <- subset(NbS, idents = nw_object)
NbS <- FindVariableFeatures(NbS, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = TRUE)
ScaleData(object = NbS, verbose = FALSE)
NbS <- RunPCA(NbS, npcs = 50)
ElbowPlot(NbS, ndims = 50)
NbS <- RunUMAP(NbS, dims = 1:50)
NbS <- FindNeighbors(NbS, dims = 1:50)
NbS <- FindClusters(NbS, resolution = 2)

#2nd round cleaning
nw_object <- c(2,4,6,7) 
NbS <- subset(NbS, idents = nw_object)
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
rename <- c("TcasGA2-TC031415" = "B-H1",
            "TcasGA2-TC032830" = "ey"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(NbS[["RNA"]]) <- rnames
saveRDS(NbS, "/tribolium_medulla_NBs.rds")

# remove low quality clusters (that do not integrate properly among the three libraries) from musca_OL.rds for visualization
# remove clusters 2,5,6,21,23,64 because they have low counts based on VlnPlot and check DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE, split.by = "orig.ident", raster=F) + NoLegend()
# clusters that have low counts but they form clear groups of cells in all three objects where kept
rm_clust <- c(2,5,6,21,23,64)
int_trib_filt <- subset(larvaOL.integrated, idents = setdiff(levels(Idents(larvaOL.integrated)),as.character(rm_clust)))
DimPlot(int_trib_filt, reduction = "umap", label = TRUE, raster=F, shuffle = T) + NoLegend()
saveRDS(int_trib_filt,"/trib_OL.rds")

#rename known terminal selectors gene names from D. melanogaster inside the nasonia_OL objects, even if not 1:1 orthologs, based on FLybase orthologs, orthoDB and blastp, to plot them easier
tribolium <-readRDS("/tribolium_OL.rds")
#### rename selectors nasonia ####
DefaultAssay(tribolium) <- "RNA"
# Extract RNA assay
rna_assay <- tribolium[["RNA"]]
# Get current feature (gene) names
rnames <- rownames(rna_assay)
# Your rename map
rename <- c("TcasGA2-TC003196" = "acj6",
            "TcasGA2-TC031040" = "ara",
            "TcasGA2-TC003238" = "Awh",
            "TcasGA2-TC003627" = "bab1",
            "TcasGA2-TC031038" = "bab2",
            "TcasGA2-TC007335" = "CG32532",
            "TcasGA2-TC034566" = "CG4328",
            "myc" = "dm",
            "TcasGA2-TC003191" = "Ets65A",
            "TcasGA2-TC032830" = "ey",
            "TcasGA2-TC032157" = "foxo",
            "TcasGA2-TC002783" = "klu",
            "TcasGA2-TC030900" = "Lin29",
            "TcasGA2-TC032451" = "mirr",
            "TcasGA2-TC006824" = "pdm3",
            "TcasGA2-TC031853" = "rn",
            "TcasGA2-TC013501" = "salm",
            "TcasGA2-TC016205" = "sim",
            "TcasGA2-TC033194" = "tj",
            "TcasGA2-TC007409" = "toy",
            "TcasGA2-TC012322" = "tsh",
            "TcasGA2-TC014798" = "vfl",
            "TcasGA2-TC011311" = "vvl"
)
# Replace where there's a match
rnames[rnames %in% names(rename)] <- rename[rnames[rnames %in% names(rename)]]
# Assign back *only* to the RNA assay
rownames(tribolium[["RNA"]]) <- rnames

saveRDS(tribolium,"/tribolium_OL.rds")
