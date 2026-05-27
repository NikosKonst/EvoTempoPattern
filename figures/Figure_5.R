# This R file contains all the information for generation of Figure 3
# NbS inputs are provided in the GEO, accession number GSE333138, except melanogaster_medulla_NBs.rds that is provided in the current githubpage.
# requires the custom order function in heatmap.R to run before call it to generate heatmaps

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#### Figure 5B,5C feature plots ####
##Virilis
Virilis <- readRDS("/virilis_medulla_NBs.rds")
FeaturePlot(Virilis,
            features = "rna_tll",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)
FeaturePlot(Virilis,
            features = "rna_CG13894",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_vir_tll.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig5C_vir_CG13894.png", units = "cm", width = 15, height = 15, bg = "white")

## Nasonia
Nasonia_nb <-readRDS("/nasonia_medulla_NBs.rds")
FeaturePlot(Nasonia,
            features = "rna_Sox21b",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5B_nas_sox21b_6.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(Nasonia,
            features = "rna_LOC100678366",
            cols = c("lightgrey", "#228B22"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5B_nas_sox_6.png", units = "cm", width = 15, height = 15, bg = "white")##9C27B0

FeaturePlot(Nasonia_nb,
            features = "rna_LOC100122379", #ken:LOC100122379
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6)
ggsave("/plots/Fig5C_ken.png", units = "cm", width = 15, height = 15, bg = "white")

## Musca
musca <-readRDS("/musca_medulla_NBs.rds") #131800466 ey, LOC101888960:Oaz fails 0 but looks temp
FeaturePlot(musca_nb,
            features = "rna_tin",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6)
ggsave("/plots/Fig5C_tin.png", units = "cm", width = 15, height = 15, bg = "white")

## Cricket
Cricket <- readRDS("gryllus_medulla_NBs.rds")
FeaturePlot(Cricket,
            features = "rna_Sp1",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_crick_sp1.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(Cricket,
            features = "rna_tup",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_crick_tup.png", units = "cm", width = 15, height = 15, bg = "white")

## Tribolium
Tribolium <- readRDS("/tribolium_medulla_NBs.rds")
FeaturePlot(Tribolium,
            features = "rna_Sox21b",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5B_trib_sox21b.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(Tribolium,
            features = "rna_D",
            cols = c("lightgrey", "#228B22"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5B_trib_D.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(Tribolium,
            features = "rna_Sox14",
            cols = c("lightgrey", "#228B22"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5B_trib_sox14_6.png", units = "cm", width = 15, height = 15, bg = "white")

FeaturePlot(Tribolium,
            features = "rna_Optix",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_trib_Optix.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(Tribolium,
            features = "rna_btd",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_trib_btd_6.png", units = "cm", width = 15, height = 15, bg = "white")

FeaturePlot(Tribolium,
            features = "rna_TcasGA2-TC011697",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",
            pt.size = 6,order=T)
ggsave("/plots/Fig5C_sp1_6.png", units = "cm", width = 15, height = 15, bg = "white")

#### Figure 5 heatmaps ####
library(smoother)
library(splus2R)
library(stringr)
library(zoo) # for moving average smoothing
library(viridis)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(dplyr)
library(scales)
library(pheatmap)

Nb_mel <- readRDS("/melanogaster_medulla_NBs.rds")
a <- c("slp1","slp2","B-H1","B-H2")
b <- c("Oaz","hbn","B-H1")
Nb<- Nb_mel

Virilis_NB <- readRDS("/virilis_medulla_NBs.rds")
a <- c("slp1","slp2","B-H1","B-H2")
b <- c("B-H1","CG13894","tll")
Nb<- Virilis_NB

musca_nb <-readRDS("/musca_medulla_NBs.rds")
a <- c("slp1","LOC101888975","B-H1","LOC101892452") #
b <- c("tin", "B-H1", "tll")
Nb<- musca_nb

Tribolium_nb <- readRDS("/tribolium_medulla_NBs.rds")
a <- c("Sox14","D","Sox21b","B-H1")
b <- c("B-H1","Optix","tll")
c <- c("dmrt99B","hth","btd","TcasGA2-TC011697") #TcasGA2-TC011697:sp1
Nb<- Tribolium_nb

Nasonia_nb <-readRDS("/nasonia_medulla_NBs.rds")
a <- c("LOC100678366","Sox21b","B-H1","LOC100122379","tll") #LOC100678366:sox, LOC100122379:ken
Nb<- Nasonia_nb

Cricket_nb <- readRDS("/gryllus_medulla_NBs.rds")
a <- c("hth","Sp1","B-H1","tup")
Nb<- Cricket_nb

Cloeon_nb <- readRDS("/cloeon_medulla_NBs.rds")
a <- c("Oaz","hbn","B-H1") #oaz:L

Bombyx_NB <- readRDS("/konstantina/Figures/Figure2/bombyx_renamed_tTF_NB.rds")
a <- c("Oaz","hbn","BH1") #oaz:L
Nb<- Bombyx_NB

##
##run if object cannot be converted in sce class and then NB_sce <- as.SingleCellExperiment(Nb, assay = "RNA")
Nb[["RNA"]] <- as(Nb[["RNA"]], Class="Assay")
# otherwise run directly
NB_sce <- as.SingleCellExperiment(Nb, assay = "RNA")
# run slingshot to estimate pseudotime
rd = Nb@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1
fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

#custom order function in heatmap.R, example for musca, the same was run for all objects
ordered_heat_max(
  object = NB_sce,
  features = a,
  window_size = 0.05,             # optional: window size for smoothing
  alpha_val = 0.5,                # optional: alpha value for smoothing
  span.value = 21,
  byMax = FALSE,
  outfolder = "/konstantina/Figures/Figure5",
  filename = "musca_BH2_heatmap"
)




