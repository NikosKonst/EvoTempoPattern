# This R file contains all the information for generation of Figure 2B and 2D
# NbS inputs are provided in the GEO, accession number GSE333138, except melanogaster_medulla_NBs.rds that is provided in the current githubpage.
# requires the custom order function in heatmap.R to run before call it in 2D to generate heatmaps

library(smoother)
library(splus2R)
library(stringr)
library(zoo) # for moving average smoothing
library(viridis)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(pheatmap)
library(dplyr)
library(scales)

#### Figure 2D ####
# medulla NBs of all species
Nb_mel <- readRDS("/melanogaster_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "opa", "erm","Oaz", "ey","hbn","D","B-H1")

Virilis_NB <- readRDS("/virilis_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "opa", "erm","Oaz", "ey","hbn","D","B-H1")

musca_nb <-readRDS("/musca_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "LOC101894095", "erm","LOC101888960", "ey","hbn","D","B-H1") #LOC101894095:opa, "LOC101888960:oaz

Aedes_nb <- readRDS("/aedes_medulla_NBs.rds") 
a <- c("hth", "dmrt99B", "opa", "erm","Oaz", "ey","hbn","D","B-H1") #renamed LOC with Dmel homolog:LOC5575091:hth,LOC5573207:Oaz, LOC5570655:B-H1 

Bombyx_NB <- readRDS("/bombyx_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "opa", "erm","Oaz","ey","hbn","D","LOC101747128") #B-H LOC101747128

Tribolium_nb <- readRDS("/tribolium_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "opa", "erm","TcasGA2-TC004682", "ey","hbn","D","B-H1")

Nasonia_nb <-readRDS("/nasonia_medulla_NBs.rds")
a <- c("hth", "LOC116416021","opa", "erm","PRAS40","ey","hbn","LOC100678366","B-H1") #LOC116416021 dmrt, LOC100678366:sox/D, Pras40: Lobe/Oaz

Cricket_nb <- readRDS("/gryllus_medulla_NBs.rds")
a <- c("hth","Oaz","GBIM-00865","ey","hbn","D","B-H1") #GBIM-00865 opa

Cloeon_nb <- readRDS("/cloeon_medulla_NBs.rds")
a <- c("hth", "dmrt99B", "opa", "erm","Oaz", "ey","hbn","D","B-H1") 


#### pseudotime ####
Nb<- Nb_mel

Nb<- Virilis_NB

Nb<- musca_nb

Nb<- Aedes_nb

Nb <- Bombyx_NB

Nb<- Tribolium_nb

Nb<- Nasonia_nb

Nb<- Cricket_nb

Nb<- Cloeon_nb

##run if object cannot be converted in sce class
Nb[["RNA"]] <- as(Nb[["RNA"]], Class="Assay")

# convert suerat objects in sce class
NB_sce <- as.SingleCellExperiment(Nb, assay = "RNA")

# run slingshot to estimate pseudotime
rd = Nb@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)

colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

#### generate heatmaps ####
#custom order function in heatmap.R, example for cloeon, the same was run for all objects
ordered_heat_max(
  object = NB_sce,
  features = a,             # you provide the list of genes defined below each object
  window_size = 0.05,             # optional: window size for smoothing
  alpha_val = 0.5,                # optional: alpha value for smoothing
  span.value = 21,
  byMax = FALSE,
  outfolder = "/plots/Figure2D",
  filename = "cloeon"
)

#### Figure 2B ####
#run the LayeredFeaturePlot function in LayeredFeaturePlot.R 
LayeredFeaturePlot(Virilis_NB, features = c("rna_hth","rna_ey","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6,
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_vir_layer_hth_ey_BH.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

LayeredFeaturePlot(musca_nb, features = c("rna_hth","rna_LOC101894095","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_musca_layer_hth_opa_BH_.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

# hth:LOC5575091, BH:LOC5570655
LayeredFeaturePlot(Aedes_nb, features = c("rna_LOC5575091","rna_opa","rna_LOC5570655"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_aedes_layer_hth_opa_bH.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)


LayeredFeaturePlot(Bombyx_NB, features = c("rna_hth","rna_opa","rna_scro"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_bombyx_layer_hth_opa_scro.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

LayeredFeaturePlot(Tribolium_nb, features = c("rna_hth","rna_opa","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_tribolium_layer_hth_opa_bH.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

LayeredFeaturePlot(Nasonia_nb, features = c("rna_hth","rna_ey","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_nasonia_layer_hth_ey_bh.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

LayeredFeaturePlot(Cricket_nb, features = c("rna_hth","rna_ey","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_cricket_layer_hth_ey_bh.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)

LayeredFeaturePlot(Cloeon_nb, features = c("rna_hth","rna_ey","rna_B-H1"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7"),pt.size = 6 , 
                   pt.brightness = 1, reduction = "umap")
ggsave("/plots/Fig2B_cloeon_layer_hth_ey_bh.png", units = "cm", width = 20, height = 25, bg = "white",dpi=300)


