# This R file contains all the information for generation of Figure 3
# NbS inputs are provided in the GEO, accession number GSE333138, except melanogaster_medulla_NBs.rds that is provided in the current githubpage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

musca <-readRDS("/musca_medulla_NBs.rds")
Nasonia <-readRDS("/nasonia_medulla_NBs.rds")
Aedes <- readRDS("/aedes_medulla_NBs.rds")
Tribolium <- readRDS("/triboliummedulla_NBs.rds.rds")
Virilis <- readRDS("/virilis_medulla_NBs.rds")
Cricket <- readRDS("/gryllus_medulla_NBs.rds")

####virilis####
FeaturePlot(Virilis,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Virilis,
            features = "rna_hbn",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)

FeaturePlot(Virilis,
            features = "rna_tll",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
ggsave("/plots/Fig3_vir_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_vir_hbn.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_vir_tll.png", units = "cm", width = 15, height = 15, bg = "white")

####musca####
musca <- musca_nb
FeaturePlot(musca,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(musca,
            features = "rna_tll",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)

FeaturePlot(musca,
            features = "rna_erm",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(musca,
            features = "rna_B-H1",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(musca,
            features = "rna_B-H1",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)
ggsave("/plots/Fig3_mus_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_mus_tll.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_mus_erm.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_mus_BH.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_mus_BH_green.png", units = "cm", width = 15, height = 15, bg = "white")

####aedes####
FeaturePlot(Aedes,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Aedes,
            features = "rna_slp2",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Aedes,
            features = "rna_hbn",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(Aedes,
            features = "rna_B-H1",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
ggsave("/plots/Fig3_ae_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_ae_slp.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_ae_hbn.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_ae_BH.png", units = "cm", width = 15, height = 15, bg = "white")

####tribolium####
FeaturePlot(Tribolium,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Tribolium,
            features = "rna_hbn",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Tribolium,
            features = "rna_B-H1",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)

FeaturePlot(Tribolium,
            features = "rna_B-H1",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Tribolium,
            features = "rna_opa",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(Tribolium,
            features = "rna_tll",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)

ggsave("/plots/Fig3_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_trib_hbn.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_trib_BH.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("plots/Fig3_trib_opa.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_trib_tll.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("plots/Fig3_trib_BH_green.png", units = "cm", width = 15, height = 15, bg = "white")

####nasonia####
FeaturePlot(Nasonia,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Nasonia,
            features = "rna_opa",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,order=T)

FeaturePlot(Nasonia,
            features = "rna_ey",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(Nasonia,
            features = "rna_scro",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
FeaturePlot(Nasonia,
            features = "rna_B-H1",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,order=T)
ggsave("/plots/Fig3_nas_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_nas_opa.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_nas_ey.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_nas_scro.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_nas_BH.png", units = "cm", width = 15, height = 15, bg = "white")

####cricket####
FeaturePlot(Cricket,
            features = "rna_hth",
            cols = c("lightgrey", "#228B22"),
            pt.size = 6,reduction="umap",order=T)
FeaturePlot(Cricket,
            features = "rna_D",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,reduction="umap",order=T)
FeaturePlot(Cricket,
            features = "rna_ey",
            cols = c("lightgrey", "#8E24AA"),
            pt.size = 6,reduction="umap",order=T)
ggsave("/plots/Fig3_cricket_hth.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_cricket_D.png", units = "cm", width = 15, height = 15, bg = "white")
ggsave("/plots/Fig3_cricket_ey.png", units = "cm", width = 15, height = 15, bg = "white")


