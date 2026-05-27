# This R file contains all the information for generation of Figure6 and Supplementary Figure S13
# whole optic lobe data as inputs are provided in the GEO, accession number GSE333138, except for melanogaster data OL.combined_published.annotations.rds from Konstantinides et al. 2022,Nature

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#### mel, load the file  OL.combined_published.annotations.rds from Konstantinides et al. 2022,Nature doi: 10.1038/s41586-022-04564 ####
mel <- readRDS("/data/old_data/datasets/Neuronal_Diversity_suppl_scripts_data/OL.combined_published.annotations.rds")
mel <- UpdateSeuratObject(mel)
musca <- readRDS("/musca_OL.rds")

#### Figure 6A ####
FeaturePlot(musca,
            features = "rna_ab",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",cells = "A")
ggsave("/plots/Fig6A_ab_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_ab",
            cols = c("lightgrey", "#8E24AA"),reduction="umap")
ggsave("/plots/Fig6A_ab_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_tin",
            cols = c("lightgrey", "#8E24AA"),reduction="umap")
ggsave("/plots/Fig6A_tin_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_tin",
            cols = c("lightgrey", "#8E24AA"),reduction="umap")
ggsave("/plots/Fig6A/tin_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_Poxm",
            cols = c("lightgrey", "#8E24AA"),reduction="umap")
ggsave("/plots/Fig6A_Poxm_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_Poxm",
            cols = c("lightgrey", "#8E24AA"),reduction="umap",pt.size=2) #also checked poxn and is not expressed in the melanogaster optic lobe data

#### Figure 6B ####
FeaturePlot(musca,
            features = "rna_Lim3",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6B_Lim3_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_Lim3",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6B_Lim3_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_tup",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6B_tup_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_tup",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6B_tup_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

#### Figure 6C ####
FeaturePlot(musca,
            features = "rna_Ets65A",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6C_ets65_musl_with_ey.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_Ets65A",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6C_Ets65A_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_fd59A",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6C_fd59A_like_mus_with_ey.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = c("rna_fd59A"),
            reduction="umap")
ggsave("/plots/Fig6C_fd59A_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_Dll",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6C_Dll_like_mus_with_ey.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)
FeaturePlot(mel,
            features = "rna_Dll",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/Fig6C_Dll_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

#### Figure S13 ####
# musca plots #
FeaturePlot(musca,
            features = "rna_ey",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_ey.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_vvl",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_vvl_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_kn",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_kn_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_erm",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_erm_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_TfAP-2",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_TfAP-2_mus_with_ey.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(musca,
            features = "rna_svp",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_svp_musca.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

# mel plots #
FeaturePlot(mel,
            features = "rna_ey",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_ey_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(mel,
            features = "rna_vvl",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_vvl_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(mel,
            features = "rna_kn",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_kn_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(mel,
            features = "rna_erm",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_erm_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(mel,
            features = "rna_TfAP-2",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_TfAP-2_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)

FeaturePlot(mel,
            features = "rna_svp",
            cols = c("lightgrey", "#000080"),reduction="umap")
ggsave("/plots/FigS13_svp_mel.png", units = "cm", width = 25, height = 20, bg = "white",dpi=300)


