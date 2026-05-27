# This R file contains all the information for generation of Supplementary Figure S9
# NbS inputs are provided in the GEO, accession number GSE333138
# tfs lists for each species are coming from Table_S2
# res_lists are provided in the current github page

suppressPackageStartupMessages({
  library(Seurat)
  library(PseudotimeDE)
  library(SingleCellExperiment)
  library(slingshot)
  library(tibble)
  library(dplyr)
  library(doParallel)
  library(ggplot2)
  library(scales)
  library(S4Vectors)
})

####virilis####
NbS <- readRDS("/virilis_medulla_NBs.rds")
tfs = read.csv("/Dvir_TFS.csv",  header = F)
res_list <- readRDS("/virilis_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS

keep= c("B-H1","B-H2","CG13894","D","dmrt99B","elB","ey","gcm","hbn","hth","L","LOC6623224","LOC6629252","LOC6629705","LOC6631725","LOC6636954","opa","slp1","slp2","tll")
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)+ scale_x_reverse()
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/virilis/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_B-H1","rna_B-H2","rna_CG13894","rna_D","rna_dmrt99B","rna_elB","rna_ey","rna_gcm","rna_hbn","rna_hth","rna_L","rna_LOC6623224","rna_LOC6629252","rna_LOC6629705","rna_LOC6631725","rna_LOC6636954","rna_opa","rna_slp1","rna_slp2","rna_tll")

FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/virilis/umap.png", units = "cm", width = 60, height = 60, bg = "white")

####musca####
NbS <- readRDS("/musca_medulla_NBs.rds")
tfs <- read.csv("/Mdom_TFS.csv",  header = F)
res_list <- readRDS("/musca_res_list.rds")


NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime
quantile_matrix = get_quantile_matrix(pseudotime, intervals = 20)
# its important to use sce not NB (seurat object), otherwise it is accesing the seurat integrated object
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("chinmo","ap","B-H1","Bx","cwo","D","dmrt99B","erm","hbn","hth","LOC101888960","LOC101888975","LOC101894095","LOC131800466","LOC131801858","LOC131804155","LOC131804616","LOC131806478","luna","sd","slp1","tin","tll","ey")
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/musca/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_ap","rna_chinmo","rna_B-H1","rna_Bx","rna_D","rna_dmrt99B","rna_erm","rna_hbn","rna_hth","rna_LOC101888960","rna_LOC101888975","rna_LOC101894095","rna_LOC131800466","rna_LOC131801858","rna_LOC131804155","rna_LOC131804616","rna_LOC131806478","rna_luna","rna_sd","rna_slp1","rna_tin","rna_tll","rna_ey")
FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/musca/umap.png", units = "cm", width = 15, height = 15, bg = "white")

####bombyx####
NbS <- readRDS("/bombyx_medulla_NBs.rds")
tfs = read.csv("/Bmori_TFS.csv",  header = F)
res_list <- readRDS("/bombyx_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime
quantile_matrix = get_quantile_matrix(pseudotime, intervals = 20)
# its important to use sce not NB (seurat object), otherwise it is accesing the seurat integrated object
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("Bmap-A","LOC100750224","LOC101739648","LOC101741268","LOC101741704","LOC101744414","LOC101747109","Optix","Sox102F","Sox21b","dmrt99B","elB","hbn","hth","opa","scro","so","toy","ey","slp2","LOC101747128","D","tll")
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/bombyx/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_Bmap-A","rna_LOC100750224","rna_LOC101739648","rna_LOC101741268","rna_LOC101741704","rna_LOC101744414","rna_LOC101747109","rna_Optix","rna_Sox102F","rna_Sox21b","rna_dmrt99B","rna_elB","rna_hbn","rna_hth","rna_opa","rna_scro","rna_so","rna_toy","rna_ey","rna_slp2","rna_LOC101747128","rna_D","rna_tll")

FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/bombyx/umap.png", units = "cm", width = 60, height = 60, bg = "white")

####tribolium####
NbS <- readRDS("/tribolium_medulla_NBs.rds")
tfs = read.csv("/Tcast_TFS.csv",  header = F)
res_list <- readRDS("/tribolium_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1


fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime

quantile_matrix = get_quantile_matrix(pseudotime, intervals = 20)
# its important to use sce not NB (seurat object), otherwise it is accesing the seurat integrated object
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("D","Optix","Sox14","Sox21b","TcasGA2-TC011697","TcasGA2-TC014730","TcasGA2-TC031415","TcasGA2-TC032830","TcasGA2-TC032832","btd","dmrt99B","dan","hbn","erm","hth","opa","slp1","tll","so")
keep= c("ab","nerfin-1")
keep= c("TcasGA2-TC004682")
keep=c("scro")
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/tribolium/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_scro","rna_D","rna_Optix","rna_Sox14","rna_Sox21b","rna_TcasGA2-TC011697","rna_TcasGA2-TC014730","rna_TcasGA2-TC031415","rna_TcasGA2-TC032830","rna_TcasGA2-TC032832","rna_btd","rna_dmrt99B","rna_dan","rna_hbn","rna_erm","rna_hth","rna_opa","rna_slp1","rna_tll","rna_so","rna_ab","rna_dan")
keep=c("rna_scro")
FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/tribolium/umap.png", units = "cm", width = 15, height = 15, bg = "white")
FeaturePlot(NbS, features = dpn, pt.size = 2, combine=T)

#### nasonia ####
NbS <- readRDS("/nasonia_medulla_NBs.rds")
tfs = read.csv("/Nvit_TFS.csv",  header = F)
res_list <- readRDS("/nasonia_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]
keep= c("tll","Acf","CG42741","Dll","ey","FoxP","grh","grn","hbn","LOC100117267","LOC100118340","LOC100122379","LOC100119062","opa","slp2","LOC100679002","LOC100678366","pdm3","scro","Sox14","Sox21b")
keep= c("LOC100679231","cas")
keep= c("PRAS40")

for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/nasonia/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_hth","rna_tll","rna_Acf","rna_CG42741","rna_Dll","rna_ey","rna_FoxP","rna_grh","rna_grn","rna_hbn","rna_LOC100117267","rna_LOC100118340","rna_LOC100122379","rna_LOC100119062","rna_opa","rna_slp2","rna_LOC100679002","rna_LOC100678366","rna_pdm3","rna_scro","rna_Sox14","rna_Sox21b","rna_LOC100679231","rna_cas","rna_PRAS40","rna_slp2")
FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/nasonia/umap.png", units = "cm", width = 60, height = 60, bg = "white")

#### aedes ####
NbS <- readRDS("/aedes_medulla_NBs.rds")
tfs = read.csv("/Aaeg_TFS.csv",  header = F)
res_list <- readRDS("/aedes_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("LOC5565423","LOC5568335","LOC5568335","LOC5570655","LOC5573206","LOC5573206","LOC5580130","LOC5566890","LOC5566891","LOC5565439","tup","cas","CG6769","CG13367","CG17202","CG42741","D","dmrt99B","erm","ey","hbn","LOC5573207","LOC5575090","LOC5575091","opa","scro","slp2","Sox15","Sp1","tll")
keep=c("CG13367","CG17202","CG42741","D","dmrt99B","erm","ey","hbn","LOC5573207","LOC5575090","LOC5575091","opa","scro","slp2","Sox15","Sp1","tll")
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/aedes/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_LOC5565423","rna_LOC5568335","rna_LOC5568335","rna_LOC5570655","rna_LOC5573206","rna_LOC5573206","rna_LOC5580130","rna_LOC5566890","rna_LOC5566891","rna_LOC5565439","rna_tup","rna_cas","rna_CG6769","rna_CG13367","rna_CG17202","rna_CG42741","rna_D","rna_dmrt99B","rna_erm","rna_ey","rna_hbn","rna_LOC5573207","rna_LOC5575090","rna_LOC5575091","rna_opa","rna_scro","rna_slp2","rna_Sox15","rna_Sp1","rna_tll")

FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("/FigureS9/aedes/umap.png", units = "cm", width = 60, height = 60, bg = "white")

#### cricket ####
NbS <- readRDS("/gryllus_medulla_NBs.rds")
tfs = read.csv("/Gbim_TFS.csv",  header = F)
res_list <- readRDS("/gryllus_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

# run slingshot and pseudotimeDE
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")

rd <- NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = "UMAP", clusterLabels = "cl")

# invert pseudotime (using the slingshot result)
pt_inv <- max(colData(fit_ori_umap)$slingPseudotime_1, na.rm = TRUE) -
  colData(fit_ori_umap)$slingPseudotime_1

LPS_ori_tbl_umap <- tibble(
  cell = colnames(fit_ori_umap),
  pseudotime = rescale(pt_inv, na.rm = TRUE)
)
Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("dve","edl","ey","GBIM-00856","GBIM-01738","GBIM-04903","GBIM-05211","GBIM-06054","GBIM-10435","GBIM-12506","GBIM-16025","GBIM-16201","hbn","Sp1","tup")

for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/gryllus/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_dve","rna_edl","rna_ey","rna_GBIM-00856","rna_GBIM-01738","rna_GBIM-04903","rna_GBIM-05211","rna_GBIM-06054","rna_GBIM-10435","rna_GBIM-12506","rna_GBIM-16025","rna_GBIM-16201","rna_hbn","rna_Sp1","rna_tup")

FeaturePlot(NbS, features = keep, pt.size = 2, combine=T, reduction = "umap")
ggsave("/FigureS9/gryllus/umap.png", units = "cm", width = 60, height = 60, bg = "white")

####cloeon####
NbS <- readRDS("/cloeon_medulla_NBs.rds")
tfs = read.csv("/Cdipt_TFS.csv",  header = F)
res_list <- readRDS("/cloeon_res_list.rds")

NbS[["RNA"]] <- as(NbS[["RNA"]], Class="Assay")
NB_sce <- as.SingleCellExperiment(NbS, assay = "RNA")
# run slingshot and pseudotimeDE
rd = NbS@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

Nb <- NbS
pseudotime = LPS_ori_tbl_umap$pseudotime

sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

keep= c("CLODIP-2-CD00329","CLODIP-2-CD02704","CLODIP-2-CD07852","CLODIP-2-CD08459","CLODIP-2-CD09832","CLODIP-2-CD10553","CLODIP-2-CD10783","CLODIP-2-CD11133","CLODIP-2-CD12980","CLODIP-2-CD14302","CLODIP-2-CD15157","CLODIP-2-CD15283","CLODIP-2-CD15498","CLODIP-2-CD15498","CLODIP-2-CD15498","elB","grh","hbn","luna","noc","opa","slp2","CLODIP-2-CD05887")

for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  plot_info <- ggplot_build(p)
  
  p <- p +
    geom_point(alpha = 0.8, size = 1.2) +
    annotate(
      geom = 'text',
      x = 0.25,
      y = max(plot_info$data[[1]]$y),
      label = paste('p =', res$para.pv),
      color = "red",
      size = 4
    ) +
    theme(
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = paste0("/FigureS9/cloeon/", gene, ".png"),
    plot = p,
    width = 15,
    height = 15,
    units = "cm",
    dpi = 300,
    bg = "white"
  )
  
  
}

keep= c("rna_CLODIP-2-CD00329","rna_CLODIP-2-CD02704","rna_CLODIP-2-CD07852","rna_CLODIP-2-CD08459","rna_CLODIP-2-CD09832","rna_CLODIP-2-CD10553","rna_CLODIP-2-CD10783","rna_CLODIP-2-CD11133","rna_CLODIP-2-CD12980","rna_CLODIP-2-CD14302","rna_CLODIP-2-CD15157","rna_CLODIP-2-CD15283","rna_CLODIP-2-CD15498","rna_CLODIP-2-CD15498","rna_CLODIP-2-CD15498","rna_elB","rna_grh","rna_hbn","rna_luna","rna_noc","rna_opa","rna_slp2","rna_CLODIP-2-CD05887")

FeaturePlot(NbS, features = keep, pt.size = 2, combine=T)
ggsave("FigureS9/cloeon/umap.png", units = "cm", width = 60, height = 60, bg = "white")
