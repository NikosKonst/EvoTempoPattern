# This R script performs trajectory analysis of species-specific medulla NB clusters
# using slingshot, followed by dynamic transcription factor (TF) expression analysis
# with PseudotimeDE.
#
# Input clusters are obtained from the corresponding species_OL_medulla_NBs.R workflow.
#
# An additional script, plot_selection.R, is required for selecting TFs with
# temporal expression patterns.
#
# The script also includes code to:
# - plot temporal TF expression profiles,
# - visualize any TF present in the PseudotimeDE output (res_list.rds),
# - and annotate plots with the PseudotimeDE parametric p-value.
#
# NbS inputs are provided in the GEO, accession number GSE333138
# tfs lists for each species are coming from Table_S2
# res_lists are provided in the current github page
#
# Workflow overview:
#
# STEP1  Load species-specific medulla_NBs.rds objects
# STEP2  Prepare SingleCellExperiment objects and infer pseudotime using slingshot
# STEP3  Generate random subsamples for PseudotimeDE
# STEP4  Load species-specific TF gene lists
# STEP5  Run PseudotimeDE and save results to res_list.rds
# STEP6  Run functions from plot_selection.R to identify temporally expressed TFs
# STEP7  Filter TFs based on genes whose expression reaches near-zero levels across pseudotime and expression criteria
# STEP8  Plot selected genes and annotate plots with PseudotimeDE parametric p-values

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

# Run the workflow separately for each species.
# Load only one medulla_NBs.rds file at a time.
#STEP1 input:
NbS<- readRDS("/melanogaster_medulla_NBs.rds")
NbS<- readRDS("/virilis_medulla_NBs.rds")
NbS<- readRDS("/musca_medulla_NBs.rds")
NbS<- readRDS("/aedes_medulla_NBs.rds")
NbS<- readRDS("/bombyx_medulla_NBs.rds")
NbS<- readRDS("/tribolium_medulla_NBs.rds")
NbS<- readRDS("/nasonia_medulla_NBs.rds")
NbS<- readRDS("/gryllus_medulla_NBs.rds")
NbS<- readRDS("/cloeon_medulla_NBs.rds")

# only one medulla_NBs.rds file at a time for the whole workflow
# make seurat_object as SingleCellExperiment
Nb<-NbS
#if required for some objects first run Nb[["RNA"]] <- as(Nb[["RNA"]], Class="Assay"), and then NB_sce <- as.SingleCellExperiment(Nb, assay = "RNA")
Nb[["RNA"]] <- as(Nb[["RNA"]], Class="Assay")

NB_sce <- as.SingleCellExperiment(Nb, assay = "RNA")

# run slingshot and pseudotimeDE, for pseudotimeDE check https://github.com/SONGDONGYUAN1994/PseudotimeDE
#STEP2
rd = Nb@reductions$umap@cell.embeddings
reducedDims(NB_sce) <- SimpleList(UMAP = rd)
colData(NB_sce)$cl <- 1

#run slingshot
#### Pseudotime assessment ####
fit_ori_umap <- slingshot(NB_sce, reducedDim = 'UMAP', clusterLabels = "cl")
LPS_ori_tbl_umap <- tibble(cell = colnames(NB_sce), pseudotime = rescale(colData(fit_ori_umap)$slingPseudotime_1))

#STEP3
set.seed(123)
## Set the cores for parallelization.
options(mc.cores = 6)

## Number of subsmaples
n = 100

## Generate random subsamples for PseudotimeDE
LPS_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(NB_sce)[2]), size = 0.8*dim(NB_sce)[2], replace = FALSE)
})

LPS_sub_tbl_umap <- mclapply(LPS_index, function(x, sce) {
  
  sce <- sce[, x]
  seu <- as.Seurat(sce)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu<- RunPCA(seu)
  seu = RunUMAP(seu, dims = 1:50)
  rd <-  seu@reductions$umap@cell.embeddings[, 1:2]
  reducedDims(sce) <- SimpleList(UMAP = rd)
  
  fit <-  slingshot(sce, reducedDim = 'UMAP')
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, LPS_ori_tbl_umap, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = NB_sce)

#STEP4
#provide the list of genes of interest, here lists with all TFs from each species. Available in Table_S2. Use the list respective to the species name in STEP1 
tfs = read.csv("/Dvir_TFS.csv",  header = F)
tfs = read.csv("/Mdom_TFS.csv",  header = F)
tfs = read.csv("/Aaeg_TFS.csv",  header = F)
tfs = read.csv("/Bmori_TFS.csv",  header = F)
tfs = read.csv("/Tcast_TFS.csv",  header = F)
tfs = read.csv("/Nvit_TFS.csv",  header = F)
tfs = read.csv("/Gbim_TFS.csv",  header = F)
tfs = read.csv("/Cdipt_TFS.csv",  header = F)
tfs = read.csv("/Dmel_TFS.csv",  header = F) # available from https://github.com/NikosKonst/larva_scSeq2022 Flybase_TFs.csv

#STEP5
#Generate the output of pseudotimeDE for each TF and store it in the rest_list.rds file
#### PseudotimeDE output ####
#input
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]
res_list = list()
for(i in sample_genes){
  res<- PseudotimeDE::runPseudotimeDE(gene.vec = i,
                                      ori.tbl = LPS_ori_tbl_umap,
                                      sub.tbl = LPS_sub_tbl_umap[1:100], ## To save time, use 100 subsamples
                                      mat = NB_sce, ## You can also use a matrix or SeuratObj as the input
                                      model = "nb",mc.cores = 6)
  res_list[[i]] = res
  print(i)
}
# the following files are provided in the current github page
saveRDS(res_list, "/virilis_res_list.rds")
saveRDS(res_list, "/musca_res_list.rds")
saveRDS(res_list, "/aedes_res_list.rds")
saveRDS(res_list, "/bombyx_res_list.rds")
saveRDS(res_list, "/tribolium_res_list.rds")
saveRDS(res_list, "/nasonia_res_list.rds")
saveRDS(res_list, "/gryllus_res_list.rds")
saveRDS(res_list, "/cloeon_res_list.rds")
saveRDS(res_list, "/melanogaster_res_list.rds")

#### select TFs with temporal expression ####
#STEP6
# inputs
pseudotime = LPS_ori_tbl_umap$pseudotime
quantile_matrix = get_quantile_matrix(pseudotime, intervals = 20)
# Use the SingleCellExperiment object (NB_sce) rather than the Seurat object (Nb),
# otherwise the integrated Seurat assay may be accessed instead of raw counts.
sample_genes = tfs$V1[tfs$V1%in%rownames(NB_sce)]

# run functions in plot_selection.R, call them here
reach0 = reach0filter(sample_genes, quantile_matrix$original,count_matrix= Nb@assays$RNA$counts)
gene_exp_filter = gene_exp_filter(sample_genes, assay(NB_sce, 'counts'))

#STEP7
#ensure temporality (reaching 0) and expression by keeping genes that pass both filters
keep = gene_exp_filter=='pass' & reach0=='pass' 
keep = sample_genes[keep]

#STEP8
dir.create(data.dir =("/keep"))

# Save plots for the 'keep' genes
### try generating the plots here directlu, in order to save memory ###
for (i in keep) {
  gene = i
  res = res_list[[gene]]
  
  # Generate the plot
  p <- PseudotimeDE::plotCurve(gene.vec = gene,
                               ori.tbl = LPS_ori_tbl_umap,
                               mat = as.matrix(Nb@assays$RNA$counts),
                               model.fit = res$gam.fit)
  
  # print parametric p-value annotation from the res_list file on each plot
  plot_info = ggplot_build(p)
  p = p + annotate(geom = 'text', x = c(0.25), y = c(max(plot_info$data[[1]]$y)), 
                   label = c(paste('para=', res$para.pv)), color = "red")
  
  # Save the plot
  ggsave(p, filename = paste0("/keep/", gene, '.pdf'),
         width = 8.27, height = 11.69, units = "in")
}

