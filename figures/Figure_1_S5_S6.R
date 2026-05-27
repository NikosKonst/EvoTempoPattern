# This R file contains all the information for generation of Figure1 and Supplementary Figure S5 and S6.
# whole optic lobe data as inputs are provided in the GEO, accession number GSE333138
# requires function LayerdFeaturePlot.R 

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

Virilis <- readRDS("/virilis_OL.rds")
Musca <-readRDS("/musca_OL.rds")
Aedes <- readRDS("/aedes_OL.rds")
Bombyx <- readRDS("/bombyx_OL.rds")
Tribolium <- readRDS("/tribolium_OL.rds")
Nasonia <-readRDS("/nasonia_OL.rds")
Cricket <- readRDS("/gryllus_OL.rds")
Cloeon <- readRDS("/cloeon_OL.rds")
Cloeon<-FindClusters(object = Cloeon, resolution = 8)


#### Figure 1E
#run the LayeredFeaturePlot function
#Virilis
LayeredFeaturePlot(Virilis, features = c("rna_mira","rna_nSyb","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 10,reduction = "umap")

ggsave("/plots/1E_virilis_mira_nsyb_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Musca
#LOC101897603=repo
LayeredFeaturePlot(Musca, features = c("rna_mira","rna_brp","rna_LOC101897603"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 2, reduction = "umap")
ggsave("/plots/1E_Musca_mira_brp_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Aedes
LayeredFeaturePlot(Aedes, features = c("rna_mira","rna_nSyb","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 10,reduction = "umap")
ggsave("/plots/1E_Aedes_mira_nsyb_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Bombyx
LayeredFeaturePlot(Bombyx, features = c("rna_mira","rna_brp","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 10,reduction = "umap")
ggsave("/plots/1E_Bombyx_mira_brp_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Tribolium
#TcasGA2-TC007824=mira
LayeredFeaturePlot(Tribolium, features = c("rna_TcasGA2-TC007824","rna_nSyb","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 5,reduction = "umap")
ggsave("/plots/1E_tribolium_mira_nsyb_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Nasonia
LayeredFeaturePlot(Nasonia, features = c("rna_mira","rna_nSyb","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 10, reduction = "umap")
ggsave("/plots/1E_nasonia_mira_nsyb_repo.png", units = "cm", width = 45, height = 30, bg = "white")

#Cricket
#GBIM-06848mira GBIM-15830fne, eaat1-axo:glia markers
LayeredFeaturePlot(Cricket, features = c("rna_GBIM-06848","rna_GBIM-15830","rna_Eaat1","rna_axo"), 
                   colors = c("#F0E442", "#0072B2","#CC79A7","#D55E00"),pt.size = 2 , 
                   pt.brightness = 5,reduction = "umap.cca")
ggsave("/plots/1E_cricket_mira_fne_eaat1_axo.png", units = "cm", width = 45, height = 30, bg = "white")

#Cloeon
LayeredFeaturePlot(Cloeon, features = c("rna_dpn","rna_nSyb","rna_repo"), 
                   colors = c("#F0E442", "#0072B2", "#CC79A7"),pt.size = 2 , 
                   pt.brightness = 5,reduction = "umap")
ggsave("/plots/1E_Cloeon_dpn_nsyb_repo.png", units = "cm", width = 45, height = 30, bg = "white")


#### Figure 1F ####

#virilis
DimPlot(Virilis,label = T, raster=F)+NoLegend()
lamina_cells<-c("55","38", "44", "34","58","60","36","35","83","7")
lobula_plate_cells<-c("31", "22", "15","25")
glia<-c("54","68","84","77","37","88","43","19","47","69","87")
central_brain <-c("10","33")
unknown_cells<-c("10", "33")
neuroblasts <- c("16")

#Musca
DimPlot(Musca,label = T, raster=F)+NoLegend()
lamina_cells<-c("72","74", "18","6")
lobula_plate_cells<-c("83","54","41","34","30","67","26","46","59")
glia<-c("32","79","21","25","73","77","4","64","28")
unknown_cells<-c("11", "85", "61","80", "58", "24", "22","15")
central_brain <-c("51")
neuroblasts <- c("20")

#Aedes
DimPlot(Aedes,label = T, raster=F)+NoLegend()
lamina_cells<-c("42","64", "63", "62","19")
lobula_plate_cells<-c("14", "29", "46","8","30","33")
glia<-c("45","11","22","77","4","65","15","0","39","44")
unknown_cells<-c("70", "69", "18","24", "28", "17", "1","57","2","6","32","10","9","20","66","23","13")
central_brain <-c("3")
neuroblasts <- c("21")

#Bombyx
DimPlot(Bombyx_orig,label = T, raster=F) + NoLegend()
lamina_cells<-c("50","13")
lobula_plate_cells<-c("61", "32")
glia<-c("33","17","51","54","11","14","38","71","63","41","18","20","59","31","10","19","22","6","72","1","47","42","23","41","56")
neuroblasts <- c("48")
central_brain <-c("8","12","16","49","46","52")

#Tribolium
DimPlot(Tribolium,label = T, raster=F)
lamina_cells<-c("38","44", "39", "45","27","78","51")
lobula_plate_cells<-c("15", "29", "94")
glia<-c("24","90","124","47","95","104","65","71","68","16","20","87","14","13","84","63","89","108","56","58","106")
neuroblasts <- c("25","61")
unknown_cells<-c("0", "1","10", "9","17", "18")
medulla<-c("50","54","38","44","39","45","26","34","31","32","29","94","30","48","77","37","55","4","12","43","97","40","103","8","125","81","52","41","49","92","91","62","75","100","110","127","128","107","113","111","116","22","32","28","73","57","59","53","98","72","114","33","7","35","82","80","74","123","120","99")
#central_brain <-c("3")

#Nasonia
lamina_cells<-c("41", "52", "45", "44")
lobula_plate_cells<-c("15", "16", "26")
glia<-c("31", "6", "7", "8", "28", "50","53")
unknown_cells<-c("37", "47", "58", "57", "27", "46","14","13","39","20","18","40","35","30")
central_brain<-c("43","12")
neuroblasts <- c("23")

#Cricket
DimPlot(Cricket,label = T, raster=F) + NoLegend()
FeaturePlot(Cloeon, features=c("rna_pros"))
lamina_cells<-c("75","42","49","66","32","85","46","45","98","59","48")
glia<-c("79","61","65","39","55","60","36","74","44","91","23","63","30","18","7","17","5","13","20","37","57","51","70","68")
unknown_cells<-c("10","12","9","3","33","103","29","27","1","4","2","0","76","78","19","28","26")
neuroblasts <- c("56")
central_brain <-c("6","8","11","16")

#Cloeon from resolution 8 cluster2 Nbs
DimPlot(Cloeon,label = T, raster=F) + NoLegend()
lamina_cells<-c("55","103","90","31","72","44","33","26","36","18","99","129")
lobula_plate_cells<-c("19","12","39","115","66","46","7","105","71","57","22","38","56","42","107")
glia<-c("64","84","47","76","133","4","3","8","79","10","48","113","68","27","82","135","96","121")
unknown_cells<-c("29","126","119")
neuroblasts <- c("2")
central_brain <-c("1","41")


#example for Musca
Neuropil<-as.character(Musca@active.ident)

for (i in 1:length(Neuropil)) {
  if (Neuropil[i]%in%lamina_cells) {
    Neuropil[i]="lamina"
  }
  else { if (Neuropil[i]%in%lobula_plate_cells) {
    Neuropil[i]="lobula_plate"
  }
    else { if (Neuropil[i]%in%glia) {
      Neuropil[i]="glia"
    }
      else { if (Neuropil[i]%in%unknown_cells) {
        Neuropil[i]="unknown"
      }
        else { if (Neuropil[i]%in%neuroblasts) {
          Neuropil[i]="medulla_NBs"
        }
          else { if (Neuropil[i]%in%central_brain) {
            Neuropil[i]="central_brain"
          }
            else {
              Neuropil[i]="medulla"
            }
          }
        }
      }
    }
  }}

Neuropil1<-as.factor(Neuropil)

Musca$Neuropil<-Neuropil1

DimPlot(
  Musca,
  group.by = "Neuropil",
  cols = c(
    "lamina"        = "#CC79A7",  # yellow #56B4E9
    "lobula_plate"  = "#F0E442",  # blue
    "glia"          = "#009E73",  # green
    "unknown"       = "#999999",  # grey
    "central_brain" = "#8E24AA",  # magent
    "medulla_NBs"   = "#56B4E9",  # orange/red" #56B4E9
    "medulla"       = "#0072B2"   # purple #0072B2
  ), pt.size = 2 ,
  raster = FALSE
)+
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) #+ NoLegend()

ggsave("/plots/1F_virilis_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Musca_colour_neuropil_with_ey.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Aedes_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Bombyx_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Nasonia_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Tribolium_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Cloeon_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_Cricket_colour_neuropil.png", units = "cm", width = 45, height = 25, bg = "white")
ggsave("/plots/1F_legend.png", units = "cm", width = 45, height = 25, bg = "white")

#### Figure S5 ####
FeaturePlot(Virilis, features = c("rna_shg","rna_dpn","rna_wor","rna_ase","rna_repo"), combine=T) 
ggsave("/plots/S5_virilis_shg_dpn_wor_ase_repo.png", units = "cm", width = 40, height = 30, bg = "white")

#LOC101890262 wor, LOC101897603 repo
FeaturePlot(Musca, features = c("rna_shg","rna_dpn","rna_LOC101890262","rna_ase","rna_elav","rna_LOC101897603"), min.cutoff = c("1", "q0", "q0","q0","q0","q0"),combine=T) 
ggsave("/plots/S5_Musca_shg_dpn_wor_ase_repo_with_EY.png", units = "cm", width = 40, height = 30, bg = "white")

#LOC5570443 wor
FeaturePlot(Aedes, features = c("rna_shg","rna_dpn","rna_LOC5570443","rna_ase","rna_repo"), combine=T) 
ggsave("/plots/S5_aedes_shg_dpn_wor_ase_repo.png", units = "cm", width = 40, height = 30, bg = "white")

#LOC101744234 wor
FeaturePlot(Bombyx, features = c("rna_shg","rna_dpn","rna_LOC101744234,","rna_ase","rna_repo"), min.cutoff = c("q0", "q0", "q0","q0","1"), combine=T) 
ggsave("/plots/S5_bombyx_shg_dpn_wor_ase_repo.png", units = "cm", width = 40, height = 30, bg = "white")

#rna_TcasGA2-TC014474:wor, ase:TcasGA2-TC008437
FeaturePlot(Tribolium, features = c("rna_shg","rna_dpn","rna_TcasGA2-TC014474","rna_TcasGA2-TC008437","rna_repo"), min.cutoff = c("1", "q0", "q0","q0","1"),combine=T) 
ggsave("/plots/S5_tribolium_shg_dpn_wor_ase_repo.png", units = "cm", width = 40, height = 30, bg = "white")

#LOC100114879 wor, no ase 
FeaturePlot(Nasonia, features = c("rna_shg","rna_dpn","rna_LOC100114879","rna_repo"), min.cutoff = c("q0", "q0", "q0","q0"),combine=T) 
ggsave("/plots/S5_nasonia_shg_dpn_wor_repo.png", units = "cm", width = 40, height = 30, bg = "white")

#GBIM-06054: dpn
FeaturePlot(Cricket, features = c("rna_shg","rna_GBIM-06054"), min.cutoff = c("q0", "q0"),combine=T,raster=F) 
ggsave("/plots/S5_cricket_shg_dpn.png", units = "cm", width = 40, height = 30, bg = "white")

#CLODIP-2-CD00134 shg
FeaturePlot(Cloeon, features = c("rna_CLODIP-2-CD00134","rna_dpn","rna_repo"),min.cutoff = c("q0", "q0", ,"1"),combine=T,raster=F) 
ggsave("/plots/S5_cloeon_shg_dpn_elav_repo.png", units = "cm", width = 40, height = 30, bg = "white")

FeaturePlot(Cricket, features = c("rna_shg","rna_Eaat1","rna_axo"), min.cutoff = c("q0", "q0", "q0"),combine=T,raster=F) 
ggsave("/plots/S5_cricket_shg_eaat1_ax0.png", units = "cm", width = 40, height = 30, bg = "white")


#### Figure S6 ####
# gcm LOC100114913
FeaturePlot(Nasonia, features = c("rna_LOC100114913","rna_eya","rna_sim","rna_tll","rna_dac","rna_acj6")) 
ggsave("plots/S6_nasonia_gcm_dac_sim_tll_eya_acj6.png", units = "cm", width = 40, height = 30, bg = "white")

# dac FBpp0080433, MDOA005900-PA LOC101894733 (info uniprot), eya FBpp0078964 MDOA009994-PA (LOC101893172 vectorbase NOT EXPR),MDOA001006-PB (LOC101894691 EXPRESSED uniprot)
FeaturePlot(Musca, features = c("rna_gcm","rna_LOC101894691","rna_sim","rna_tll","rna_dac","rna_acj6")) 
ggsave("/plots/S6_Musca_lamina_lob_pl.png", units = "cm", width = 40, height = 30, bg = "white")

#AAEL000631-PB, uniprot-> 5564850=gcm
FeaturePlot(Aedes, features = c("rna_LOC5564850","rna_eya","rna_sim","rna_tll","rna_dac","rna_acj6")) 
ggsave("/plots/S6_aedes_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")

#gcm:TcasGA2-TC014730, sim:TcasGA2-TC016205, acj6:TcasGA2-TC003196
FeaturePlot(Tribolium, features = c("rna_TcasGA2-TC014730","rna_eya","rna_TcasGA2-TC016205","rna_tll","rna_dac","rna_TcasGA2-TC003196")) 
ggsave("/plots/S6_tribolium_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")

#gcm_LOC101742696, acj6:100750224
FeaturePlot(Bombyx, features = c("rna_LOC101742696","rna_eya","rna_sim","rna_tll","rna_dac","rna_LOC100750224")) 
ggsave("/plots/S6_bombyx_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")

FeaturePlot(Virilis, features = c("rna_gcm","rna_eya","rna_sim","rna_tll","rna_dac","rna_acj6")) 
ggsave("/plots/S6_virilis_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")

#gcm_CLODIP_2_CD06246 tll:CLODIP-2-CD05887
FeaturePlot(Cloeon, features = c("rna_CLODIP-2-CD06246","rna_eya","rna_sim","rna_CLODIP-2-CD05887","rna_dac","rna_acj6")) 
ggsave("/plots/S6_Cloeon_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")

#tll:GBIM-04903
FeaturePlot(Cricket, features = c("rna_gcm","rna_eya","rna_sim","rna_GBIM-04903","rna_acj6"), raster = F) # no dac
ggsave("/plots/S6_Cricket_lamina_lob.png", units = "cm", width = 40, height = 30, bg = "white")
