library(ggplot2)
library(tidyverse)

setwd("~/Seafile/Konstantinides_lab/Konstantinides_lab_papers/drafts/Konstantina_paper/Figures/FigureS2/GeneExt_extensions/")

Dvir<-read.table("tmp_Virilis_merged_subsample100m_orphans/extensions.tsv")
Musca<-read.table("tmp_Musca_merged_subsample100m_orphans/extensions.tsv")
Aedes<-read.table("tmp_Aedes_merged_subsample100m_orphans/extensions.tsv", header = F)
Bombyx<-read.table("tmp_Bombyx_merged_subsample100m_orphans/extensions.tsv")
Tribolium<-read.table("tmp_Tribolium_merged_GCA_annot_subsample100m_orphans/extensions.tsv")
Nasonia<-read.table("tmp_Nasonia_merged_subsample100m_orphans/extensions.tsv")
Gryllus<-read.table("tmp_Gryllus_merged_subsample100m_orphans/extensions.tsv")
Cloeon<-read.table("tmp_Cloeon_merged_subsample100m_orphans/extensions.tsv")


boxplot(Dvir$V3, Musca$V3, Aedes$V3, Bombyx$V3, Tribolium$V3, Nasonia$V3, Gryllus$V3, Cloeon$V3)


df <- tibble(
  Species = rep(c("Dvir","Musca","Aedes","Bombyx","Tribolium","Nasonia","Gryllus","Cloeon"), 
                times = c(length(Dvir$V3), length(Musca$V3), length(Aedes$V3),
                          length(Bombyx$V3), length(Tribolium$V3), length(Nasonia$V3),
                          length(Gryllus$V3), length(Cloeon$V3))),
  Value = c(Dvir$V3, Musca$V3, Aedes$V3, Bombyx$V3, Tribolium$V3, Nasonia$V3, Gryllus$V3, Cloeon$V3)
)

df$Species <- factor(df$Species, 
                     levels = c("Dvir","Musca","Aedes","Bombyx","Tribolium","Nasonia","Gryllus","Cloeon"))

ggplot(df, aes(x = Species, y = Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
#  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  theme_classic() +
  labs(x = "", y = "Your measurement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p<-ggplot(df, aes(x = Species, y = Value)) +
  geom_violin(trim = FALSE, fill = "#AA87DE", color = "grey30", width = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white") +
  theme_classic() +
  labs(x = "", y = "Your measurement") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("GeneExt_extensions6.png", units = "cm", width = 20, height = 15, bg = "white")
