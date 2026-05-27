# This R file contains all the information for generation of Figure S7
#orthologs_assigned.txt is based on Table_S2 and provided in the current github page.

library(ggplot2)

df <- read_tsv("orthologs_assigned.txt", col_names = T)

names(df)[1] <- "Category"

df_long <- df %>%
  pivot_longer(-Category, names_to = "Species", values_to = "Count")

df_long$Species <- factor(df_long$Species,
                          levels = c("DROVI", "MUSDO", "AEDAE", "BOMMO", 
                                     "TRICA", "NASVI", "GRYBI", "CLOEON", "DAPPU"))

df_long$Category <- factor(df_long$Category,
                           levels = c("not_found", "many_to_many", 
                                      "many_to_one", "one_to_many", "one_to_one"))


okabe <- palette.colors(palette = "Okabe-Ito")

mycols <- okabe[c(9,8,4,3,6)]  # pick distinct ones

p<-ggplot(df_long, aes(x = Species, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_classic() +
  labs(x = "", y = "Number of genes", fill = "Orthology class") +
  scale_fill_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("orthology_vs_Drosophila.png", units = "cm", width = 20, height = 15, bg = "white")
