# This R file contains all the information for generation of Figure 4B

library(ggplot2)

# -------------------------
# Data (updated with 4 categories) based on events described in the tree of Figure 4A, 
# early is considered whatever would be in the window of hth-erm-opa of melanogaster, middle the ey-hbn-slp and the rest are considered as late 
# -------------------------
df <- data.frame(
  period = rep(c("Early", "Mid", "Late"), each = 4),
  cat    = rep(c("#AADD87", "#008000", "#800000", "#DE8787"), times = 3),
  n      = c(
    7, 3, 5, 2,   # Early
    1, 2, 2, 0,   # Mid
    6, 15, 3, 0   # Late
  )
)

# Factor order
df$period <- factor(df$period, levels = c("Early", "Mid", "Late"))
df$cat    <- factor(df$cat, levels = c("#AADD87", "#008000", "#800000", "#DE8787"))

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(x = period, y = n, fill = cat)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "#AADD87"    = "#AADD87",
      "#008000" = "#008000",
      "#800000" = "#800000",
      "#DE8787" = "#DE8787"   # new colour
    ),
    name = NULL
  ) +
  labs(x = NULL, y = "n") +
  theme_classic(base_size = 13)

p

ggsave("/plots/Fig4B_quantification_changes.png", units = "cm", width = 15, height = 15, bg = "white")


