library(ggplot2)
library(rlang)
library(ggrepel)
library(cowplot)
LayeredFeaturePlot <- function(object, features = NULL, assay = NULL, 
                               reduction = "tsne", cells.use = NULL, pt.size = 1, 
                               colors = NULL, min.cutoff = NULL, max.cutoff = NULL, 
                               fig.ratio = 6, label = FALSE, repel = FALSE, bold = FALSE, 
                               background.color = "white", pt.brightness = 10) {
  
  ##### Internal helpers #####
  exp_encolor <- function(data, features, pt.brightness) {
    data <- data[, features, drop = FALSE]
    norm_data <- apply(data, 2, function(x) {
      if (max(x, na.rm = TRUE) == 0 || all(is.na(x))) return(rep(0, length(x)))
      norm_to <- max(x, na.rm = TRUE) / (100 - pt.brightness)
      return(x / norm_to)
    })
    ldata <- reshape(data.frame(norm_data), times = features, varying = features, 
                     timevar = "type", idvar = "row.names", v.names = "l", 
                     direction = "long", new.row.names = seq(nrow(data) * length(features)))
    if (is.null(colors)) {
      type_hue <- (90 + (360 / length(features)) * (seq_along(features) - 1)) %% 360
    } else {
      type_hue <- sapply(colors, function(x) rgb2hsv(col2rgb(x))["h", ] * 360)
    }
    names(type_hue) <- features
    ldata$h <- type_hue[ldata$type]
    ldata$s <- ifelse(ldata$l > 0, 100, 0)
    ldata$l <- ldata$l + pt.brightness
    ldata$color <- apply(ldata[, c("h", "s", "l")], 1, function(x) {
      h <- min(max(as.numeric(x[1]), 0), 360)
      s <- min(max(as.numeric(x[2]), 0), 100)
      l <- min(max(as.numeric(x[3]), 0), 100)
      hsv(h/360, s/100, l/100)
    })
    ldata$r <- strtoi(substr(ldata$color, 2, 3), 16)
    ldata$g <- strtoi(substr(ldata$color, 4, 5), 16)
    ldata$b <- strtoi(substr(ldata$color, 6, 7), 16)
    min.col <- hsv(1, 0, pt.brightness/100)
    for (ch in c("r", "g", "b")) ldata[[ch]][ldata$color == min.col] <- ldata[[ch]][ldata$color == min.col]/2
    pal <- data.frame(
      r = sapply(tapply(ldata$r, ldata$row.names, sum), function(x) min(x, 255)),
      g = sapply(tapply(ldata$g, ldata$row.names, sum), function(x) min(x, 255)),
      b = sapply(tapply(ldata$b, ldata$row.names, sum), function(x) min(x, 255))
    )
    apply(pal, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  }
  
  # ✅ Updated make_legend block
  make_legend <- function(data, features, colors = NULL, min.col) {
    dummy <- data[, features, drop = FALSE]
    dummy[nrow(dummy) + 1, ] <- 0
    dummy$x <- c(1:nrow(dummy))
    dummy$y <- c(1:nrow(dummy))
    
    if (is.null(colors)) {
      type_hue <- rep(90, length(features))
      angle_step <- 360 / length(features)
      times <- seq(length(features)) - 1
      type_hue <- (type_hue + (angle_step * times)) %% 360
    } else {
      type_hue <- sapply(colors, function(x) rgb2hsv(col2rgb(x))["h", ] * 360)
    }
    
    names(type_hue) <- features
    type_col <- sapply(type_hue, function(x) hsv(x/360, s = 1, v = 1))
    legends <- list()
    
    for (i in features) {
      gene <- sym(i)
      legends[[i]] <- get_legend(
        ggplot(dummy, aes(x = x, y = y, color = !!gene)) + 
          geom_point() +
          scale_color_gradient(low = min.col, high = type_col[i]) +
          labs(color = NULL) +  # ← remove gene name
          theme(legend.title = element_blank(),
                legend.text = element_text(size = 10))
      )
    }
    
    # align legends 
    legend_width <- max(sapply(legends, function(x) {
      as.numeric(grid::convertWidth(sum(x$widths), "cm"))
    }))
    
    for (i in seq_along(legends)) {
      legends[[i]]$widths <- grid::unit(rep(legend_width / length(legends[[i]]$widths),
                                            length(legends[[i]]$widths)), "cm")
    }
    
    combined_legend <- cowplot::plot_grid(plotlist = legends, ncol = 1, align = "v")
    
    return(combined_legend)
  }
  
  cut_values <- function(data, min.cutoff = NULL, max.cutoff = NULL) {
    data_mod <- apply(data, 2, function(x) {
      values <- x; lcutoff <- -Inf; hcutoff <- Inf
      if (!is.null(min.cutoff)) {
        if (substr(min.cutoff, 1, 1) == "q") {
          qval <- as.numeric(substr(min.cutoff, 2, nchar(min.cutoff))) / 100
          lcutoff <- quantile(x, probs = qval, na.rm = TRUE)
        } else lcutoff <- as.numeric(min.cutoff)
      }
      if (!is.null(max.cutoff)) {
        if (substr(max.cutoff, 1, 1) == "q") {
          qval <- as.numeric(substr(max.cutoff, 2, nchar(max.cutoff))) / 100
          hcutoff <- quantile(x, probs = qval, na.rm = TRUE)
        } else hcutoff <- as.numeric(max.cutoff)
      }
      values <- ifelse(values < lcutoff, lcutoff, values)
      values <- ifelse(values > hcutoff, hcutoff, values)
      values
    })
    as.data.frame(data_mod)
  }
  
  ### Handle vector cutoffs
  if (length(min.cutoff) == 1) min.cutoff <- rep(min.cutoff, length(features))
  if (length(max.cutoff) == 1) max.cutoff <- rep(max.cutoff, length(features))
  
  if (is.null(cells.use)) cells.use <- rownames(object@meta.data)
  if (!is.null(colors) && length(colors) != length(features))
    stop(paste(length(features), "colors required, but", length(colors), "given."))
  
  reductionobj <- object@reductions[[tolower(reduction)]]
  cellemb <- as.data.frame(reductionobj@cell.embeddings)
  
  if (!is.null(assay)) {
    assaykey <- object@assays[[assay]]@key
    featurekey <- paste0(assaykey, features)
  } else featurekey <- features
  
  featuresexp <- FetchData(object, vars = featurekey, cells = cells.use)
  featurekey_ori <- featurekey
  featurekey <- make.names(featurekey)
  names(featurekey_ori) <- featurekey
  colnames(featuresexp) <- make.names(colnames(featuresexp))
  
  # ✅ Per-feature cutoff handling
  for (i in seq_along(featurekey)) {
    featuresexp[, i] <- cut_values(
      featuresexp[, i, drop = FALSE],
      min.cutoff = min.cutoff[min(i, length(min.cutoff))],
      max.cutoff = max.cutoff[min(i, length(max.cutoff))]
    )[, 1]
  }
  
  featuresexp$color <- exp_encolor(featuresexp, features = featurekey, pt.brightness = pt.brightness)
  
  plotdf <- merge(cellemb, featuresexp, by = "row.names")
  axis_1 <- sym(colnames(cellemb)[1]); axis_2 <- sym(colnames(cellemb)[2])
  
  base_plot <- ggplot(plotdf, aes(x = !!axis_1, y = !!axis_2)) +
    geom_point(color = plotdf$color, size = pt.size) +
    labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
    theme_dark() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = background.color),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  if (label) {
    cellemb$ident <- object@active.ident
    label_pos <- data.frame(
      x = tapply(cellemb[,1], cellemb$ident, mean),
      y = tapply(cellemb[,2], cellemb$ident, mean),
      label = names(tapply(cellemb[,2], cellemb$ident, mean))
    )
    label_layer <- if (repel)
      if (bold) geom_text_repel(data = label_pos, aes(x, y, label = label), fontface = "bold") 
    else geom_text_repel(data = label_pos, aes(x, y, label = label))
    else if (bold)
      geom_text(data = label_pos, aes(x, y, label = label), fontface = "bold")
    else geom_text(data = label_pos, aes(x, y, label = label))
    base_plot <- base_plot + label_layer
  }
  
  min.col <- hsv(1, 0, pt.brightness/100)
  fig_legends <- make_legend(featuresexp, featurekey, colors, min.col = min.col)
  plot_grid(base_plot, fig_legends, rel_widths = c(fig.ratio, 1), ncol = 2)
}



okabe_ito <- c(
  "Orange"          = "#E69F00",
  "Sky Blue"        = "#56B4E9",
  "Bluish Green"    = "#009E73",
  "Yellow"          = "#F0E442",
  "Blue"            = "#0072B2",
  "Vermillion"      = "#D55E00",
  "Reddish Purple"  = "#CC79A7",
  "Black"           = "#000000"
)

# Visualize the palette
library(scales)
show_col(okabe_ito, labels = TRUE)
