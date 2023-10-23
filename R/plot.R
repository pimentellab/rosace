#' @import dplyr
#' @import ggplot2
#' @import colorspace
#' @import ggpubr
#' @import utils
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate Score Heatmap
#'
#' `scoreHeatmap` generates a series of heatmaps for visualizing the
#'  score patterns across different positions.
#'
#' @param data Scores data frame. Expected to have columns containing information about 
#' position, control amino acid, mutated amino acid, mutation type, and score. 
#' @param savedir Character string specifying the directory to save plots.
#' @param ctrl.name String specifying the control mutation name. Default is `synonymous`.
#' @param pos.col Column name in `data` for mutation positions. Default is `position`.
#' @param wt.col Column name in `data` for wildtype amino acids. Default is `wildtype`.
#' @param mut.col Column name in `data` for mutated amino acids. Default is `mutation`.
#' @param type.col Column name in `data` for mutation types. Default is `type`.
#' @param score.col Column name in `data` for mutation scores. Default is `mean`. 
#' @param aa.order Character vector defining the order of amino acid mutations in the y-axis.
#' Default to using mutations from `data` in alphabetical order.
#' @param npos Integer specifying the number of positions per subplot. Default is `100`.
#' @param ncol Integer specifying the number of columns of subplots. Default is `1`.
#' @param pos.step Integer specifying the steps between x-axis labels. Default is `5`.
#' @param x.text Numeric value for x-axis text size. Default is `6`.
#' @param y.text Numeric value for y-axis text size. Default is `3`.
#' @param seq.text Numeric value for wildtype sequence text size. Default is `1.1`
#' @param c.pallete Character string or vector defining the color palette. Default is `'RdBu'`.
#' @param c.scale List of parameters definning the score color scale.
#' @param ht Numeric value for the height of the saved plot. Default is `11`.
#' @param wd Numeric value for the width of the saved plot. Default is `8.5`.
#' @param name Character string specifying the base name of the saved file.
#' @param savepdf Logical indicating whether to also save a PDF version of the plot. 
#' Default is `TRUE`.
#' @param show Logical indicating whether or not to display the plot in the viewer. 
#' Default is `FALSE`.
#'
#' @return NULL.
#'
#' @examples
#' \dontrun{
#' scoreHeatmap(scores.data, savedir = "./plots/", name = "Heatmap", savepdf = TRUE)
#' }
#'
#' @export
#'
scoreHeatmap <- function(data,
                         savedir,
                         ctrl.name = "synonymous",
                         pos.col = "position",
                         wt.col = "wildtype",
                         mut.col = "mutation",
                         type.col = "type",
                         score.col = "mean",
                         aa.order = NA,
                         npos = 100,
                         ncol = 1,
                         pos.step = 5,
                         x.text = 6,
                         y.text = 3,
                         seq.text = 1.1,
                         c.pallete = 'RdBu',
                         c.scale = list(),
                         ht = 11,
                         wd = 8.5,
                         name = "Heatmap",
                         savepdf = TRUE,
                         show = FALSE
){
  
  # determine certain plot properties
  if (is.na(aa.order)) {
    aa.order = unique(data[[mut.col]])
  }
  pos.order <- levels(factor(data[[pos.col]]))
  npanel <- length(pos.order) %/% npos +1
  nrow <- npanel / ncol
  c.default <- list(palette = c.pallete, mid = 0, rev=TRUE, na.value = '#E5E5E5')
  c.args <- modifyList(c.default, c.scale)
  
  # parse the positions
  starts <- seq(1, length(pos.order), by = npos)
  if (length(starts) > 1) {
    ends <- c((starts - 1)[2:length(starts)], length(pos.order))
  } else {
    ends <- length(pos.order)
  }
  
  plot_list <- lapply(1:length(starts), function(i) {
    
    start <- starts[i]
    end <- ends[i]
    legend = ifelse(i == length(starts), 'right', 'none')
    
    # subset data for positions in range
    sub_data <- data[data[[pos.col]] %in% pos.order[start:end], ]
    sub_pos.order <- levels(factor(sub_data[[pos.col]]))
    
    # create subplot heatmaps
    p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
      y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
      geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
      do.call(scale_fill_continuous_divergingx, c.args) + 
      scale_color_manual(values = c(NA,'green'), name = ctrl.name) +
      scale_x_discrete(breaks =  sub_pos.order[seq(1, length(sub_pos.order), by = pos.step)]) +
      coord_fixed(ratio = 1, clip = "off") +
      theme(
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
        axis.title.y = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = x.text),
        axis.text.y = element_text(size = y.text),
        axis.ticks = element_blank(),
        legend.position= legend,
        legend.box = "horizontal"
      ) +
      geom_text(data = sub_data, 
                aes(x = factor(.data[[pos.col]]), label = factor(.data[[wt.col]]), y = Inf), 
                vjust = -1, check_overlap = TRUE, size = seq.text) +
      labs(y = "Mutation", x = "Position")
    
    
    return(p)
  })
  
  p_all <- do.call(ggarrange, c(plot_list, list(nrow = nrow, ncol = ncol)))
  
  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  ggsave(file.path(savedir, paste0(name,".png")), plot = p_all, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p_all, height = ht, width = wd)
  } 
  
  if (show) {
    cat(paste("Showing the first", as.character(npos) , "positions.",
    "Full figure can be found in the saved directory."))
    print(plot_list[[1]])
  }
}

#' Generate Score Violin Plot
#'
#'`scoreVlnplot` generates a series of violin plots for visualizing the 
#' distribution of scores across different positions.
#'
#' @param data Scores data frame. Expected to have columns containing information about 
#' position, control amino acid, mutated amino acid, mutation type, and score.
#' @param savedir Character string specifying the directory to save plots.
#' @param pos.col Column name in `data` for mutation positions. Default is `position`.
#' @param wt.col Column name in `data` for wildtype amino acids. Default is `wildtype`.
#' @param score.col Column name in `data` for mutation scores. Default is `mean`.
#' @param jitter Logical indicating whether to add jitter points to violin plot. 
#' Default is `FALSE`.
#' @param c.fill String indicating the fill color for the violin plots.
#' @param npos Integer specifying the number of positions per subplot. Default is `50`.
#' @param ncol Integer specifying the number of columns of subplots. Default is `1`.
#' @param pos.step Integer specifying the steps between x-axis labels. Default is `5`.
#' @param x.text Numeric value for x-axis text size. Default is `8`.
#' @param y.text Numeric value for x-axis text size. Default is `8`.
#' @param seq.text Numeric value for wildtype sequence text size. Default is `3`
#' @param pt.size Numeric value specifying the size of jitter points. Default is `0.2`
#' @param ht Numeric value for the height of the saved plot. Default is `30`.
#' @param wd Numeric value for the width of the saved plot. Default is `15`.
#' @param savedir A string providing the directory path to save plots.
#' @param name Character string specifying the base name of the saved file.
#' @param savepdf Logical indicating whether to also save a PDF version of the plot. 
#' Default is `TRUE`.
#' @param show Logical indicating whether or not to display the plot in the viewer. 
#' Default is `FALSE`.
#'
#' @return NULL.
#'
#' @examples
#' \dontrun{
#' scoreVlnplot(data = scores.data, jitter = TRUE, c.fill = "lightpink",
#'             savedir = "./plots/", name = "ViolinPlot")
#' }
#' 
#' @export
#'
scoreVlnplot <- function(data,
                         savedir,
                         pos.col = "position",
                         wt.col = "wildtype",
                         score.col = "mean",
                         jitter = FALSE,
                         c.fill = "lightblue",
                         npos = 50,
                         ncol = 1,
                         pos.step = 1,
                         x.text = 8,
                         y.text = 8,
                         seq.text = 3,
                         pt.size = 0.2,
                         ht = 30,
                         wd = 15,
                         name = "Violinplot",
                         savepdf = TRUE,
                         show = FALSE
){
  
  # determine certain plot properties
  pos.order <- levels(factor(data[[pos.col]]))
  npanel <- length(pos.order) %/% npos +1
  nrow <- npanel / ncol
  
  # parse the positions
  starts <- seq(1, length(pos.order), by = npos)
  if (length(starts) > 1) {
    ends <- c((starts - 1)[2:length(starts)], length(pos.order))
  } else {
    ends <- length(pos.order)
  }
  
  plot_list <- lapply(1:length(starts), function(i) {
    
    start <- starts[i]
    end <- ends[i]
    legend = ifelse(i == length(starts), 'right', 'none')
    
    # subset data for positions in range
    sub_data <- data[data[[pos.col]] %in% pos.order[start:end], ]
    sub_pos.order <- levels(factor(sub_data[[pos.col]]))
    
    # create subplot violinplots
    p <- ggplot(sub_data, aes(x = factor(.data[[pos.col]]), y = .data[[score.col]])) + 
      geom_violin(scale = "count", fill = c.fill, colour = c.fill) +
      scale_x_discrete(breaks =  sub_pos.order[seq(1, length(sub_pos.order), by = pos.step)]) +
      coord_fixed(ratio = 1, clip = "off") +
      theme(
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
        axis.title.y = element_text(margin = margin(t = 3, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = x.text),
        axis.text.y = element_text(size = y.text),
        axis.ticks = element_blank()
      ) +
      geom_text(data = sub_data, aes(x = factor(.data[[pos.col]]), 
                label = factor(.data[[wt.col]]), y = Inf), 
                vjust = -1, check_overlap = TRUE, size = seq.text) +
      labs(y = "Score", x = "Position")
    
    if (jitter) {
      p <- p+geom_jitter(height = 0, width = 0.1, size = pt.size) 
    }
    
    return(p)
  })
  
  p_all <- do.call(ggarrange, c(plot_list, list(nrow = nrow, ncol = ncol)))
  
  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  ggsave(file.path(savedir, paste0(name,".png")), plot = p_all, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p_all, height = ht, width = wd)
  } 
  
  if (show) {
    cat(paste("Showing the first", as.character(npos) , "positions.", 
    "Full figure can be found in the saved directory."))
    print(plot_list[[1]])
  }
}

#' Generate Score Density Plot
#'
#'`scoreDensity` generates a density plot for visualizing the distribution of scores 
#' across different mutation types.
#'
#' @param data Scores data frame. Expected to have columns containing information about 
#' position, control amino acid, mutated amino acid, mutation type, and score. 
#' @param savedir Character string specifying the directory to save plots.
#' @param type.col Column name in `data` for mutation types. Default is `type`.
#' @param score.col Column name in `data` for mutation scores. Default is `mean`.  
#' @param hist Logical indicating whether to plot count histogram or density. 
#' Default is `FALSE`.
#' @param nbins Numeric value specifying the number of bins for the histogram. 
#' Default is `30`.
#' @param c.fill Vector indicating the fill color for mutation types.
#' @param alpha Numeric value between 0-1 indicating the fill transparency. 
#' Default is `0.5`
#' @param x.text Numeric value for x-axis text size. Default is `10`.
#' @param y.text Numeric value for x-axis text size. Default is `10`.
#' @param scale.free Logical indicating whether to make score range proportional.
#' Default is `FALSE`.
#' @param space.free Logical indicating whether to make panel heights variable. 
#' Default is `FALSE`.
#' @param ht Numeric value for the height of the saved plot. Default is `10`.
#' @param wd Numeric value for the width of the saved plot. Default is `8`.
#' @param name Character string specifying the base name of the saved file.
#' @param savepdf Logical indicating whether to also save a PDF version of the plot. 
#' Default is `TRUE`.
#' @param show Logical indicating whether or not to display the plot in the viewer. 
#' Default is `TRUE`.
#'
#' @return NULL.
#'
#' @examples
#' \dontrun{
#' scoreDensity(data = scores.data, scale.free = TRUE,
#'             savedir = "./plots/", name = "DensityPlot")
#' }
#' 
#' @export
#'
scoreDensity <- function(data,
                         savedir,
                         type.col = "type",
                         score.col = "mean",
                         hist = FALSE,
                         nbins = 30,
                         c.fill = c('#FF7575', 'lightgreen', "#7298BF"),
                         alpha = 0.5,
                         x.text = 10,
                         y.text = 10,
                         scale.free = FALSE,
                         space.free = FALSE,
                         ht = 10,
                         wd = 8,
                         name = "DensityPlot",
                         savepdf = TRUE,
                         show = TRUE
){
  
  sc <- ifelse(scale.free, "free_y", "fixed")
  sp <- ifelse(space.free, "free_y", "fixed")
  if (length(c.fill) != length(levels(factor(data[[type.col]])))) {
    warning("Length of color vector does not match the number of mutation types.")
  }
  
  p <- ggplot(data, aes(x = .data[[score.col]], fill = .data[[type.col]])) +
    facet_grid(.data[[type.col]] ~ ., scales = sc, space = sp) +
    theme_minimal() +
    scale_fill_manual(values = c.fill) +
    theme(
      strip.text.y = element_blank(),  
      strip.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(size = x.text),
      axis.text.y = element_text(size = y.text)
    )
  
  if (hist) {
    p <- p + geom_histogram(aes(alpha = 0.5), bins = nbins, position="dodge") + 
      labs(x = "Score", y = "Count") +
      scale_alpha(guide = "none")
  }
  else{
    p <- p +geom_density(alpha = alpha) + labs(x = "Score", y = "Density")
  }
  
  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  ggsave(file.path(savedir, paste0(name,".png")), plot = p, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p, height = ht, width = wd)
  } 
  
  if (show) {
    print(p)
  }
}