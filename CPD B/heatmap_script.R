pacman::p_load(data.table,
               magrittr,
               tidyverse,
               pheatmap,
               RColorBrewer,
               gridExtra)

#raw_data <- fread("Cpd B 3ad vs 3 HC_normalizedCountsWithAnnotations.txt")
#  write.csv(raw_data, "raw_data.csv")
raw_data <- fread("raw_data.csv")
  data_matrix <- raw_data %>% 
    select(-c(GeneID, Alias, Description))

#==========================================
#==========================================
AD <- c(1, 3, 5)
HC <- c(7, 9, 11)
AD_treated <- c(2, 4, 6)
HC_treated <- c(8, 10, 12)

#Total ===================================
generate.heatmap <- function(data_matrix){
  DF <- column_to_rownames(data_matrix, "Symbol")
  
  metadata <- data.frame(1:12) %>%
    rename(Trial = X1.12) %>%
    mutate(Group = ifelse(Trial %in% AD, "AD", 
                              ifelse(Trial %in% HC, "HC",
                                     ifelse(Trial %in% AD_treated, "AD Treated", "HC Treated")))) %>% 
    column_to_rownames(., "Trial")
  
  annotation_colors <- list(Group = c("#a6cee3", "#1f78b4", "#FB9A99", "#E31A1C"))
    names(annotation_colors$Group) <- unique(metadata$Group)
  
  heatmap <- pheatmap(DF,
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      cluster_cols = TRUE, cluster_rows = TRUE,
                      scale = "row",
                      cex = 1,
                      color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      border_color = FALSE,
                      annotation_col = metadata,
                      annotation_colors = annotation_colors)
  
  return(heatmap)
}
#Pheatmap Save======================================
save.pheatmap <- function(x, filename, width=1080, height=2160) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#Subset ===============================================
generate.grouped.heatmap <- function(data_matrix, group1, group2){
  trials = append(group1, group2)
  
  DF <- column_to_rownames(data_matrix, "Symbol")
  DF <- DF[, trials]
  
  metadata <- data.frame(Trial = 1:12) %>%
    mutate(Group = ifelse(Trial %in% AD, "AD", 
                          ifelse(Trial %in% HC, "HC",
                                 ifelse(Trial %in% AD_treated, "AD Treated", "HC Treated")))) %>%
    filter(Trial %in% colnames(DF)) %>% 
    column_to_rownames(., "Trial")
  
  annotation_colors <- list(Group = c("#a6cee3", "#1f78b4", "#FB9A99", "#E31A1C"))
  names(annotation_colors$Group) <- c("AD Treated", "AD", "HC Treated", "HC")
  
  heatmap <- pheatmap(DF,
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      cluster_cols = TRUE, cluster_rows = TRUE,
                      scale = "row",
                      cex = 1,
                      color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
                      clustering_distance_rows = "euclidean",
                      clustering_method = "complete",
                      border_color = FALSE,
                      annotation_col = metadata,
                      annotation_colors = annotation_colors)
  
  return(heatmap)
}

temp <- generate.grouped.heatmap(data_matrix, HC, AD)
  save.pheatmap(temp, "HCxAD.png")
temp <- generate.grouped.heatmap(data_matrix, HC, HC_treated)
  save.pheatmap(temp, "HCxHC_Treated.png")
temp <- generate.grouped.heatmap(data_matrix, AD, AD_treated)
  save.pheatmap(temp, "ADxAD_Treated.png")
temp <- generate.grouped.heatmap(data_matrix, HC, AD_treated)
  save.pheatmap(temp, "HCxAD_Treated.png")

