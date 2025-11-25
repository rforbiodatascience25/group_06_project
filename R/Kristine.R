
# Loading libraries
library(tidyverse)
library(pheatmap)


# Remove genomes without country information:
genomes_country <- genome_metadata_aug %>%
  filter(!is.na(Country))



############# HEATMAPS ##############

# Defining function for heatmap based on species, genus, family etc...
heatmap_abundance_per_country <- function(tree_group){
  
  # Check that the column exists
  if (!tree_group %in% names(genomes_country)) {
    stop(paste("Column", tree_group, "not found in genomes_country"))
  }
  
  # Count how many times each group appears in each country
  group_freq <- genomes_country %>%
    group_by(Country, .data[[tree_group]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = .data[[tree_group]],
      values_from = count,
      values_fill = 0
    )
  
  
  # Convert to matrix, using Country as rownames
  group_freq_mat <- group_freq %>%
    column_to_rownames("Country") %>%
    as.matrix()
  
  # Remove zero-variance columns
  group_freq_mat <- group_freq_mat[, apply(group_freq_mat, 2, sd) > 0, drop = FALSE]

  # Make heatmap
  pheatmap(
    group_freq_mat,
    cluster_cols = TRUE,
    scale = "column",
    main = paste0(tree_group, " frequency per country")
  )
  
  # alternative without scaling:
  pheatmap(group_freq_mat,
           cluster_cols = TRUE,
           scale = "none",
           main = "Species frequency per country")
}

heatmap_abundance_per_country(tree_group = "Species")
heatmap_abundance_per_country(tree_group = "Genus")
heatmap_abundance_per_country(tree_group = "Family")
heatmap_abundance_per_country(tree_group = "Order")
heatmap_abundance_per_country(tree_group = "Class")
heatmap_abundance_per_country(tree_group = "Phylum")


### Normalize to the amount of samples per country first.
# Defining function for heatmap based on species, genus, family etc... (normalized by countries sample size)
heatmap_abundance_per_country_normalized <- function(tree_group){
  
  # Check that the column exists
  if (!tree_group %in% names(genomes_country)) {
    stop(paste("Column", tree_group, "not found in genomes_country"))
  }
  
  # Count how many times each group appears in each country
  norm_freq <- genomes_country %>%
    group_by(Country, .data[[tree_group]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    left_join(genomes_country %>%
                count(Country, name = "n_samples")) %>%
    mutate(freq_norm = count / n_samples) %>%
    select(Country, .data[[tree_group]], freq_norm) %>%
    pivot_wider(names_from = .data[[tree_group]],
                values_from = freq_norm,
                values_fill = 0)
  
  # Convert to matrix, using Country as rownames
  group_freq_mat <- norm_freq %>%
    column_to_rownames("Country") %>%
    as.matrix()
  
  # Remove zero-variance columns
  group_freq_mat <- group_freq_mat[, apply(group_freq_mat, 2, sd) > 0, drop = FALSE]
  
  # alternative without scaling:
  pheatmap(group_freq_mat,
           cluster_cols = TRUE,
           scale = "none",
           main = paste0(tree_group, " frequency per country (normalized)")
  )
}

heatmap_abundance_per_country_normalized(tree_group = "Species")
heatmap_abundance_per_country_normalized(tree_group = "Genus")
heatmap_abundance_per_country_normalized(tree_group = "Family")
heatmap_abundance_per_country_normalized(tree_group = "Order")
heatmap_abundance_per_country_normalized(tree_group = "Class")
heatmap_abundance_per_country_normalized(tree_group = "Phylum")



############# Multiple Correspondence Analysis ##############

MCA_abundance_per_country(tree_group = "Species")
MCA_abundance_per_country(tree_group = "Genus")
MCA_abundance_per_country(tree_group = "Family")
MCA_abundance_per_country(tree_group = "Order")
MCA_abundance_per_country(tree_group = "Class")
MCA_abundance_per_country(tree_group = "Phylum")

# 
MCA_abundance_per_country <- function(tree_group, norm = FALSE) {
  
  # Loading libraries:
  library(FactoMineR)
  library(factoextra)
  
  if(norm == TRUE){
    
    # Count how many times each group appears in each country
    group_freq <- genomes_country %>%
      group_by(Country, .data[[tree_group]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      left_join(genomes_country %>%
                  count(Country, name = "n_samples")) %>%
      mutate(freq_norm = count / n_samples) %>%
      select(Country, .data[[tree_group]], freq_norm) %>%
      pivot_wider(names_from = .data[[tree_group]],
                  values_from = freq_norm,
                  values_fill = 0)
  }
  
  else{
    # Count how many times each group appears in each country
    group_freq <- genomes_country %>%
      group_by(Country, .data[[tree_group]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(
        names_from = .data[[tree_group]],
        values_from = count,
        values_fill = 0
      )
  }
  
  # Convert numeric species frequencies into categorical bins
  group_freq_cat <- group_freq |> 
    mutate(across(-Country, ~case_when(
      . == 0 ~ "Absent",
      . == 1 ~ "Low",
      . <= 3 ~ "Medium",
      TRUE ~ "High"
    ))) |>
    mutate(across(everything(), as.factor))
  
  # Use column NAME instead of index
  mca_res <- FactoMineR::MCA(X = group_freq_cat,
                             quali.sup = "Country",
                             graph = FALSE)

  # Look at eigenvalues (explained variance)
  eig.val <- get_eigenvalue(mca_res)
  print(paste0("Eigenvalues", head(eig.val)))

  # Visualize MCA outputs of variables
  factoextra::fviz_mca_var(
    X = mca_res,
    repel = TRUE,
    labelsize = 2.5,
    title = "MCA: Species categories contributions"
  )

}


MCA_abundance_per_country(tree_group = "Phylum", norm = TRUE)




MCA_abundance_per_country_2 <- function(tree_group, norm = FALSE) {
  
  # Loading libraries:
  library(FactoMineR)
  library(factoextra)
  
  if(norm == TRUE){
    
    # Count how many times each group appears in each country
    group_freq <- genomes_country %>%
      group_by(Country, .data[[tree_group]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      left_join(genomes_country %>%
                  count(Country, name = "n_samples")) %>%
      mutate(freq_norm = count / n_samples) %>%
      select(Country, .data[[tree_group]], freq_norm) %>%
      pivot_wider(names_from = .data[[tree_group]],
                  values_from = freq_norm,
                  values_fill = 0)
  }
  
  else{
    # Count how many times each group appears in each country
    group_freq <- genomes_country %>%
      group_by(Country, .data[[tree_group]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(
        names_from = .data[[tree_group]],
        values_from = count,
        values_fill = 0
      )
  }
  
  
  print(group_freq)

  #factor it.
  
  
  # 
  # # Use column NAME instead of index
  # mca_res <- FactoMineR::MCA(X = group_freq_cat,
  #                            quali.sup = "Country",
  #                            graph = FALSE)
  # 
  # # Look at eigenvalues (explained variance)
  # eig.val <- get_eigenvalue(mca_res)
  # print(paste0("Eigenvalues", head(eig.val)))
  # 
  # # Visualize MCA outputs of variables
  # factoextra::fviz_mca_var(
  #   X = mca_res,
  #   repel = TRUE,
  #   labelsize = 2.5,
  #   title = "MCA: Species categories contributions"
  # )
  
  print(mca_res$eig)
  
  plot(mca_res)
}


MCA_abundance_per_country_2(tree_group = "Phylum", norm = TRUE)



##################### PCA ###########################
PCA_abundance_per_country <- function(tree_group) {
  
  # Install libraries
  library(broom)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)  # For non-overlapping labels
  
  # Count how many times each group appears in each country (and normalize)
  group_freq_norm <- genomes_country %>%
    group_by(Country, .data[[tree_group]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    left_join(genomes_country %>%
                count(Country, name = "n_samples")) %>%
    mutate(freq_norm = count / n_samples) %>%
    select(Country, .data[[tree_group]], freq_norm) %>%
    pivot_wider(names_from = .data[[tree_group]],
                values_from = freq_norm,
                values_fill = 0)

  print("Normalized frequency table:")
  print(group_freq_norm)
  
  # Extract country names
  country_names <- group_freq_norm$Country
  
  # Remove Country column and run PCA on numeric data only
  pca_fit <- group_freq_norm %>% 
    select(-Country) %>%  # Remove the Country column
    scale() %>%            # scale data
    prcomp()               # do PCA
  
  # Get variance explained
  variance_explained <- pca_fit %>%
    tidy(matrix = "eigenvalues")
  
  print("Variance explained:")
  print(variance_explained)
  
  # Get rotation matrix (loadings)
  rotation_matrix <- pca_fit %>%
    tidy(matrix = "rotation")
  
  print("Rotation matrix (loadings):")
  print(rotation_matrix)
  
  # Get PCA scores and add country names back
  pca_scores <- pca_fit %>%
    augment(group_freq_norm %>% select(-Country)) %>%
    mutate(Country = country_names)
  
  # Calculate % variance for axis labels
  pc1_var <- round(variance_explained$percent[1] * 100, 1)
  pc2_var <- round(variance_explained$percent[2] * 100, 1)
  
  # Prepare loadings for biplot (scale them for visualization)
  loadings_biplot <- rotation_matrix %>%
    filter(PC %in% c(1, 2)) %>%
    pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
    mutate(
      PC1 = PC1 * 3,  # Scale factor for visibility (adjust as needed)
      PC2 = PC2 * 3
    )
  
  # Plot 1: PCA biplot with countries AND species loadings
  p1 <- ggplot() +
    # Add country points and labels
    geom_point(data = pca_scores, 
               aes(x = .fittedPC1, y = .fittedPC2, color = Country), 
               size = 4,
               show.legend = TRUE) +
    geom_text(data = pca_scores, 
              aes(x = .fittedPC1, y = .fittedPC2, label = Country, color = Country),
              vjust = -1,
              size = 4,
              show.legend = FALSE) +
    # Add loading arrows
    geom_segment(data = loadings_biplot,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 color = "gray50",
                 linewidth = 0.7) +
    geom_text_repel(data = loadings_biplot,
                    aes(x = PC1, y = PC2, label = column),
                    color = "gray10",
                    size = 3.5,
                    fontface = "italic",
                    box.padding = 0.5,      # Space around labels
                    point.padding = 0.3,    # Space from arrow tips
                    segment.color = "gray60",  # Connection lines
                    segment.size = 0.3,
                    max.overlaps = Inf) +   # Allow all labels
    labs(title = paste0("PCA Biplot: Countries + ", tree_group),
         x = paste0("PC1 (", pc1_var, "%)"),
         y = paste0("PC2 (", pc2_var, "%)")) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
  
  
  # Plot 2: Scree plot (variance explained)
  p2 <- ggplot(variance_explained, aes(x = PC, y = percent)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Variance Explained by Each PC",
         x = "Principal Component",
         y = "Variance Explained (%)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  # Saving plots
  ggsave(p1, "../results/VarianceExplained.png")
  ggsave(p2, "../results/VarianceExplained.png")

}


# Run the function
PCA_abundance_per_country(tree_group = "Phylum")

