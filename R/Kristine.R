
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

#dependencies:
# install.packages(c("tidyverse", "FactoMineR", "factoextra"))

MCA_abundance_per_country <- function(tree_group) {
  
  # Loading libraries:
  library(FactoMineR)
  library(factoextra)
  
  # Count how many times each group appears in each country
  group_freq <- genomes_country %>%
    group_by(Country, .data[[tree_group]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = .data[[tree_group]],
      values_from = count,
      values_fill = 0
    )
  
  # Convert numeric species frequencies into categorical bins
  group_freq_cat <- group_freq |> 
    mutate(across(-Country, ~case_when(
      . == 0 ~ "Absent",
      . == 1 ~ "Low",
      . <= 3 ~ "Medium",
      TRUE ~ "High"
    ))) |>
    mutate(across(everything(), as.factor))
  
  # The fix → use column NAME instead of index
  mca_res <- FactoMineR::MCA(X = group_freq_cat,
                             quali.sup = "Country",
                             graph = FALSE)
  
  # Using factorextra package to extract the results for visualization.

  # # Visualize MCA outputs of individuals
  # factoextra::fviz_mca_ind(
  #   X = mca_res,
  #   habillage = group_freq_cat$Country,
  #   addEllipses = TRUE,
  #   repel = TRUE,
  #   title = "MCA: Countries based on bacteria presence categories"
  # )
  # 
  # # Visualize MCA outputs of variables
  # factoextra::fviz_mca_var(
  #   X = mca_res,
  #   repel = TRUE,
  #   title = "MCA: Species categories contributions"
  # )
  
  # Visualize MCA biplot of individuals and variables
  factoextra::fviz_mca_biplot(
    X = mca_res,
    repel = TRUE,
    habillage = group_freq_cat$Country,
    palette = "jco",
    # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    title = "MCA Biplot: Countries + Species"
  )
}

MCA_abundance_per_country(tree_group = "Species")
MCA_abundance_per_country(tree_group = "Genus")
MCA_abundance_per_country(tree_group = "Family")
MCA_abundance_per_country(tree_group = "Order")
MCA_abundance_per_country(tree_group = "Class")
MCA_abundance_per_country(tree_group = "Phylum")


# second version - use some of this in the original above. 
MCA_abundance_per_country <- function(tree_group) {
  
  # Loading libraries:
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  
  # Count how many times each group appears in each country
  group_freq <- genomes_country %>%
    group_by(Country, .data[[tree_group]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = .data[[tree_group]],
      values_from = count,
      values_fill = 0
    )
  
  # Convert numeric species frequencies into categorical bins
  group_freq_cat <- group_freq |> 
    mutate(across(-Country, ~case_when(
      . == 0 ~ "Absent",
      . == 1 ~ "Low",
      . <= 3 ~ "Medium",
      TRUE ~ "High"
    ))) |>
    mutate(across(everything(), as.factor))
  
  # Run MCA
  mca_res <- FactoMineR::MCA(X = group_freq_cat,
                             quali.sup = "Country",
                             graph = FALSE)
  
  # Extract individual (country) coordinates
  ind_coord <- as.data.frame(mca_res$ind$coord[, 1:2])
  colnames(ind_coord) <- c("Dim1", "Dim2")
  ind_coord$Country <- group_freq_cat$Country
  
  # Extract variable coordinates
  var_coord <- as.data.frame(mca_res$var$coord[, 1:2])
  colnames(var_coord) <- c("Dim1", "Dim2")
  var_coord$label <- rownames(var_coord)
  
  # Extract the category (Absent, Low, Medium, High) from variable names
  var_coord$category <- sapply(strsplit(var_coord$label, "_"), function(x) x[length(x)])
  
  # Define colors for each category
  category_colors <- c("Absent" = "gray70", 
                       "Low" = "#FED976", 
                       "Medium" = "#FD8D3C", 
                       "High" = "#E31A1C")
  
  var_coord$color <- category_colors[var_coord$category]
  
  # Get eigenvalues for axis labels
  eig <- get_eigenvalue(mca_res)
  
  # Build plot from scratch
  p <- ggplot() +
    # Add country points
    geom_point(data = ind_coord, 
               aes(x = Dim1, y = Dim2, color = Country),
               size = 3) +
    geom_text(data = ind_coord,
              aes(x = Dim1, y = Dim2, label = Country, color = Country),
              vjust = -1, size = 3.5) +
    # Add variable points colored by abundance
    geom_point(data = var_coord, 
               aes(x = Dim1, y = Dim2),
               color = var_coord$color,
               shape = 17,  # Triangle
               size = 2) +
    geom_text(data = var_coord, 
              aes(x = Dim1, y = Dim2, label = label),
              color = var_coord$color,
              size = 2.5, 
              vjust = -0.7) +
    labs(title = "MCA Biplot: Countries + Species",
         x = paste0("Dim 1 (", round(eig[1, 2], 1), "%)"),
         y = paste0("Dim 2 (", round(eig[2, 2], 1), "%)")) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}



############# OTHER STUFF ###############

# #### FROM CHAT - the rest of the analysis #### - OBS! Not changed to be for higher up the tree.
# # 6. Species presence/absence (just 0/1)
# pa_mat <- df_clean %>%
#   mutate(value = 1) %>%
#   group_by(Country, Species) %>%
#   summarise(value = 1, .groups = "drop") %>%
#   pivot_wider(
#     names_from = Species,
#     values_from = value,
#     values_fill = 0
#   ) %>%
#   column_to_rownames("Country")
# 
# # Keep only species that vary
# pa_mat <- pa_mat[, apply(pa_mat, 2, sd) > 0]
# 
# # 7. Heatmap (presence/absence only)
# pheatmap(
#   pa_mat,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   main = "Presence/absence of species by country"
# )
# 
# # 8. Statistical test (chi-square/Fisher)
# # Which species differ between countries?
# species_list <- colnames(pa_mat)
# test_results <- lapply(species_list, function(sp) {
#   tbl <- table(df_clean$Country, df_clean$Species == sp)
#   if(any(tbl < 5)) {
#     p <- fisher.test(tbl)$p.value
#   } else {
#     p <- chisq.test(tbl)$p.value
#   }
#   data.frame(Species = sp, p_value = p)
# })
# 
# species_stats <- bind_rows(test_results) %>%
#   arrange(p_value)
# 
# # View species most associated with a country
# head(species_stats, 20)
# 
# # Optional 9. PCoA using genome-level presence/absence
# # (each genome is a sample)
# genome_pa_mat <- df_clean %>%
#   mutate(value = 1) %>%
#   group_by(Genome, Country, Species) %>%
#   summarise(value = 1, .groups = "drop") %>%
#   pivot_wider(
#     names_from = Species,
#     values_from = value,
#     values_fill = 0
#   )
# 
# meta <- genome_pa_mat[, c("Genome", "Country")]
# genome_pa <- genome_pa_mat %>%
#   select(-Genome, -Country)
# 
# # Jaccard distance
# dist_jacc <- vegdist(genome_pa, method = "jaccard", binary = TRUE)
# 
# # PCoA
# pcoa <- cmdscale(dist_jacc, eig = TRUE, k = 2)
# 
# pcoa_df <- data.frame(
#   PC1 = pcoa$points[,1],
#   PC2 = pcoa$points[,2],
#   Country = meta$Country
# )
# 
# # Plot
# library(ggplot2)
# ggplot(pcoa_df, aes(PC1, PC2, color = Country)) +
#   geom_point(size = 3, alpha = 0.7) +
#   theme_minimal() +
#   labs(title = "PCoA of genomes (Jaccard distance)")