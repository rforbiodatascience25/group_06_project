
### Analysis 4.

#' Normalize data by sample size within each country group
#'
#' @param data A data frame to normalize.
#' @param tree_group String indicating the lineage level.
#' @return A wide data frame of normalized frequencies, one row per country.
normalize <- function(data, tree_group) {
  
  # Summarize counts per country × tree_group
  group_freq <- data %>%
    dplyr::group_by(Country, .data[[tree_group]]) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(
      data %>% dplyr::count(Country, name = "n_samples"),
      by = "Country"
    ) %>%
    dplyr::mutate(freq_norm = count / n_samples) %>%
    dplyr::select(Country, .data[[tree_group]], freq_norm) %>%
    tidyr::pivot_wider(
      names_from  = .data[[tree_group]],
      values_from = freq_norm,
      values_fill = 0
    )
  
  # Aggregate duplicates to ensure one row per Country
  group_freq <- group_freq %>%
    dplyr::group_by(Country) %>%
    dplyr::summarise(dplyr::across(everything(), mean), .groups = "drop")
  
  return(group_freq)
}


#' Create a heatmap of abundancy in each country based on the selected tree_group
#'
#' @param tree_group String indicating the lineage level.
#' @return A heatmap of the abundancy for the selected tree_group per country.
heatmap_abundance_per_country <- function(tree_group){
  
  if (!tree_group %in% names(genomes_country)) {
    stop(paste("Column", tree_group, "not found in genomes_country"))
  }
  
  group_freq <- normalize(data = genomes_country,
                          tree_group = tree_group)
  
  group_freq_mat <- group_freq |>
    column_to_rownames("Country") |> 
    as.matrix()
  
  group_freq_mat <- group_freq_mat[, apply(group_freq_mat, 2, sd) > 0, drop = FALSE]
  
  pheatmap::pheatmap(
    mat = group_freq_mat,
    cluster_cols = TRUE,
    scale = "none",
    main = paste0(tree_group, " frequency per country (normalized)")
  )
}

#' Counting the species (or another tree_group) per country.
#'
#' @param data A data frame to normalize.
#' @param tree_group String indicating the lineage level.
#' @return A wide data frame of frequencies, one row per country.
counting_species <- function(data, tree_group){
  data |> 
    group_by(Country, .data[[tree_group]]) |>
    summarise(count = n(), .groups = "drop") |> 
    pivot_wider(names_from = .data[[tree_group]],
                values_from = count,
                values_fill = 0)
}


