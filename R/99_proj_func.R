
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


#' Perform Fisher's exact test between two countries
#'
#' @param data A data frame containing at least `country_var` and `outcome_var`.
#' @param country_var Name of the country variable (string).
#' @param outcome_var Name of the binary outcome variable (string).
#' @param country1 Name of the first country (string).
#' @param country2 Name of the second country (string).
#'
#' @return A list with countries compared, contingency table, p-value, odds ratio,
#'   confidence interval, and method description.
fisher_compare_countries <- function(data, country_var, outcome_var, country1, country2) {
  
  country_col <- data[[country_var]]
  outcome_col <- data[[outcome_var]]
  
  subset_data <- data[country_col %in% c(country1, country2), ]
  tab <- table(subset_data[[country_var]], subset_data[[outcome_var]])
  
  fisher_res <- stats::fisher.test(tab)
  
  list(
    countries  = c(country1, country2),
    table      = tab,
    p_value    = fisher_res$p.value,
    odds_ratio = fisher_res$estimate,
    conf_int   = fisher_res$conf.int,
    method     = fisher_res$method
  )
}


#' Compute Fisher's exact test results for multiple country pairs
#'
#' @param data A data frame containing at least `country_var` and `outcome_var`.
#' @param country_var Name of the country variable (string).
#' @param outcome_var Name of the binary outcome variable (string).
#' @param pairs A list of length-2 character vectors, each giving a country pair.
#'
#' @return A tibble with one row per pair, containing countries, p-value,
#'   odds ratio, and confidence interval.
compute_fisher_results <- function(data, country_var, outcome_var, pairs) {
  
  # Run Fisher test for each pair
  results <- lapply(pairs, function(p) {
    fisher_compare_countries(
      data        = data,
      country_var = country_var,
      outcome_var = outcome_var,
      country1    = p[1],
      country2    = p[2]
    )
  })
  
  # Convert each result into a tibble row
  rows <- lapply(results, function(res) {
    tibble::tibble(
      country1   = res$countries[1],
      country2   = res$countries[2],
      p_value    = round(res$p_value, 4),
      odds_ratio = round(as.numeric(res$odds_ratio), 3),
      CI_low     = round(res$conf_int[1], 3),
      CI_high    = round(res$conf_int[2], 3)
    )
  })
  
  # Combine into a table
  dplyr::bind_rows(rows)
}

