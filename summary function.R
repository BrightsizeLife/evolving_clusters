summarize_evo_clustering_results <- function(
    results,
    burn_in_frac = 0.1,
    print_summaries = TRUE
) {
  # 'results' is assumed to be the object returned by your main evo_clustering_wss_ratio()
  #   results$cluster_dist_df
  #   results$part_stats_df
  #   results$cluster_centers_df
  #   (and possibly other elements)
  #
  # 'burn_in_frac' is the fraction of generations to ignore as 'burn-in'.
  #   For example, if burn_in_frac=0.1 and you have 50 generations total,
  #   we ignore the first 5 generations for the summaries/plots.
  
  # Load required packages for data manipulation and plotting.
  # (If your environment already has them, you can skip the library calls or do them outside)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  
  library(dplyr)
  library(ggplot2)
  
  #############################
  # 1) Helper: Summarize Part Stats
  #############################
  summarize_part_stats <- function(part_stats_df, burn_in_frac) {
    if (!all(c("generation","ratio","r2") %in% names(part_stats_df))) {
      stop("part_stats_df must contain at least 'generation', 'ratio', and 'r2'.")
    }
    max_gen <- max(part_stats_df$generation, na.rm = TRUE)
    burn_in_gen <- floor(burn_in_frac * max_gen)
    
    # Filter out the burn-in portion
    df_filt <- part_stats_df %>%
      filter(generation > burn_in_gen)
    
    # Summarize best & average ratio, best & average R^2, etc. by generation
    summary_df <- df_filt %>%
      group_by(generation) %>%
      summarise(
        best_ratio = min(ratio, na.rm = TRUE),
        avg_ratio  = mean(ratio, na.rm = TRUE),
        best_r2    = max(r2, na.rm = TRUE),
        avg_r2     = mean(r2, na.rm = TRUE),
        .groups = "drop"
      )
    
    return(summary_df)
  }
  
  #############################
  # 2) Plot Best & Average Ratio Over Generations
  #############################
  plot_best_avg_ratio <- function(summary_df) {
    p <- ggplot(summary_df, aes(x = generation)) +
      geom_line(aes(y = best_ratio,  color = "Best ratio"), size = 1) +
      geom_line(aes(y = avg_ratio,   color = "Avg ratio"),  size = 1) +
      labs(
        x = "Generation",
        y = "WSS/TSS Ratio",
        color = "Metric"
      ) +
      ggtitle("Best vs. Average Ratio Over Generations (post burn-in)") +
      theme_minimal()
    return(p)
  }
  
  #############################
  # 3) Plot Distribution of k Over Generations
  #############################
  plot_k_distribution <- function(cluster_dist_df, burn_in_frac) {
    if (!all(c("generation","k","count") %in% names(cluster_dist_df))) {
      stop("cluster_dist_df must contain 'generation', 'k', and 'count'.")
    }
    max_gen <- max(cluster_dist_df$generation, na.rm = TRUE)
    burn_in_gen <- floor(burn_in_frac * max_gen)
    
    df_filt <- cluster_dist_df %>%
      filter(generation > burn_in_gen)
    
    # We can plot how many individuals had each k by generation
    # There are many ways to do this. We'll do a line plot.
    
    p <- ggplot(df_filt, aes(x = generation, y = count, color = factor(k))) +
      geom_line(size = 1) +
      labs(
        x = "Generation",
        y = "Count of Individuals",
        color = "Number of Clusters k"
      ) +
      ggtitle("Distribution of k across Generations (post burn-in)") +
      theme_minimal()
    
    return(p)
  }
  
  #############################
  # 4) Create Summaries & Plots
  #############################
  
  # Check that the necessary data frames exist
  if (!"part_stats_df" %in% names(results)) {
    stop("results must contain 'part_stats_df' with (generation, ratio, r2).")
  }
  if (!"cluster_dist_df" %in% names(results)) {
    stop("results must contain 'cluster_dist_df' with (generation, k, count).")
  }
  
  part_stats_df    <- results$part_stats_df
  cluster_dist_df  <- results$cluster_dist_df
  
  # a) Summarize part_stats
  part_stats_summary <- summarize_part_stats(part_stats_df, burn_in_frac)
  
  # b) Create the ratio plot
  best_avg_ratio_plot <- plot_best_avg_ratio(part_stats_summary)
  
  # c) Create the k distribution plot
  k_distribution_plot <- plot_k_distribution(cluster_dist_df, burn_in_frac)
  
  #############################
  # 5) (Optional) Print a Text Summary
  #############################
  if (print_summaries) {
    cat("=== Part Stats Summary (post burn-in) ===\n")
    print(part_stats_summary)
    
    # Show the last row (final generation or near-final) as a quick snapshot
    cat("\n--- Last Generation Stats (post burn-in) ---\n")
    last_gen <- max(part_stats_summary$generation)
    last_row <- dplyr::filter(part_stats_summary, generation == last_gen)
    print(last_row)
    cat("\n")
  }
  
  #############################
  # 6) Return Everything
  #############################
  return(list(
    part_stats_summary   = part_stats_summary,
    best_avg_ratio_plot  = best_avg_ratio_plot,
    k_distribution_plot  = k_distribution_plot
  ))
}
