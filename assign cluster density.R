assign_clusters_density <- function(
    data,
    pre_scaled = FALSE,
    results,
    burn_in_frac = 0.1
) {
  library(dplyr)
  library(stats)   # for density()
  
  #get data right
  if (pre_scaled == FALSE) {
    data <- scale_data(data)
  }
  
  # --- 1) Basic checks ---
  if (!"cluster_dist_df" %in% names(results)) {
    stop("results must contain 'cluster_dist_df' with columns (generation, k, count).")
  }
  if (!"cluster_centers_df" %in% names(results)) {
    stop("results must contain 'cluster_centers_df' with columns (generation, individual, cluster, dimension, value).")
  }
  
  
  
  
  cluster_dist_df    <- results$cluster_dist_df
  cluster_centers_df <- results$cluster_centers_df
  
  # We need the max generation or total generations from your run
  # If not stored, we can guess from the cluster_dist_df
  max_gen <- max(cluster_dist_df$generation, na.rm = TRUE)
  burn_in_gen <- floor(burn_in_frac * max_gen)
  
  # --- 2) Identify the "dominant" k after burn-in ---
  # Filter out generations <= burn_in_gen
  cluster_dist_post <- cluster_dist_df %>%
    filter(generation > burn_in_gen)
  
  # Summarize how many times each k appears across all relevant generations
  k_freq <- cluster_dist_post %>%
    group_by(k) %>%
    summarise(
      total_count = sum(count),  # total individuals with that k
      .groups = "drop"
    ) %>%
    arrange(desc(total_count))
  
  if (nrow(k_freq) == 0) {
    stop("No post-burn-in solutions found. Check burn_in_frac or generation data.")
  }
  
  # The top row is the "dominant k"
  k_star <- k_freq$k[1]
  # (If there's a tie, we pick the first. You can refine if needed.)
  
  message("Dominant k after burn-in is: ", k_star)
  
  # --- 3) Gather cluster centers for all individuals with k_star (post burn-in) ---
  # We also filter out generation <= burn_in_gen
  # Then we only keep rows where cluster_centers_df$cluster <= k_star
  #  (since cluster is in [1..k] for that individual)
  #  but we also must join to something that tells us the individual's k
  #  or we rely on the fact that "cluster" never exceed k for that individual.
  #
  # We'll do a left join with cluster_dist or a separate approach:
  # Usually, you'd also have "k" in cluster_centers_df if you stored it.
  # If not, we can store it in part_stats_df or something else.
  #
  # For simplicity, let's assume cluster_centers_df also had "k" or we left-join from part_stats_df.
  # We'll do a small hack: "k" might be in "part_stats_df" by (generation, individual).
  # Let's assume you have something like results$part_stats_df too.
  #
  # We'll assume cluster_centers_df has a column 'k' or we do an inner join with part_stats_df by (generation, individual).
  
  if (!"part_stats_df" %in% names(results)) {
    stop("results must contain 'part_stats_df' so we can know each individual's k.")
  }
  
  part_stats_df <- results$part_stats_df
  
  # part_stats_df: (generation, individual, k, possibly ratio, etc.)
  # cluster_centers_df: (generation, individual, cluster, dimension, value)
  # We'll combine them to filter only individuals with k=k_star
  centers_k_df <- cluster_centers_df %>%
    inner_join(part_stats_df %>% select(generation, individual, k),
               by = c("generation","individual")) %>%
    filter(generation > burn_in_gen,
           k == k_star)
  
  # centers_k_df should now have: generation, individual, k, cluster, dimension, value
  # where k == k_star, cluster in [1..k_star].
  # We'll *pool* all the center values for each cluster index, dimension.
  
  # --- 4) Build 1D kernel density for each cluster i in [1..k_star], for each dimension d ---
  # We must find out how many dims we have
  n_dims <- max(centers_k_df$dimension, na.rm = TRUE)
  
  # We'll store them in a list-of-lists: density_list[[i]][[d]] = density object
  density_list <- vector("list", k_star)
  for (i in seq_len(k_star)) {
    # For dimension-level sublists
    dens_dim_list <- vector("list", n_dims)
    for (d in seq_len(n_dims)) {
      # subset the center coords where cluster = i and dimension = d
      vals <- centers_k_df %>%
        filter(cluster == i, dimension == d) %>%
        pull(value)
      
      # If no values found for some reason, we handle that
      if (length(vals) < 2) {
        # fallback: if there's only one or zero points,
        # we artificially create a density near that point, or skip
        dens_dim_list[[d]] <- NULL
        next
      }
      
      # Build a 1D kernel density
      dens_dim_list[[d]] <- density(vals, na.rm = TRUE)
    }
    density_list[[i]] <- dens_dim_list
  }
  
  # --- 5) For each row in 'data', compute membership across k_star clusters ---
  # We'll assume 'data' is a data.frame or matrix of numeric columns (n x n_dims).
  
  # Coerce 'data' to matrix for convenience
  X <- as.matrix(data)
  if (ncol(X) != n_dims) {
    stop("Number of columns in 'data' does not match the dimension used in cluster centers.")
  }
  
  # We'll create an output tibble with one row per data point
  # columns: cluster_1_membership, ..., cluster_k*_membership, cluster_assignment
  n_rows <- nrow(X)
  membership_mat <- matrix(0, nrow = n_rows, ncol = k_star)
  
  # define a small helper to approximate the density at x from a 'density' object
  dens_approx <- function(dens_obj, x) {
    if (is.null(dens_obj)) {
      return(0)
    }
    # otherwise dens_obj is a density object (a list)
    approx_y <- approx(dens_obj$x, dens_obj$y, xout = x)$y
    if (is.na(approx_y)) approx_y <- 0
    return(approx_y)
  }
  
  
  for (row_idx in seq_len(n_rows)) {
    # The data point's coordinates
    point <- X[row_idx, ]  # length = n_dims
    
    # For each cluster i, compute product of densities across dims
    for (i in seq_len(k_star)) {
      prod_score <- 1
      for (d in seq_len(n_dims)) {
        dens_obj <- density_list[[i]][[d]]
        dens_val <- dens_approx(dens_obj, point[d])
        prod_score <- prod_score * dens_val
      }
      membership_mat[row_idx, i] <- prod_score
    }
  }
  
  # membership_mat now has unnormalized scores for each cluster
  row_sums <- rowSums(membership_mat)
  row_sums[row_sums < 1e-15] <- 1e-15  # avoid dividing by zero
  
  membership_normed <- membership_mat / row_sums
  
  # cluster assignment is the cluster with highest membership
  cluster_assignment <- max.col(membership_normed)
  
  # --- 6) Build final tibble output ---
  # e.g. cluster_1_membership, cluster_2_membership, ..., cluster_k*_membership, cluster_assignment
  out_df <- as.data.frame(membership_normed)
  colnames(out_df) <- paste0("cluster_", seq_len(k_star), "_membership")
  
  out_df$cluster_assignment <- cluster_assignment
  
  # Return
  return(out_df)
}


##########TEST##############
cluster_assignments <- assign_clusters_density(
  data = data,
  pre_scaled = FALSE,
  results = evo_results,
  burn_in_frac = .10
)


#explore
tibble(
  data,
  species = iris$Species,
  cluster = cluster_assignments$cluster_assignment
) -> cluster_w_data

cluster_w_data %>% group_by(cluster, species) %>% tally()


cluster_w_data %>% ggplot(
  aes(x = Sepal.Length, y = Sepal.Width, color = factor(cluster))) +
  geom_point() +
  facet_grid(. ~ species)


cluster_w_data %>% ggplot(
  aes(x = Sepal.Length, y = Sepal.Width, color = factor(species))) +
  geom_point() +
  facet_grid(. ~ cluster)


cluster_w_data %>% group_by(species) %>% summarise()