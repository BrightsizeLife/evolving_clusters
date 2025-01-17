library(parallel)
library(doParallel)

run_chains_in_parallel <- function(
    data, 
    n_chains = 4,
    ...
) {
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  
  results_list <- foreach(chain_idx = 1:n_chains) %dopar% {
    # each chain in parallel
    set.seed(100 + chain_idx)
    chain_res <- evo_clustering_wss_ratio(data = data, ...)
    
    chain_res$cluster_dist_df$chain    <- chain_idx
    chain_res$cluster_centers_df$chain <- chain_idx
    chain_res$part_stats_df$chain      <- chain_idx
    
    chain_res
  }
  
  stopCluster(cl)
  
  # combine
  final_populations <- lapply(results_list, `[[`, "final_population")
  cluster_dist_df   <- do.call(rbind, lapply(results_list, `[[`, "cluster_dist_df"))
  cluster_centers_df<- do.call(rbind, lapply(results_list, `[[`, "cluster_centers_df"))
  part_stats_df     <- do.call(rbind, lapply(results_list, `[[`, "part_stats_df"))
  
  list(
    final_populations = final_populations,
    cluster_dist_df   = cluster_dist_df,
    cluster_centers_df= cluster_centers_df,
    part_stats_df     = part_stats_df
  )
}
