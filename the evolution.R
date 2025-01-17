
data <- iris[-5]

##########################################################
# scale_data
# - Purpose: Optionally scale each column of data 
#   to the range [0,1].
# - If scale=FALSE, returns data unchanged.
##########################################################
scale_data <- function(data, scale = TRUE) {
  if (!scale) {
    return(data)
  }
  
  # For each column, transform to [0,1]
  min_vals <- apply(data, 2, min)
  max_vals <- apply(data, 2, max)
  
  # Avoid zero-division in case min==max
  scaled <- mapply(
    function(col, mn, mx) {
      rng <- mx - mn
      if (rng < 1e-15) rng <- 1e-15
      (col - mn) / rng
    },
    data, min_vals, max_vals
  )
  
  # Convert back to data frame
  scaled_data <- as.data.frame(scaled)
  return(scaled_data)
}

scaled_data <- scale_data(data)



##########################################################
# initialize_population
# - Purpose: Create the initial generation of individuals.
#   Each individual has:
#    - k : the number of clusters (sampled between min_centers and max_centers)
#    - dna: a numeric vector storing all cluster centers 
#           (k centers x n_dims per center).
#
#   For example, if k=3 and n_dims=2, we have 3*2=6 DNA values.
##########################################################
initialize_population <- function(
    pop_size     = 100,   # number of individuals in the population
    min_centers  = 2,     # min number of clusters per individual
    max_centers  = 5,     # max number of clusters per individual
    n_dims       = 2      # number of dimensions in the data
) {
  population <- vector("list", length = pop_size)
  
  for (i in seq_len(pop_size)) {
    # Randomly choose how many clusters this individual has
    k <- sample(min_centers:max_centers, 1)
    
    # For each center, we need n_dims values
    # So total DNA length = k * n_dims.
    # We'll initialize them uniformly in [0,1].
    dna <- runif(k * n_dims, min = 0, max = 1)
    
    population[[i]] <- list(
      k   = k,   # how many clusters
      dna = dna  # the cluster centers, flattened
    )
  }
  
  return(population)
}


initialize_population <- function(
    pop_size     = 100,   # number of individuals in the population
    min_centers  = 2,     # min number of clusters per individual
    max_centers  = 5,     # max number of clusters per individual
    n_dims       = 2      # number of dimensions in the data
) {
  population <- vector("list", length = pop_size)
  
  for (i in seq_len(pop_size)) {
    k <- sample(min_centers:max_centers, 1)
    dna <- runif(k * n_dims, min = 0, max = 1)
    
    pop_id <- paste0("init_", i)  # or something unique
    population[[i]] <- list(
      k       = k,
      dna     = dna,
      parents = NULL,   # no parents => they are the “Adam/Eve” generation
      id      = pop_id
    )
  }
  return(population)
}



individuals = initialize_population(n_dims = ncol(scaled_data))



##########################################################
# calc_wss_tss_ratio
# - Purpose: Given the full dataset and an individual's 
#   cluster centers, compute the WSS/TSS ratio.
#
#   Steps:
#     1) Reconstruct cluster centers from DNA (k x n_dims).
#     2) Assign each data point to the nearest center.
#     3) Compute TSS (sum of squared distances to overall mean).
#     4) Compute WSS (sum of squared distances from points to 
#        their cluster means).
#     5) ratio = WSS / TSS (0 <= ratio <= 1).
##########################################################
calc_wss_tss_ratio <- function(data, individual) {
  # Data dimensions
  n_dims <- ncol(data)
  n_rows <- nrow(data)
  
  # Extract cluster centers from individual's DNA
  k <- individual$k
  centers <- matrix(individual$dna, nrow = k, ncol = n_dims, byrow = TRUE)
  
  # 1) Assign each row in data to its nearest center
  #    We'll compute squared distances to each center:
  dist_mat <- sapply(seq_len(k), function(ci) {
    rowSums( (data - matrix(centers[ci, ],
                            nrow = n_rows,
                            ncol = n_dims,
                            byrow = TRUE))^2 )
  })
  
  # For each data point, find index of cluster with the minimum distance
  # max.col() trick: negative distances => min
  cluster_assign <- max.col(-dist_mat)
  
  # 2) Compute TSS (total sum of squares) relative to the overall mean
  overall_mean <- colMeans(data)
  TSS <- sum(rowSums( (data - matrix(overall_mean,
                                     nrow = n_rows,
                                     ncol = n_dims,
                                     byrow = TRUE))^2 ))
  
  # 3) Compute WSS by summing squared distances to cluster centroids
  WSS <- 0
  for (c_idx in seq_len(k)) {
    # subset of data assigned to cluster c_idx
    cluster_points <- (cluster_assign == c_idx)
    if (!any(cluster_points)) next  # if no points assigned, skip
    
    sub_data <- data[cluster_points, , drop = FALSE]
    sub_rows <- nrow(sub_data)
    center_c <- centers[c_idx, ]
    
    # sum of squared distances
    WSS_c <- sum(rowSums( (sub_data - matrix(center_c,
                                             nrow = sub_rows,
                                             ncol = n_dims,
                                             byrow = TRUE))^2 ))
    WSS <- WSS + WSS_c
  }
  
  # 4) ratio = WSS / TSS
  #    If TSS is extremely small (degenerate case), set ratio ~ 1 or some large number
  if (TSS < 1e-15) {
    return(1)
  } else {
    return(WSS / TSS)
  }
}

calc_wss_tss_ratio(scaled_data, individual = individuals[[1]])



##########################################################
# expose_and_evaluate_wss
# - Purpose: For each individual, sample a portion of the data
#   if exposure_prop < 1, then compute the WSS/TSS ratio.
# - Returns a numeric vector of fitness values (ratio).
##########################################################

calc_fitness <- function(data, individual, penalty = 0.065) {
  # ratio = WSS/TSS
  ratio <- calc_wss_tss_ratio(data, individual)
  
  # k
  k_val <- individual$k
  
  # Weighted sum: ratio + penalty * k
  # You can tune 'penalty' to suit your data scale and preference
  return(ratio + penalty * k_val)
}

expose_and_evaluate_wss <- function(population, data, exposure_prop = 0.1, penalty = 0.05) {
  N <- nrow(data)
  n_sample <- ceiling(exposure_prop * N)
  
  fitness_values <- numeric(length(population))
  
  for (i in seq_along(population)) {
    # sample subset
    sampled_rows <- sample(seq_len(N), n_sample, replace = TRUE)
    sub_data <- data[sampled_rows, , drop = FALSE]
    
    # compute fitness (ratio + penalty * k)
    fitness_values[i] <- calc_fitness(sub_data, population[[i]], penalty = penalty)
  }
  
  return(fitness_values)
}



##########################################################
# select_population
# - Purpose: Select the best-performing individuals 
#   based on their fitness (here, the ratio WSS/TSS).
# - smaller is better => we sort ascending, keep top `top_prop`.
##########################################################
select_population <- function(population, fitness_values, top_prop = 0.5) {
  pop_size <- length(population)
  n_keep   <- ceiling(top_prop * pop_size)
  
  # sort indices by ascending ratio
  sorted_idx <- order(fitness_values, decreasing = FALSE)
  
  # keep top n_keep
  keep_idx   <- sorted_idx[seq_len(n_keep)]
  selected   <- population[keep_idx]
  
  return(selected)
}



##########################################################
# select_population
# - Purpose: Select the best-performing individuals 
#   based on their fitness (here, the ratio WSS/TSS).
# - smaller is better => we sort ascending, keep top `top_prop`.
##########################################################
select_population <- function(population, fitness_values, top_prop = 0.5) {
  pop_size <- length(population)
  n_keep   <- ceiling(top_prop * pop_size)
  
  # sort indices by ascending ratio
  sorted_idx <- order(fitness_values, decreasing = FALSE)
  
  # keep top n_keep
  keep_idx   <- sorted_idx[seq_len(n_keep)]
  selected   <- population[keep_idx]
  
  return(selected)
}
##########################################################
# reproduce_dna
# - Purpose: Combine "DNA" from two parents to produce 
#   an offspring DNA vector.
##########################################################
reproduce_dna <- function(parent1, parent2, method = c("random", "average")) {
  method <- match.arg(method)
  dna1 <- parent1$dna
  dna2 <- parent2$dna
  
  # ensure both have same length
  stopifnot(length(dna1) == length(dna2))
  
  if (method == "random") {
    # random gene-by-gene
    mask <- runif(length(dna1)) < 0.5
    child_dna <- ifelse(mask, dna1, dna2)
  } else if (method == "average") {
    # simple average
    child_dna <- (dna1 + dna2) / 2
  }
  
  return(child_dna)
}

##########################################################
# reproduce_dna_cluster_level
# - Purpose: Combine parent's DNA at the *cluster* level 
#            rather than gene-by-gene.
# - method options:
#   "random_cluster"  : For each cluster, entirely use 
#                       that cluster's center from Parent 1 or Parent 2
#   "average_cluster" : For each cluster, average the entire cluster center
#   "random_gene"     : (optional) The old gene-by-gene approach
##########################################################
reproduce_dna_cluster_level <- function(parent1, parent2, 
                                        method = c("random_cluster", 
                                                   "average_cluster", 
                                                   "random_gene")) {
  method <- match.arg(method)
  
  dna1 <- parent1$dna
  dna2 <- parent2$dna
  
  stopifnot(length(dna1) == length(dna2))  # both parents must have same DNA length
  stopifnot(parent1$k == parent2$k)        # same k (number of clusters)
  
  k <- parent1$k
  total_len <- length(dna1)
  # total_len = k * d, where d is number of dimensions
  d <- total_len / k
  if ((d %% 1) != 0) {
    stop("DNA length not divisible by k—check your input.")
  }
  d <- as.integer(d)
  
  # We'll store child's DNA in this vector
  child_dna <- numeric(total_len)
  
  # Reshape dna1, dna2 into k x d matrices (one row per cluster)
  mat1 <- matrix(dna1, nrow = k, ncol = d, byrow = TRUE)
  mat2 <- matrix(dna2, nrow = k, ncol = d, byrow = TRUE)
  
  if (method == "random_cluster") {
    # For each cluster c, pick entire row from parent1 or parent2
    for (c_idx in seq_len(k)) {
      if (runif(1) < 0.5) {
        child_dna[((c_idx-1)*d + 1):(c_idx*d)] <- mat1[c_idx, ]
      } else {
        child_dna[((c_idx-1)*d + 1):(c_idx*d)] <- mat2[c_idx, ]
      }
    }
    
  } else if (method == "average_cluster") {
    # For each cluster c, average entire row from parent1 and parent2
    for (c_idx in seq_len(k)) {
      avg_center <- (mat1[c_idx, ] + mat2[c_idx, ]) / 2
      child_dna[((c_idx-1)*d + 1):(c_idx*d)] <- avg_center
    }
    
  } else if (method == "random_gene") {
    # Old gene-by-gene approach
    mask <- runif(total_len) < 0.5
    child_dna <- ifelse(mask, dna1, dna2)
  }
  
  return(child_dna)
}


##########################################################
# mate_and_reproduce
# - Purpose: Among the selected individuals, 
#   randomly pair them (within the same k),
#   produce offspring, and add random mutations.
##########################################################

mate_and_reproduce <- function(
    selected_pop,
    reproduction_method = "random_cluster",  
    random_mutation_sd = 0.1,
    offspring_delta = 2
) {
  # Group by k
  same_k <- split(selected_pop, sapply(selected_pop, function(x) x$k))
  
  new_generation <- list()
  
  for (k_val in names(same_k)) {
    group <- same_k[[k_val]]
    
    if (length(group) < 2) {
      new_generation <- c(new_generation, group)
      next
    }
    
    # Shuffle for random pairing
    shuffle_idx <- sample(seq_along(group))
    group <- group[shuffle_idx]
    
    # We'll do a manual loop to pair, skipping siblings
    used <- rep(FALSE, length(group))   # track who is used
    
    for (i in seq_along(group)) {
      if (used[i]) next  # already paired
      
      parent1 <- group[[i]]
      best_pair_index <- NA
      
      # find a partner that is not a sibling
      for (j in (i+1):length(group)) {
        if (j > length(group)) break  # or next
        if (!used[j]) {
          parent2 <- group[[j]]
          
          # Check if they share the same parents => siblings
          siblings <- FALSE
          if (!is.null(parent1$parents) && !is.null(parent2$parents)) {
            # If any overlap in parents => siblings
            # (depending on your definition, you may require both parents same, or just one, etc.)
            if (length(intersect(parent1$parents, parent2$parents)) > 0) {
              siblings <- TRUE
            }
          }
          
          if (!siblings) {
            best_pair_index <- j
            break
          }
        }
      }
      
      if (is.na(best_pair_index)) {
        # no non-sibling found, so parent1 remains unpaired
        new_generation <- c(new_generation, list(parent1))
        used[i] <- TRUE
      } else {
        # pair found
        parent2 <- group[[best_pair_index]]
        used[i] <- TRUE
        used[best_pair_index] <- TRUE
        
        # produce offspring
        n_offspring <- rpois(1, lambda = offspring_delta)
        if (n_offspring < 1) n_offspring <- 1
        
        for (o in seq_len(n_offspring)) {
          child_dna <- reproduce_dna_cluster_level(
            parent1,
            parent2,
            method = reproduction_method
          )
          
          # mutation
          child_dna <- child_dna + rnorm(length(child_dna), mean = 0, sd = random_mutation_sd)
          child_dna <- pmax(pmin(child_dna, 1), 0)
          
          # create new child => store parents, assign unique id
          # For example, just use temp id from global counter or
          # store gen in the environment, etc.
          child_id <- paste0("child_", sample.int(999999, 1))  
          
          child <- list(
            k       = parent1$k,
            dna     = child_dna,
            parents = c(parent1$id, parent2$id),
            id      = child_id
          )
          
          new_generation <- c(new_generation, list(child))
        }
      }
    }
  }
  
  return(new_generation)
}




##########################################################
# evolve_one_generation
# - Purpose: Conduct one generation's evolution.
#   Steps:
#    1) Evaluate each individual -> get WSS/TSS ratio
#    2) Select top performers
#    3) Mate & Reproduce -> new_generation
##########################################################
evolve_one_generation <- function(
    population,
    data,
    exposure_prop      = 0.1,
    top_prop           = 0.5,
    reproduction_method= "random",
    random_mutation_sd = 0.1,
    offspring_delta    = 2
) {
  # 1) Evaluate -> get WSS/TSS ratio
  fitness_values <- expose_and_evaluate_wss(population, data, exposure_prop)
  
  # 2) Select
  selected_pop <- select_population(population, fitness_values, top_prop)
  
  # 3) Mate & Reproduce
  new_generation <- mate_and_reproduce(
    selected_pop          = selected_pop,
    reproduction_method   = reproduction_method,
    random_mutation_sd    = random_mutation_sd,
    offspring_delta       = offspring_delta
  )
  
  return(new_generation)
  
  #STORE CLUSTER CENTERS #### 
  # For each individual in population:
  #   1) read their k, 
  #   2) reshape the dna into (k x n_dims)
  #   3) store a row for each cluster & dimension
  
  store_cluster_centers <- function(population, generation) {
    # We'll build a list of rows:
    rows <- list()
    idx <- 1
    
    for (i in seq_along(population)) {
      indiv   <- population[[i]]
      k       <- indiv$k
      dna     <- indiv$dna
      n_dims  <- length(dna) / k
      
      mat <- matrix(dna, nrow = k, ncol = n_dims, byrow = TRUE)
      
      for (c_idx in seq_len(k)) {
        # if we do wide format, we might do:
        # row: (generation, individual=i, cluster=c_idx, dim_1=..., dim_2=..., etc.)
        # but let's show a LONG approach:
        for (d_idx in seq_len(n_dims)) {
          rows[[idx]] <- data.frame(
            generation = generation,
            individual = i,
            cluster    = c_idx,
            dimension  = d_idx,
            value      = mat[c_idx, d_idx]
          )
          idx <- idx + 1
        }
      }
    }
    
    # combine into one data frame
    df <- do.call(rbind, rows)
    return(df)
  }
  
}


##########################################################
# compute_wss_bss_tss
# - Purpose: For a single 'individual' (which has k, dna),
#   compute WSS, TSS, BSS, ratio = WSS/TSS, and R^2 = BSS/TSS.
#
#   data: matrix or data.frame of size [n x d]
#   individual: list with:
#       $k   (number of clusters)
#       $dna (flattened cluster centers, length = k*d)
#
#   Returns a list:
#   list(
#     WSS   = ...,
#     TSS   = ...,
#     BSS   = ...,
#     ratio = WSS/TSS,
#     r2    = BSS/TSS
#   )
##########################################################
compute_wss_bss_tss <- function(data, individual) {
  # Convert data to matrix if needed
  data <- as.matrix(data)
  n_rows <- nrow(data)
  n_dims <- ncol(data)
  
  k <- individual$k
  dna <- individual$dna
  
  # Reshape dna -> centers matrix [k x d]
  centers <- matrix(dna, nrow = k, ncol = n_dims, byrow = TRUE)
  
  # 1) Assign each row to nearest center
  dist_mat <- sapply(seq_len(k), function(ci) {
    rowSums((data - 
               matrix(centers[ci, ], nrow = n_rows, ncol = n_dims, byrow = TRUE))^2)
  })
  # for each point, the cluster with minimum distance
  cluster_assign <- max.col(-dist_mat)
  
  # 2) Compute TSS = sum of squared distances of data points to the overall mean
  overall_mean <- colMeans(data)
  TSS <- sum(rowSums((data - matrix(overall_mean, 
                                    nrow = n_rows, 
                                    ncol = n_dims, 
                                    byrow = TRUE))^2))
  
  # 3) Compute WSS by summing the (point - cluster_center)^2
  WSS <- 0
  for (c_idx in seq_len(k)) {
    # subset of data for cluster c_idx
    is_c <- (cluster_assign == c_idx)
    if (!any(is_c)) next  # if no points assigned, skip
    sub_data <- data[is_c, , drop = FALSE]
    center_c <- centers[c_idx, ]
    
    WSS_c <- sum(rowSums((sub_data - 
                            matrix(center_c, nrow = nrow(sub_data), 
                                   ncol = ncol(sub_data), byrow = TRUE))^2))
    WSS <- WSS + WSS_c
  }
  
  # 4) BSS = TSS - WSS
  BSS <- TSS - WSS
  
  # 5) ratio = WSS/TSS (can be >1 if centers are far from data, especially if random)
  ratio_val <- if (TSS < 1e-12) 1 else WSS / TSS
  
  # 6) R^2 = BSS / TSS
  r2_val <- if (TSS < 1e-12) 0 else BSS / TSS
  
  return(list(
    WSS   = WSS,
    TSS   = TSS,
    BSS   = BSS,
    ratio = ratio_val,
    r2    = r2_val
  ))
}


#######STORE CLUSTERSSTATS ############

store_partition_stats <- function(population, data, generation) {
  # We'll loop over each individual, compute WSS/TSS, BSS, R2, etc.
  # Return a data frame with one row per individual
  rows <- vector("list", length(population))
  
  for (i in seq_along(population)) {
    indiv <- population[[i]]
    ratio <- calc_wss_tss_ratio(data, indiv)       # already computing WSS/TSS
    
    # if you'd like to break out wss, tss, etc.:
    stats <- compute_wss_bss_tss(data, indiv)      # you'd define a function
    # where compute_wss_bss_tss() returns a list: list(WSS=..., TSS=..., BSS=..., ratio=...)
    # or you can adapt calc_wss_tss_ratio() to return more details.
    
    # example:
    wss <- stats$WSS
    tss <- stats$TSS
    bss <- stats$BSS
    ratio_val <- stats$ratio  # or wss/tss
    r2_val    <- bss / tss
    
    rows[[i]] <- data.frame(
      generation = generation,
      individual = i,
      k          = indiv$k,
      wss        = wss,
      bss        = bss,
      tss        = tss,
      ratio      = ratio_val,
      r2         = r2_val
    )
  }
  
  df <- do.call(rbind, rows)
  return(df)
}



##########################################################
# evo_clustering_wss_ratio
# - Purpose: Perform evolutionary clustering using 
#   WSS/TSS ratio as the fitness measure.
#
#   data             = your dataset
#   scale            = whether to scale data to [0,1]
#   min_centers, max_centers = range for #clusters
#   n_gen_1          = population size of the 1st generation
#   n_gen            = total generations to evolve
#   exposure_prop    = fraction of data to use for each individual
#   top_prop         = fraction of population to keep each gen
#   reproduction_method = how to combine DNA ("random" or "average")
#   random_mutation_sd = mutation std dev
#   offspring_delta  = lambda for Poisson #offspring
#   verbose          = print progress messages
##########################################################
evo_clustering_wss_ratio <- function(
    data,
    scale               = TRUE,
    min_centers         = 2,
    max_centers         = 5,
    n_gen_1             = 300,
    n_gen               = 50,
    exposure_prop       = 0.2,
    top_prop            = 0.1,
    reproduction_method = "random_cluster",  # e.g., "average_cluster" or "random_gene"
    random_mutation_sd  = 0.1,
    offspring_delta     = 15,
    verbose             = TRUE
) {
  # 1) Scale the data if requested
  data <- scale_data(data, scale)
  
  # 2) Initialize population
  population <- initialize_population(
    pop_size    = n_gen_1,
    min_centers = min_centers,
    max_centers = max_centers,
    n_dims      = ncol(data)
  )
  
  # 3) Prepare data structures to store diagnostics each generation
  cluster_dist_list        <- vector("list", length = n_gen)  # distribution of k
  cluster_centers_df_list  <- vector("list", length = n_gen)  # cluster centers
  part_stats_df_list       <- vector("list", length = n_gen)  # wss, bss, tss, ratio, etc.
  
  # 4) Evolve for n_gen generations
  for (g in seq_len(n_gen)) {
    # Evaluate current population -> WSS/TSS ratio
    fitness_values <- expose_and_evaluate_wss(population, data, exposure_prop)
    
    # Select top individuals
    selected_pop <- select_population(population, fitness_values, top_prop)
    
    # Mate & Reproduce => new generation
    population <- mate_and_reproduce(
      selected_pop          = selected_pop,
      reproduction_method   = reproduction_method,
      random_mutation_sd    = random_mutation_sd,
      offspring_delta       = offspring_delta
    )
    
    # 4a) Distribution of k
    k_vals <- sapply(population, function(ind) ind$k)
    k_tab  <- table(k_vals)  
    dist_k_df <- data.frame(
      generation = g,
      k = as.integer(names(k_tab)),
      count = as.integer(k_tab)
    )
    cluster_dist_list[[g]] <- dist_k_df
    
    # 4b) Store cluster centers in a tidy data frame
    # (assuming store_cluster_centers(population, g) returns a df with columns 
    #  like: generation, individual, cluster, dimension, value)
    cluster_centers_df_list[[g]] <- store_cluster_centers(population, g)
    
    # 4c) Store partition stats: WSS, TSS, ratio, etc.
    # (assuming store_partition_stats(population, data, g) returns a df 
    #  with columns: generation, individual, k, wss, bss, tss, ratio, r2, ...)
    part_stats_df_list[[g]] <- store_partition_stats(population, data, g)
    
    # Print distribution of k
    cat(sprintf("Gen %d | Distribution of k: %s\n", 
                g, 
                paste(names(k_tab), k_tab, sep="=", collapse=", ")))
    
    # (Optional) Print best & avg ratio
    if (verbose) {
      best_ratio <- min(fitness_values)
      avg_ratio  <- mean(fitness_values)
      cat(sprintf(
        "Generation %d | Pop size: %d | Best ratio: %.4f | Avg ratio: %.4f\n\n",
        g, length(population), best_ratio, avg_ratio
      ))
    }
  }
  
  # 5) Combine all generational data frames into one
  cluster_dist_df    <- do.call(rbind, cluster_dist_list)
  cluster_centers_df <- do.call(rbind, cluster_centers_df_list)
  part_stats_df      <- do.call(rbind, part_stats_df_list)
  
  # 6) Return final population + all diagnostics
  return(list(
    final_population     = population,
    cluster_dist_df      = cluster_dist_df,
    cluster_centers_df   = cluster_centers_df,
    part_stats_df        = part_stats_df
  ))
}





evo_results <- evo_clustering_wss_ratio(n_gen = 100, 
                                        data = data,
                                        offspring_delta = 20, 
                                        n_gen_1 = 1000,
                                        top_prop = 0.1,
                                        exposure_prop = 0.2,
                                        random_mutation_sd = 0.1,
                                        min_centers = 2,
                                        max_centers = 5)



########EXPLORING THE RESULTS#########


library(tidyverse)

evo_results$cluster_centers_df %>%
  group_by(individual) %>%
  mutate(clusters = max(cluster)) %>% 
  group_by(clusters, dimension, cluster) %>% 
  summarise(median = median(value),
            lower_89 = quantile(value, .055),
            upper_89 = quantile(value, .945)
  ) %>% 
  ggplot(aes(x = median, xmin = lower_89,
             xmax = upper_89,
             y = dimension, color = factor(cluster))) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  #geom_point(position = position_dodge(width = 0.5)) +
  facet_wrap( . ~ clusters, scales = "free", ncol = 1)



evo_results$part_stats_df %>% 
  ggplot(aes(x = generation, y = r2)) +
  geom_smooth()



evo_results$cluster_centers_df %>%
  group_by(individual) %>%
  mutate(clusters = max(cluster)) %>% 
  group_by(clusters, dimension, cluster) %>% 
  summarise(median = median(value),
            lower_89 = quantile(value, .055),
            upper_89 = quantile(value, .945)
  ) %>% 
  filter(clusters == 3) %>% 
  ggplot(aes(x = median, xmin = lower_89,
             xmax = upper_89,
             y = dimension, color = factor(cluster))) +
  geom_pointrange(position = position_dodge(width = 0.5))



scaled_data %>%
  tibble(Species = iris$Species) %>% 
  pivot_longer(!Species, names_to = "dimension") %>% 
  group_by(Species, dimension) %>% 
  summarise(median = median(value),
            lower_89 = quantile(value, .055),
            upper_89 = quantile(value, .945)
  ) %>% 
  ggplot(aes(x = median, xmin = lower_89,
             xmax = upper_89,
             y = dimension, color = factor(Species))) +
  geom_pointrange(position = position_dodge(width = 0.5))
