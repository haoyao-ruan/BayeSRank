
convert_to_rank_lists <- function(rank_matrix) {
  # 1) For each column, sort the ranks and extract the names 
  #    (i.e., the row names of `rank_matrix`).
  #    This returns a character matrix.
  sorted_mat <- apply(rank_matrix, 2, function(col) names(sort(col)))
  
  # 2) Convert each column of that matrix into an element of a list,
  #    so we get a list of length = number of columns.
  lapply(seq_len(ncol(sorted_mat)), function(j) sorted_mat[, j])
}

# # Generate rank lists
# rank_lists <- convert_to_rank_lists(mydat)

convert_to_ranking <- function(ordinal_list) {
  ranking <- seq_along(ordinal_list)  # Assign ranks based on order
  names(ranking) <- ordinal_list      # Assign item names for clarity
  return(ranking[order(names(ranking))]) # Reorder to match original item order
}

rank_diff_self_vs_peers <- function(rank_mat) {
  stopifnot(is.matrix(rank_mat),
            nrow(rank_mat) == ncol(rank_mat))
  
  # 1. Copy and blank out the diagonal
  peer_mat <- rank_mat
  diag(peer_mat) <- NA      # leave only peer ranks in each row
  
  # 2. Average peer ranks for each item (NA-aware)
  peer_means <- rowMeans(peer_mat, na.rm = TRUE)
  peer_ranks <- rank(peer_means)
  
  # 3. Difference: self rank – average peer rank
  diff1 <- peer_means - diag(rank_mat) 
  diff2 <- peer_ranks - diag(rank_mat) 
  
  # keep row names (items) if they exist
  if (!is.null(rownames(rank_mat))) names(diff) <- rownames(rank_mat)
  
  return(list(meandiff=diff1, rankdiff=diff2, peer_means=peer_means, peer_ranks=peer_ranks))
}



sort_matrix_columns <- function(n, r=myrank) {
  # Initialize W with random normal values
  initW <- W <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n, ncol = n)
  
  # Sorting each column of W and storing in initW
  for (i in 1:n) {
    initW[, i] <- sort(W[, i], decreasing = TRUE)[r[, i]]
  }
  
  return(initW)
}


find_sigmas <- function(rho, perc_bias =0.5, randomRho=FALSE){
  
  if (randomRho){
    rho = runif(1, 0.1, 0.9) 
  }
  
  sig2_beta =sig2_eps = (1/rho^2 - 1)
  
  return(list(sig2_beta=sig2_beta,
              sig2_eps=sig2_eps))
}
