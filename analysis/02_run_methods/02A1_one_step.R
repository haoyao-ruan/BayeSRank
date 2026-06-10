#--------------------------
memory.limit(size=NA)
library(data.table)  # Faster data handling
library(doParallel)
library(foreach)
library(progressr)
handlers(global = TRUE)

library(RobustRankAggreg)
library(TopKLists)
source("00_functions.r")
source("biger.r")
source("BIRRA.r")

# source("BayeSRank_Rcpp.r")
# sourceCpp("BayeSRank.cpp")

load("RData/01_default_simulated_1000_seed42.RData")


# Detect available cores
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
getDoParWorkers()

methods_to_compare <- c("BayeSRank", "Mean", "Geomean", "L2norm",
                        "MC1", "MC2", "MC3", "CEMC.s", "CEMC.k",
                        "Stuart", "RRA", "min", "BIRRA", "BiGER")

M_itrns <- 20000
m_burn <- 10000
to_drop <-1:m_burn

#(sdb_worst=qnorm(1 - (1 / (2 * 50)))/2)

# Store results
n_iter <- nrow(params) * n_Sims

set.seed(42)   
results_list <- foreach( i = 1:n_iter,  

                        #.noexport = c('BayeSRank_Rcpp'),
                        .packages = c('data.table', 
                                      'stats', 'expm', 
                                      'MCMCpack', 'extraDistr', 
                                      'TopKLists', 'RobustRankAggreg')) %dopar% {
                                
                                        if (i %% 1000 == 0) {
                                          cat(sprintf("Finished %d/%d at %s\n", 
                                                      i, n_iter, format(Sys.time(), "%H:%M:%S")),
                                              file = "progress_log.txt", append = TRUE)
                                        }
                                        
  data_info <- all_datasets[[i]]
  n_teams <- data_info$n_teams
  true_rank <- data_info$true_rank
  rank_matrix <- data_info$rank
  rank_list <- data_info$rank_list
  itemnames <- rownames(rank_matrix)
  
  # Initialize storage
  corr_results <- data.table(n_teams = n_teams, 
                               rho = data_info$rho, 
                               muBeta = data_info$muBeta,
                               sigmaBeta = data_info$sigmaBeta,
                               sim_id = data_info$sim_id)
  top1_results <- data.table()
  top3_results <- data.table()
  
  # Run all methods on the **same dataset**
  # # #BayeSRank------------------------------
  (W<-initW <- matrix(rnorm(n_teams^2, mean = 0, sd = 1), n_teams, n_teams))  # Initial W with no order

  # Step 1: Sort each column of W in decreasing order
  sorted_W <- apply(W, 2, sort, decreasing = TRUE)

  # Step 2: Use mapply to map sorted_W to initW based on mydat$rank
  initW <- mapply(function(sorted_col, rank_col) {
    sorted_col[rank_col]  # Map sorted values using the ranking index
  }, as.data.frame(sorted_W), as.data.frame(rank_matrix))

  # Convert back to a matrix (if needed)
  initW <- matrix(unlist(initW), nrow = nrow(W), ncol = ncol(W))


  res_bayesrank <- BayeSRank.pre(n = n_teams,
                              r = rank_matrix,
                              M = M_itrns,
                              burnin = m_burn,
                              W = initW,
                              mu0 = rnorm(n_teams,0,1),
                              mubeta0 = rnorm(1,0,1),
                              Var_beta0 = 1,
                              Var_epsilon0 = 1,
                              minbeta = -10,
                              maxbeta = 10,
                              a=0.1,
                              b=0.1,
                              c=0.1,
                              d=0.1)


    corr_results[, BayeSRank := cor(true_rank, rank(-res_bayesrank$post_mean_mus), method="spearman")]

    top1_results[, top1_BayeSRank := calculate_top1(rank(-res_bayesrank$post_mean_mus), true_rank)]
    top3_results[, top3_BayeSRank := calculate_top3(rank(-res_bayesrank$post_mean_mus), true_rank)]


  # # Append results-------
  combined_result <- cbind(corr_results, 
                           top1_results, 
                           top3_results)
}

# Stop cluster
stopCluster(cl)


# Convert results to data.table
results <- rbindlist(results_list)

