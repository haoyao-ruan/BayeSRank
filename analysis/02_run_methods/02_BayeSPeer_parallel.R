#--------------------------
memory.limit(size=NA)
library(data.table)  # Faster data handling
library(doParallel)
library(foreach)
library(progressr)
handlers(global = TRUE)

source("00_functions.r")
source("bayesrank_core.r")

# source("BayeSRank_Rcpp.r")
# sourceCpp("BayeSRank.cpp")

#load("RData/01_default_simulated_1000_seed42.RData")
load("RData/01E_default_simulated_datasets_mix.RData")


# Detect available cores
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
getDoParWorkers()

M_itrns <- 2000
m_burn <- 1000
to_drop <-1:m_burn

#(sdb_worst=qnorm(1 - (1 / (2 * 50)))/2)

# Store results
n_iter <- nrow(params) * n_Sims

set.seed(42)   
results_list <- foreach( i = 1:n_iter, 

                        #.noexport = c('BayeSRank_Rcpp'),
                        .packages = c('data.table', 
                                      'stats', 'expm')) %dopar% {
                                
                                        if (i %% 100 == 0) {
                                          cat(sprintf("Finished %d/%d at %s\n", 
                                                      i, n_iter, format(Sys.time(), "%H:%M:%S")),
                                              file = "progress_log1.txt", append = TRUE)
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

  
  res_bayesrank <- bayesrank_2step(rank_matrix = rank_matrix,
                  M_itrns = M_itrns,
                  m_burn = m_burn,
                  step1_a = 0.1, step1_b = 0.1,
                  step1_c = 0.1, step1_d = 0.1,
                  unit_sd= 3)
  
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

# IG_0001, IG_001, IG_1
# SD_2, SD_5, SD_10
# bound 100, 1000
# HC_05, HC_1, HC_2, HC_4

# results$df <- 25

results %>% 
  group_by(n_teams, rho, muBeta, sigmaBeta)%>% summarise(mean(BayeSRank))

results %>% 
  group_by(n_teams, rho, muBeta, sigmaBeta)%>% summarise(mean(top1_BayeSRank))

results %>% 
  group_by(n_teams, rho, muBeta, sigmaBeta)%>% summarise(mean(top3_BayeSRank))
