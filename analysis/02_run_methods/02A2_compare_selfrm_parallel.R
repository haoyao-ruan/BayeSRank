#--------------------------
memory.limit(size=NA)
library(data.table)  # Faster data handling
library(foreach)
library(doParallel)
library(RobustRankAggreg)
library(TopKLists)
pacman::p_load('Rcpp', 'RcppArmadillo')

source("00_functions.r")
source("biger.r")
source("BIRRA.r")

# source("BayeSRank_Rcpp.r")
# sourceCpp("BayeSRank.cpp")

load("RData/01_default_simulated_1000_seed42.RData")
#load("RData/01E_default_simulated_datasets_mix.RData")

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
                                        
results_list.selfrm <- foreach( i = 1:(nrow(params)*n_Sims),  
                        #.noexport = c('BayeSRank_Rcpp'),
                        .packages = c('data.table', 
                                      'stats', 'expm', 
                                      'MCMCpack', 'extraDistr', 
                                      'TopKLists', 'RobustRankAggreg',
                                      'Rcpp', 'RcppArmadillo')) %dopar% {
                                
                                    #    source("BayeSRank_Rcpp.r")
                                        
  data_info <- all_datasets[[i]]
  n_teams <- data_info$n_teams
  
  rank_matrix <- data_info$rank
  
  diag(rank_matrix) <- NA
  rank_list <- convert_to_rank_lists(rank_matrix )
  
  true_rank <- data_info$true_rank
  
  itemnames <- rownames(rank_matrix)
  
  # Initialize storage
  corr_results <- data.table(n_teams = n_teams,
                               rho = data_info$rho,
                               muBeta = data_info$muBeta,
                             sigmaBeta = data_info$sigmaBeta,
                               sim_id = data_info$sim_id)
  top1_results <- data.table()
  top3_results <- data.table()
   
  # # # Run all methods on the **same dataset**
  #------------------------
  res_Borda <- Borda(rank_list)[[1]]

  corr_results[, Mean.selfrm :=   cor(true_rank, match(itemnames, res_Borda$mean), method="spearman")]
  corr_results[, Geomean.selfrm :=cor(true_rank, match(itemnames, res_Borda$geo.mean), method="spearman")]
  corr_results[, Median.selfrm := cor(true_rank, match(itemnames, res_Borda$median), method="spearman")]
  corr_results[, L2norm.selfrm :=cor(true_rank, match(itemnames, res_Borda$l2norm), method="spearman")]

  top1_results[, top1_Mean.selfrm :=calculate_top1(match(itemnames, res_Borda$mean), true_rank)]
  top1_results[, top1_Geomean.selfrm :=calculate_top1(match(itemnames, res_Borda$geo.mean), true_rank)]
  top1_results[, top1_Median.selfrm :=calculate_top1(match(itemnames, res_Borda$median), true_rank)]
  top1_results[, top1_L2norm.selfrm :=calculate_top1(match(itemnames, res_Borda$l2norm), true_rank)]

  top3_results[, top3_Mean.selfrm :=calculate_top3(match(itemnames, res_Borda$mean), true_rank)]
  top3_results[, top3_Geomean.selfrm :=calculate_top3(match(itemnames, res_Borda$geo.mean), true_rank)]
  top3_results[, top3_Median.selfrm :=calculate_top3(match(itemnames, res_Borda$median), true_rank)]
  top3_results[, top3_L2norm.selfrm :=calculate_top3(match(itemnames, res_Borda$l2norm), true_rank)]

  #------------------------
  res_mc <- MC(rank_list)
  corr_results[, MC1.selfrm :=cor(true_rank, match(itemnames, res_mc$MC1.TopK), method = "spearman")]
  corr_results[, MC2.selfrm :=cor(true_rank, match(itemnames, res_mc$MC2.TopK), method = "spearman")]
  corr_results[, MC3.selfrm :=cor(true_rank, match(itemnames, res_mc$MC3.TopK), method = "spearman")]

  top1_results[, top1_MC1.selfrm :=calculate_top1(match(itemnames, res_mc$MC1.TopK), true_rank)]
  top1_results[, top1_MC2.selfrm :=calculate_top1(match(itemnames, res_mc$MC2.TopK), true_rank)]
  top1_results[, top1_MC3.selfrm :=calculate_top1(match(itemnames, res_mc$MC3.TopK), true_rank)]
  
  top3_results[, top3_MC1.selfrm :=calculate_top3(match(itemnames, res_mc$MC1.TopK), true_rank)]
  top3_results[, top3_MC2.selfrm :=calculate_top3(match(itemnames, res_mc$MC2.TopK), true_rank)]
  top3_results[, top3_MC3.selfrm :=calculate_top3(match(itemnames, res_mc$MC3.TopK), true_rank)]


  res_cemc_s <- CEMC(rank_list, dm = "s")$TopK
  res_cemc_k <- CEMC(rank_list, dm = "k")$TopK
  corr_results[, CEMC.s.selfrm :=cor(true_rank, match(itemnames, res_cemc_s), method = "spearman")]
  corr_results[, CEMC.k.selfrm :=cor(true_rank, match(itemnames, res_cemc_k), method = "spearman")]

  top1_results[, top1_CEMC.s.selfrm :=calculate_top1(match(itemnames, res_cemc_s), true_rank)]
  top1_results[, top1_CEMC.k.selfrm :=calculate_top1(match(itemnames, res_cemc_k), true_rank)]

  top3_results[, top3_CEMC.s.selfrm :=calculate_top3(match(itemnames, res_cemc_s), true_rank)]
  top3_results[, top3_CEMC.k.selfrm :=calculate_top3(match(itemnames, res_cemc_k), true_rank)]


  res_stuart <- aggregateRanks(glist = rank_list, method = "stuart")$Name
  res_rra <- aggregateRanks(glist = rank_list, method = "RRA")$Name
  res_min <- aggregateRanks(glist = rank_list, method = "min")$Name

  corr_results[, Stuart.selfrm :=cor(true_rank, match(itemnames, res_stuart), method = "spearman")]
  corr_results[, RRA.selfrm :=cor(true_rank, match(itemnames, res_rra), method = "spearman")]
  corr_results[, min.selfrm :=cor(true_rank, match(itemnames, res_min), method = "spearman")]

  top1_results[, top1_Stuart.selfrm :=calculate_top1(match(itemnames, res_stuart),true_rank)]
  top1_results[, top1_RRA.selfrm :=calculate_top1(match(itemnames, res_rra),true_rank)]
  top1_results[, top1_min.selfrm :=calculate_top1(match(itemnames, res_min),true_rank)]

  top3_results[, top3_Stuart.selfrm :=calculate_top3(match(itemnames, res_stuart),true_rank)]
  top3_results[, top3_RRA.selfrm :=calculate_top3(match(itemnames, res_rra),true_rank)]
  top3_results[, top3_min.selfrm :=calculate_top3(match(itemnames, res_min),true_rank)]


  # BiGER----------------
  # res_biger <- BiGER(r=rank_matrix,
  #                    W=initW,
  #                    Beta0=rnorm(n_teams, 0, 1),
  #                    tau2= 5.5,
  #                    sigma_s0=varEpsilon_est,
  #                    sigma_sb0 = varBeta_est,
  #                    a=a_est,
  #                    b=b_est,
  #                    c=c_est,
  #                    d=d_est,
  #                    mu0=rnorm(n_teams,0,1),
  #                    M = M_itrns, burnin = m_burn,
  #                    verbose=1000)
  # corr_results[, BiGER.selfrm :=cor(true_rank, rank(-res_biger$post.mean.mu), method="spearman")]
  # top1_results[, top1_BiGER.selfrm :=calculate_top1(rank(-res_biger$post.mean.mu),true_rank)]
  # top3_results[, top3_BiGER.selfrm :=calculate_top3(rank(-res_biger$post.mean.mu),true_rank)]

  # BIRRA-----------------
  (res_BIRRA <- BIRRA(data=rank_matrix,
                      prior = (3/n_teams),   # Top-k prior (e.g., 0.05 = top 5% are "positive")
                      num.bin = n_teams,   # Number of bins to group data for Bayes factor estimation
                      num.iter = 10  # Max number of iterations for the algorithm to converge
  ) )
  corr_results[, BIRRA.selfrm := cor(true_rank, res_BIRRA, method="spearman")]

  top1_results[, top1_BIRRA.selfrm := calculate_top1(res_BIRRA, true_rank)]
  top3_results[, top3_BIRRA.selfrm := calculate_top3(res_BIRRA, true_rank)]


  # Append results-------
  combined_result.selfrm <- cbind(corr_results, 
                           top1_results, 
                           top3_results)
}

# Stop cluster
stopCluster(cl)

# Convert results to data.table
results_selfrm <- rbindlist(results_list.selfrm)

# Convert results to data.table
# results.selfrm <- rbindlist(results_list.selfrm)
# results.selfrm[, .(avg_corr = mean(BIRRA.selfrm)), by = .(n_teams, rho, muBeta)]
# results.selfrm[, .(avg_top1 = mean(top1_BIRRA.selfrm)), by = .(n_teams, rho, muBeta)]
# results.selfrm[, .(avg_top3 = mean(top3_BIRRA.selfrm)), by = .(n_teams, rho, muBeta)]
