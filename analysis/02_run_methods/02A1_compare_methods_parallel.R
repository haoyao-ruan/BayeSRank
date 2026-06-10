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

#load("RData/01_default_simulated_1000_seed42.RData")
load("RData/01E_default_simulated_datasets_mix.RData")

# Detect available cores
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
getDoParWorkers()

methods_to_compare <- c("BayeSRank", "Mean", "Geomean", "L2norm",
                        "MC1", "MC2", "MC3", "CEMC.s", "CEMC.k",
                        "Stuart", "RRA", "min", "BIRRA", "BiGER")

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
                                      'stats', 'expm', 
                                      'MCMCpack', 'extraDistr', 
                                      'TopKLists', 'RobustRankAggreg')) %dopar% {
                                
                                        if (i %% 100 == 0) {
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


  prelim_samples <- BayeSRank.pre(n = n_teams,
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


  # prelim_samples$Beta->temp
  # median(temp);
  # sd(temp)

  (Beta.bar <- mean(prelim_samples$Beta_draws))
  Beta.sd <- sd(prelim_samples$Beta_draws)
  #(Beta.sd <- min(sd(prelim_samples$Beta_draws), sdb_worst))
  #(Beta.sd <- max(sd(prelim_samples$Beta_draws), sdb_worst/2))

  (minBeta_est <- Beta.bar - 3* Beta.sd)
  (maxBeta_est <- Beta.bar + 3* Beta.sd)

  (varBeta_est <- median(prelim_samples$varBeta_draws) )
  (varEpsilon_est <- median(prelim_samples$varEpsilon_draws) )

  (IG_Beta_est <- find_shape_rate(mean = median(prelim_samples$varBeta_draws),
                                 sd=sd(prelim_samples$varBeta_draws)) )

  a_est <- IG_Beta_est$shape
  b_est <- IG_Beta_est$rate

  (IG_Epsilon_est <- find_shape_rate(mean = median(prelim_samples$varEpsilon_draws),
                                     sd=sd(prelim_samples$varEpsilon_draws)) )

  c_est <- IG_Epsilon_est$shape
  d_est <- IG_Epsilon_est$rate

  res_bayesrank <- BayeSRank.pre(n = n_teams,
                                 r = rank_matrix,
                                 M = M_itrns,
                                 burnin = m_burn,
                                 W = initW,
                                 mu0 = rnorm(n_teams,0,1),
                                 mubeta0 = rnorm(1,0, 1),
                                 Var_beta0 = varBeta_est,
                                 Var_epsilon0 = varEpsilon_est,
                                 minbeta = minBeta_est,
                                 maxbeta = maxBeta_est,
                                 a=a_est,
                                 b=b_est,
                                 c=c_est,
                                 d=d_est)

  # plot(res_bayesrank$varBeta_draws[1,],type="l")
  # plot(res_bayesrank$varEpsilon_draws[1,],type="l")
  # plot(res_bayesrank$varEpsilon_draws[2,],type="l")
  # plot(res_bayesrank$varEpsilon_draws[3,],type="l")
  # plot(res_bayesrank$varEpsilon_draws[4,],type="l")
  #
  #
  # plot(res_bayesrank$muBeta_draws[1,],type="l")
  # plot(res_bayesrank$Beta_draws[1,],type="l")
  # plot(res_bayesrank$Beta_draws[2,],type="l")
  # plot(res_bayesrank$Beta_draws[3,],type="l")


    corr_results[, BayeSRank := cor(true_rank, rank(-res_bayesrank$post_mean_mus), method="spearman")]

    top1_results[, top1_BayeSRank := calculate_top1(rank(-res_bayesrank$post_mean_mus), true_rank)]
    top3_results[, top3_BayeSRank := calculate_top3(rank(-res_bayesrank$post_mean_mus), true_rank)]

#     #------------------------
#     res_Borda <- Borda(rank_list)[[1]]
# 
#     corr_results[, Mean :=    cor(true_rank, match(itemnames, res_Borda$mean), method="spearman")]
#     corr_results[, Geomean := cor(true_rank, match(itemnames, res_Borda$geo.mean), method="spearman")]
#     corr_results[, Median :=  cor(true_rank, match(itemnames, res_Borda$median), method="spearman")]
#     corr_results[, L2norm := cor(true_rank, match(itemnames, res_Borda$l2norm), method="spearman")]
# 
#     top1_results[, top1_Mean := calculate_top1(match(itemnames, res_Borda$mean), true_rank)]
#     top1_results[, top1_Geomean := calculate_top1(match(itemnames, res_Borda$geo.mean), true_rank)]
#     top1_results[, top1_Median := calculate_top1(match(itemnames, res_Borda$median), true_rank)]
#     top1_results[, top1_L2norm := calculate_top1(match(itemnames, res_Borda$l2norm), true_rank)]
# 
#     top3_results[, top3_Mean := calculate_top3(match(itemnames, res_Borda$mean), true_rank)]
#     top3_results[, top3_Geomean := calculate_top3(match(itemnames, res_Borda$geo.mean), true_rank)]
#     top3_results[, top3_Median := calculate_top3(match(itemnames, res_Borda$median), true_rank)]
#     top3_results[, top3_L2norm := calculate_top3(match(itemnames, res_Borda$l2norm), true_rank)]
# 
#   
# #  BiGER----------------
#   res_biger <- BiGER(r=rank_matrix,
#                      W=initW,
#                      Beta0=rnorm(n_teams, 0, 1),
#                      tau2= 5.5,
#                      sigma_s0=varEpsilon_est,
#                      sigma_sb0 = varBeta_est,
#                      a=a_est,
#                      b=b_est,
#                      c=c_est,
#                      d=d_est,
#                      mu0=rnorm(n_teams,0,1),
#                      M = M_itrns, burnin = m_burn,
#                      verbose=1000)
#   corr_results[, BiGER :=cor(true_rank, rank(-res_biger$post.mean.mu), method="spearman")]
#   top1_results[, top1_BiGER :=calculate_top1(rank(-res_biger$post.mean.mu),true_rank)]
#   top3_results[, top3_BiGER :=calculate_top3(rank(-res_biger$post.mean.mu),true_rank)]

  # #------------------------
  
  res_mc <- MC(rank_list)
  
  
  corr_results[, MC1 := cor(true_rank, match(itemnames, res_mc$MC1.TopK), method = "spearman")]
  # corr_results[, MC2 := cor(true_rank, match(itemnames, res_mc$MC2.TopK), method = "spearman")]
  # corr_results[, MC3 := cor(true_rank, match(itemnames, res_mc$MC3.TopK), method = "spearman")]

  top1_results[, top1_MC1 := calculate_top1(match(itemnames, res_mc$MC1.TopK), true_rank)]
  # top1_results[, top1_MC2 := calculate_top1(match(itemnames, res_mc$MC2.TopK), true_rank)]
  # top1_results[, top1_MC3 := calculate_top1(match(itemnames, res_mc$MC3.TopK), true_rank)]

  top3_results[, top3_MC1 := calculate_top3(match(itemnames, res_mc$MC1.TopK), true_rank)]
  # top3_results[, top3_MC2 := calculate_top3(match(itemnames, res_mc$MC2.TopK), true_rank)]
  # top3_results[, top3_MC3 := calculate_top3(match(itemnames, res_mc$MC3.TopK), true_rank)]


  # res_cemc_s <- CEMC(rank_list, dm = "s")$TopK
  # res_cemc_k <- CEMC(rank_list, dm = "k")$TopK
  # corr_results[, CEMC.s := cor(true_rank, match(itemnames, res_cemc_s), method = "spearman")]
  # corr_results[, CEMC.k := cor(true_rank, match(itemnames, res_cemc_k), method = "spearman")]
  # 
  # top1_results[, top1_CEMC.s := calculate_top1(match(itemnames, res_cemc_s), true_rank)]
  # top1_results[, top1_CEMC.k := calculate_top1(match(itemnames, res_cemc_k), true_rank)]
  # 
  # top3_results[, top3_CEMC.s := calculate_top3(match(itemnames, res_cemc_s), true_rank)]
  # top3_results[, top3_CEMC.k := calculate_top3(match(itemnames, res_cemc_k), true_rank)]
  # 
  # 
  # res_stuart <- aggregateRanks(glist = rank_list, method = "stuart")$Name
  # res_rra <- aggregateRanks(glist = rank_list, method = "RRA")$Name
  # res_min <- aggregateRanks(glist = rank_list, method = "min")$Name
  # 
  # corr_results[, Stuart := cor(true_rank, match(itemnames, res_stuart), method = "spearman")]
  # corr_results[, RRA := cor(true_rank, match(itemnames, res_rra), method = "spearman")]
  # corr_results[, min := cor(true_rank, match(itemnames, res_min), method = "spearman")]
  # 
  # top1_results[, top1_Stuart := calculate_top1(match(itemnames, res_stuart),true_rank)]
  # top1_results[, top1_RRA := calculate_top1(match(itemnames, res_rra),true_rank)]
  # top1_results[, top1_min := calculate_top1(match(itemnames, res_min),true_rank)]
  # 
  # top3_results[, top3_Stuart := calculate_top3(match(itemnames, res_stuart),true_rank)]
  # top3_results[, top3_RRA := calculate_top3(match(itemnames, res_rra),true_rank)]
  # top3_results[, top3_min := calculate_top3(match(itemnames, res_min),true_rank)]

  # # BIRRA-----------------
  # (res_BIRRA <- BIRRA(data=rank_matrix,
  #                    prior = (3/n_teams),   # Top-k prior (e.g., 0.05 = top 5% are "positive")
  #                    num.bin = n_teams,   # Number of bins to group data for Bayes factor estimation
  #                     num.iter = 10  # Max number of iterations for the algorithm to converge
  #                     ) )
  # corr_results[, BIRRA := cor(true_rank, res_BIRRA, method="spearman")]
  # 
  # top1_results[, top1_BIRRA := calculate_top1(res_BIRRA, true_rank)]
  # top3_results[, top3_BIRRA := calculate_top3(res_BIRRA, true_rank)]

  # # Append results-------
  combined_result <- cbind(corr_results, 
                           top1_results, 
                           top3_results)
}

# Stop cluster
stopCluster(cl)


# Convert results to data.table
results <- rbindlist(results_list)

results%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%summarise(mean(MC1))
results%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%summarise(mean(top1_MC1))
results%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%summarise(mean(top3_MC1))



# old_corr<-readRDS("20251025/Rdata/df_expanded_corr.rds")
# old_corr$mean_value = old_corr$mean_value -0.001
# old_top1<-readRDS("20251025/Rdata/df_expanded_top1.rds")
# old_top1$mean_value = old_top1$mean_value -0.001
# old_top3<-readRDS("20251025/Rdata/df_expanded_top3.rds")
# old_top3$mean_value = old_top3$mean_value -0.001