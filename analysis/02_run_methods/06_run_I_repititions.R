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
  
  # # Run all methods on the **same dataset**
  # #BayeSRank------------------------------
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

  # Append results-------
  combined_result <- cbind(corr_results, 
                           top1_results, 
                           top3_results) 
}

# Stop cluster
stopCluster(cl)


# Convert results to data.table
results_1000 <- rbindlist(results_list)


# results_200 <- rbindlist(results_list)
# 
# summ50 <- results_50%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%
#   summarise(mean50=mean(BayeSRank))
# 
# summ100 <- results_100%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%
#   summarise(mean100=mean(BayeSRank))
# 
# summ200 <- results_200%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%
#   summarise(mean200=mean(BayeSRank))
# 
#  summ500 <- results_500%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%
#    summarise(mean500=mean(BayeSRank))
# 
#  summ1000 <- results_1000%>% group_by(n_teams, rho, muBeta, sigmaBeta)%>%
#    summarise(mean1000=mean(BayeSRank))
# 
# merge(summ50,summ100) %>% 
#   merge(summ200)%>% merge(summ500)%>% merge(summ1000)%>%
#   mutate(mean50=mean50*100,
#          mean100=mean100*100,
#          mean200=mean200*100,
#          mean500=mean500*100,
#          mean1000=mean1000*100)%>%
#   arrange(n_teams, rho, muBeta)->a
# 
