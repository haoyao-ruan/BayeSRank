#----------------------
memory.limit(size=NA)
source("00_functions.r")
library(data.table)
library(foreach)
library(doParallel)

# Detect available CPU cores (reserve 1 for system processes)
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Define parameter values
n_values <- c(5, 10, 15, 20, 25)
rho_values <- c(0.15, 0.35, 0.55, 0.75, 0.95)
muBeta_values <- c(-1.5, -1,-0.5,  0, 0.5,  1, 1.5, 5.0)
sigmaBeta_values <- sqrt(c(0.25, 0.5, 1, 2, 4))

# Set default values
default_params <- data.frame(n_teams = 15, rho = 0.55, muBeta = 1, sigmaBeta=sqrt(0.5))

# Varying one parameter at a time
df_n <- data.frame(n_teams= n_values, 
                   rho = default_params$rho, 
                   muBeta = default_params$muBeta, sigmaBeta=default_params$sigmaBeta)

df_rho <- data.frame(n_teams= default_params$n, 
                     rho = rho_values, 
                     muBeta = default_params$muBeta, sigmaBeta=default_params$sigmaBeta)

df_muBeta <- data.frame(n_teams= default_params$n, 
                        rho = default_params$rho, 
                        muBeta = muBeta_values, sigmaBeta=default_params$sigmaBeta)

df_sigmaBeta <- data.frame(n_teams= default_params$n, 
                        rho = default_params$rho, 
                        muBeta = default_params$muBeta, sigmaBeta=sigmaBeta_values)
# Combine all variations into one data frame
params <- unique(rbind(df_n, 
                       df_rho, 
                       df_muBeta,
                       df_sigmaBeta))%>%
  mutate(sigmaj = sqrt(find_sigmas(rho)$sig2_eps))#sigma_beta = sigma_epsilon by assumption

params <- as.data.table(params)

n_Sims <- 1000  # Number of datasets per parameter setting

set.seed(42)
# Run simulations in parallel


n_iter <- nrow(params)

  
all_datasets <- foreach(i = 1:n_iter, 
                        .combine = 'c', 
                        .packages = c('data.table')) %dopar% {
                          
    datasets_list <- vector("list", n_Sims)
    
  for (s in 1:n_Sims) {
    
    
    mydat <- gen_data_dir(n = params$n_teams[i], 
                          mu_beta = params$muBeta[i],
                          sigma_beta=params$sigmaBeta[i], 
                          sigma_epsilon = params$sigmaj[i])
    
    datasets_list[[s]] <- list(
      sim_id = s,
      n_teams = params$n_teams[i], 
      rho = params$rho[i], 
      muBeta = params$muBeta[i],
      sigmaBeta = params$sigmaBeta[i],
      true_betas = mydat$true_beta, 
      rank = mydat$rank, 
      true_rank = mydat$true_rank,
      rank_list = mydat$rank_list
    )
  }
  
  return(datasets_list)  # Return the list of datasets

}


# Stop parallel cluster
stopCluster(cl)

# Flatten the nested list structure before saving
#all_datasets <- unlist(all_datasets, recursive = FALSE)

save.image("RData/01_default_simulated_1000_seed42.RData")

print("Data simulation complete and saved to simulated_datasets.rds.")
