library(copula)
library(doParallel)

d <- 2 # dimension

# First try this
#sample_size <- c(20, 30)
#num_samples <- 3

sample_size <- c(100)
num_samples <- 1000
 
# Currently without the t copula, as the derivations of xv-CIC (approximation of leave-one-out cross-validation) for this copula are much more complicated.
cop_TYPES <-  c("clayton", "gumbel", "joe", "frank", "gaussian")

source("functions_information_criteria.R")
source("explicit_derivatives_xv_CIC.R")
source("simulation_scenario_information_criteria.R")

param_GRID <- construct_param_GRID(cop_TYPES, sample_size)
param_LIST <- split(param_GRID, seq(nrow(param_GRID)))

log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

### PARALLEL SIMULATIONS 

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Export necessary objects to workers
clusterExport(cl, c("sample_size", "num_samples", "d",
                    "copula_choice_by_tau", "copula_param_domain",
                    "copula_funct_choice", "simulate_cop_name_tau_sample_size"))

# Load required packages on workers
clusterEvalQ(cl, {
  library(copula)
  source("functions_information_criteria.R")
  source("simulation_scenario_information_criteria.R")
  source("explicit_derivatives_xv_CIC.R")
})

# Parallel processing
results_list <-  foreach(param_row = param_LIST, .packages = c("copula")) %dopar% {
  
  cop_name <- as.character(param_row$cop_name)
  true_tau <- as.numeric(param_row$tau)
  n <- as.integer(param_row$sample_size)
  
  # Create a log file for cop_name,tau and sample_size, which will help us track at each step that the simulation is being conducted
  log_file <- file.path(log_dir, paste0("cop_name_", cop_name, "_tau_", true_tau, "_sample_size_", n, ".csv"))
  
  column_names <- c("cop_name", "sample_size", "iteration", "true_param", "true_tau"
                    ,"AIC_clayton", "AIC_gumbel", "AIC_joe", "AIC_frank", "AIC_gaussian"
                    ,"xv1_clayton", "xv1_gumbel", "xv1_joe", "xv1_frank", "xv1_gaussian"
                    ,"xv_CIC_clayton", "xv_CIC_gumbel", "xv_CIC_joe", "xv_CIC_frank", "xv_CIC_gaussian"
                    
                    ,"xvK_2_clayton", "xvK_2_gumbel", "xvK_2_joe", "xvK_2_frank", "xvK_2_gaussian"
                    ,"xvK_5_clayton", "xvK_5_gumbel", "xvK_5_joe", "xvK_5_frank", "xvK_5_gaussian"
                    ,"xvK_10_clayton", "xvK_10_gumbel", "xvK_10_joe", "xvK_10_frank", "xvK_10_gaussian"
                    
                    ,"xvK_n:logn_clayton", "xvK_n:logn_gumbel", "xvK_n:logn_joe", "xvK_n:logn_frank", "xvK_n:logn_gaussian"
                    ,"xvK_n05_clayton", "xvK_n05_gumbel", "xvK_n05_joe", "xvK_n05_frank", "xvK_n05_gaussian"
                    ,"xvK_n05:logn_clayton", "xvK_n05:logn_gumbel", "xvK_n05:logn_joe", "xvK_n05:logn_frank", "xvK_n05:logn_gaussian"
                    ,"xvK_n:2_clayton", "xvK_n:2_gumbel", "xvK_n:2_joe", "xvK_n:2_frank", "xvK_n:2_gaussian"
                    )
  
  # Write headers to the CSV file
  write(paste(column_names, collapse = ","), file = log_file)
  
  cop_results <- simulate_cop_name_tau_sample_size(cop_name, cop_TYPES, true_tau, n, num_samples, d, log_file)
  return(cop_results)
}

# Stop cluster
stopCluster(cl)

