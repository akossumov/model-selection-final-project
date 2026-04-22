# Load required libraries
library(copula)
library(numDeriv)
library(dplyr)
library(survival)
library(foreach)
library(doParallel)


# All references are made to the following papers 
# (1) "Copula information criterion for model selection with two-stage maximum likelihood estimation" 2019, by Vinnie Ko and Nils Lid Hjort
# (2) "Model robust inference with two-stage maximum likelihood estimation for copulas" 2019, by Vinnie Ko and Nils Lid Hjort

###### SIMULATION SETUP & TRUE MODEL PARAMETERS
n_obs <- 1000        # Sample size per simulation 
n_simulations <- 1  

# True Data Generating Model
true_theta <- 3
true_copula <- gumbelCopula(true_theta, dim = 4)

# Calculate the true 0.7-quantiles ($q_{0.7}$ from p. 174) from the true margins
q07 <- c(
  qweibull(0.7, shape = 1.5, scale = 4),
  qweibull(0.7, shape = 2, scale = 3),
  qgamma(0.7, shape = 2, rate = 1),
  qgamma(0.7, shape = 3, rate = 1)
)

# True Joint Probability P(Y > q07) = P(Y1 > q07[1], ..., Y4 > q07[4]) = P(U1 > 0.7, ..., U4 > 0.7) for the true model
true_P_surv <- prob(true_copula, l = c(0.7, 0.7, 0.7, 0.7), u = c(1, 1, 1, 1))

###### CANDIDATE MODELS GRID
candidates_cop <- c("gumbel", "normal", "frank", "surv_clayton")
candidates_margin1 <- c("weibull", "gamma", "lnorm")
candidates_margin2 <- c("weibull", "gamma", "norm")
candidates_margin3 <- c("gamma", "weibull", "lnorm")
candidates_margin4 <- c("gamma", "weibull", "lnorm")

# Create a grid of all 324 combinations
model_grid <- expand.grid(
  Copula = candidates_cop, margin1 = candidates_margin1, margin2 = candidates_margin2, 
  margin3 = candidates_margin3, margin4 = candidates_margin4, stringsAsFactors = FALSE
)

###### HELPER FUNCTIONS

# Calculate MLE of the marginal distributions and then fit the estimated model to the data
fit_margin <- function(data_vec, dist_name) {
  n <- length(data_vec)
  # For Weibull and Gamma distributions, lower=c(1e-8, 1e-8) ensures both shape and scale are positive
  if(dist_name == "weibull") {
    nll <- function(p) -sum(dweibull(data_vec, p[1], p[2], log=TRUE))
    fit <- optim(c(1,1), nll, lower=c(1e-8, 1e-8), method="L-BFGS-B")
    u_hat <- pweibull(data_vec, fit$par[1], fit$par[2])
  } else if(dist_name == "gamma") {
    nll <- function(p) -sum(dgamma(data_vec, p[1], p[2], log=TRUE))
    fit <- optim(c(1,1), nll, lower=c(1e-8, 1e-8), method="L-BFGS-B")
    u_hat <- pgamma(data_vec, fit$par[1], fit$par[2])
  } else if(dist_name == "lnorm") {
    nll <- function(p) -sum(dlnorm(data_vec, p[1], p[2], log=TRUE))
    # First-stage estimation is data-dependent for the log-normal/normal distribution
    # Also note that since the true margins are Weibull/Gamma, log(data_vec) is always well-defined
    fit <- optim(c(mean(log(data_vec)), sd(log(data_vec))), nll, lower=c(-Inf, 1e-8), method="L-BFGS-B")
    u_hat <- plnorm(data_vec, fit$par[1], fit$par[2])
  } else if(dist_name == "norm") {
    nll <- function(p) -sum(dnorm(data_vec, p[1], p[2], log=TRUE))
    fit <- optim(c(mean(data_vec), sd(data_vec)), nll, lower=c(-Inf, 1e-8), method="L-BFGS-B")
    u_hat <- pnorm(data_vec, fit$par[1], fit$par[2])
  }
  return(list(par = fit$par, nll = fit$value, u_hat = u_hat))
}

# Calculating margin penalty $\widetilde{p}_{\alpha_{j}}^{*}$ based on p.172 of (1)

calc_margin_penalty <- function(data_vec, dist_name, param_hat) {
  n <- length(data_vec)
  p1 <- param_hat[1]
  p2 <- param_hat[2]
  
  # For the Weibull and Gamma distributions, both the shape and scale must be positive
  if(dist_name %in% c("weibull", "gamma")) {
    p1 <- max(p1, 1e-8); p2 <- max(p2, 1e-8)
    # For the Normal and Log-Normal distributions, standard deviation must be positive
  } else if(dist_name %in% c("lnorm", "norm")) {
    p2 <- max(p2, 1e-8)
  }
  
  # Initialize Score Matrix (U_alpha) and Hessian Sum (H_sum)
  # U_alpha is an n x 2 matrix (scores vector for each parameter)
  # H_sum is a 2 x 2 matrix (sum of the second derivatives)
  U_alpha <- matrix(0, nrow=n, ncol=2)
  H_sum <- matrix(0, nrow=2, ncol=2)
  
  if (dist_name == "norm") {
    mu <- p1; sig <- p2
    
    # Vectorized Analytical Gradient (Score)
    U_alpha[,1] <- (data_vec - mu) / sig^2
    U_alpha[,2] <- -1/sig + (data_vec - mu)^2 / sig^3
    
    # Analytical Hessian Sum
    H11 <- sum(rep(-1/sig^2, n))
    H22 <- sum(1/sig^2 - 3*(data_vec - mu)^2 / sig^4)
    H12 <- sum(-2*(data_vec - mu) / sig^3)
    H_sum <- matrix(c(H11, H12, H12, H22), nrow=2, ncol=2)
    
  } else if (dist_name == "lnorm") {
    mu <- p1; sig <- p2
    lx <- log(data_vec)
    
    # Vectorized Analytical Gradient (Score)
    U_alpha[,1] <- (lx - mu) / sig^2
    U_alpha[,2] <- -1/sig + (lx - mu)^2 / sig^3
    
    # Analytical Hessian Sum
    H11 <- sum(rep(-1/sig^2, n))
    H22 <- sum(1/sig^2 - 3*(lx - mu)^2 / sig^4)
    H12 <- sum(-2*(lx - mu) / sig^3)
    H_sum <- matrix(c(H11, H12, H12, H22), nrow=2, ncol=2)
    
  } else if (dist_name == "weibull") {
    k <- p1; lambda <- p2
    
    # Vectorized Analytical Gradient (Score)
    U_alpha[,1] <- -k/lambda + (k/lambda) * ((data_vec/lambda)^k)
    U_alpha[,2] <- 1/k - log(data_vec/lambda)*((data_vec/lambda)^k-1)
    
    # Analytical Hessian Sum
    H11 <- -sum(k*(k*(data_vec/lambda)^k + (data_vec/lambda)^k - 1)/lambda^2)
    H22 <- sum(-((data_vec/lambda)^k)*(log(data_vec/lambda))^2 - 1/k^2)
    H12 <- sum(((data_vec/lambda)^k + k*((data_vec/lambda)^k)*log(data_vec/lambda) - 1)/lambda)
    H_sum <- matrix(c(H11, H12, H12, H22), nrow=2, ncol=2)
    
  } else if (dist_name == "gamma") {
    # Using 'scale' parameterization 
    alpha <- p1; beta <- p2
    
    # Vectorized Analytical Gradient (Score)
    # R has built-in functions for the gamma derivatives: digamma() and trigamma()
    U_alpha[,1] <- -log(beta) - digamma(alpha) + log(data_vec)
    U_alpha[,2] <- -alpha/beta + data_vec/beta^2
    
    # Analytical Hessian Sum
    H11 <- sum(rep(-trigamma(alpha), n))
    H22 <- sum(alpha/beta^2 - 2*data_vec/beta^3)
    H12 <- sum(rep(-1/beta, n))
    H_sum <- matrix(c(H11, H12, H12, H22), nrow=2, ncol=2)
  }
  
  # Calculate K and I Matrices using the exact analytical results
  # K_alpha is the empirical covariance of the scores
  K_alpha <- t(U_alpha) %*% U_alpha / n
  
  # I_alpha is the negative expected Hessian
  I_alpha <- -(H_sum / n)
  
  # Calculate trace penalty
  tr_val <- sum(diag(solve(I_alpha + diag(1e-8, nrow(I_alpha))) %*% K_alpha))
  
  # THE FIX: Export the matrices so we can use them for the copula penalty later!
  return(list(pen = tr_val, U_alpha = U_alpha, I_alpha = I_alpha))
}

# Calculating copula penalty $\widetilde{p}_{\theta}}^{*}$ based on p.172 of (1)

calc_copula_penalty <- function(margin1_data, margin2_data, margin3_data, margin4_data, model_df, fit_cop, Y_true){
  # Combine the cached marginal matrices
  U_alpha_full <- cbind(margin1_data$U_alpha, margin2_data$U_alpha, margin3_data$U_alpha, margin4_data$U_alpha)
  
  I_alpha_full <- matrix(0, 8, 8)
  I_alpha_full[1:2, 1:2] <- margin1_data$I_alpha
  I_alpha_full[3:4, 3:4] <- margin2_data$I_alpha
  I_alpha_full[5:6, 5:6] <- margin3_data$I_alpha
  I_alpha_full[7:8, 7:8] <- margin4_data$I_alpha
  
  # Define the joint copula log-likelihood function for numDeriv
  cop_log_lik_vec <- function(p) {
    u1 <- get_u_hat(Y_true[, 1], model_def$margin1, p[1:2])
    u2 <- get_u_hat(Y_true[, 2], model_def$margin2, p[3:4])
    u3 <- get_u_hat(Y_true[, 3], model_def$margin3, p[5:6])
    u4 <- get_u_hat(Y_true[, 4], model_def$margin4, p[7:8])
    
    U_mat <- pmin(pmax(cbind(u1, u2, u3, u4), 1e-5), 1 - 1e-5)
    
    # Clamp copula boundaries so numDeriv doesn't crash
    th <- p[9]
    if (model_def$Copula == "gumbel") th <- max(th, 1.01)
    else if (model_def$Copula == "surv_clayton") th <- max(th, 0.01)
    else if (model_def$Copula == "normal") th <- pmin(pmax(th, -0.99), 0.99)
    else if (model_def$Copula == "frank" && abs(th) < 1e-4) th <- 1e-4
    
    cop <- switch(model_def$Copula,
                  "gumbel" = gumbelCopula(th, dim=4),
                  "normal" = normalCopula(th, dim=4),
                  "frank" = frankCopula(th, dim=4),
                  "surv_clayton" = rotCopula(claytonCopula(th, dim=4)))
    
    return(log(pmax(dCopula(U_mat, cop), 1e-10)))
  }
  
  all_params <- c(margin1_data$par, margin2_data$par, margin3_data$par, margin4_data$par, fit_cop@estimate)
  n_obs_loop <- nrow(Y_true)
  
  # Calculate Copula Scores (Jacobian) and Information (Hessian)
  J_cop <- try(numDeriv::jacobian(func = cop_log_lik_vec, x = all_params), silent = TRUE)
  H_cop <- try(numDeriv::hessian(func = function(p) sum(cop_log_lik_vec(p)), x = all_params), silent = TRUE)
  
  if (inherits(J_cop, "try-error") || inherits(H_cop, "try-error") || 
      any(is.na(J_cop)) || any(is.na(H_cop))) return(NA)
  
  # Extract Sub-matrices
  U_alpha_star <- J_cop[, 1:8]
  U_theta <- J_cop[, 9, drop = FALSE]
  
  I_cop <- -(H_cop / n_obs_loop)
  I_alpha_theta <- I_cop[1:8, 9, drop = FALSE]
  I_theta <- I_cop[9, 9, drop = FALSE]
  
  # Empirical Covariances (K matrices)
  K_alpha_circ <- t(U_alpha_full) %*% U_alpha_star / n_obs_loop
  K_alpha_theta <- t(U_alpha_full) %*% U_theta / n_obs_loop
  K_theta <- t(U_theta) %*% U_theta / n_obs_loop
  
  # Calculate the exact Trace Formula, p.172 of (1)
  calc_traces <- try({
    I_alpha_inv <- solve(I_alpha_full + diag(1e-8, 8))
    I_theta_inv <- solve(I_theta + diag(1e-8, 1))
    
    t1 <- sum(diag(I_alpha_inv %*% K_alpha_circ))
    t2 <- as.numeric(I_theta_inv %*% t(I_alpha_theta) %*% I_alpha_inv %*% K_alpha_theta)
    t3 <- as.numeric(I_theta_inv %*% K_theta)
    t1 - t2 + t3
  }, silent = TRUE)
  
  if (inherits(calc_traces, "try-error") || is.na(calc_traces)) return(NA)
  
  p_theta_star <- calc_traces
  return(p_theta_star)
}

# CDF of the estimated margins
get_u_hat <- function(q, dist, p) {
  if (dist == "weibull") pweibull(q, p[1], p[2])
  else if (dist == "gamma") pgamma(q, p[1], p[2])
  else if (dist == "lnorm") plnorm(q, p[1], p[2])
  else if (dist == "norm") pnorm(q, p[1], p[2])
}

# For each possible marginal model, compute MLE and its margin penalty
precompute_margin <- function(y_vec, candidates, q_val) {
  res <- list()
  for(dist in unique(candidates)) {
    # MLE for the candidate model based on the data y_vec
    fit <- try(fit_margin(y_vec, dist), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      margin_obj <- try(calc_margin_penalty(y_vec, dist, fit$par), silent = TRUE)
      u_hat_P_surv <- get_u_hat(q_val, dist, fit$par)
      # Store the U_alpha and I_alpha matrices in the cache
      res[[dist]] <- list(par = fit$par, nll = fit$nll, u_hat = fit$u_hat, 
                            pen = margin_obj$pen, U_alpha = margin_obj$U_alpha, 
                          I_alpha = margin_obj$I_alpha, u_hat_P_surv = u_hat_P_surv)
    } else {
      # Keep the EXACT same list structure, filling missing values with NA
      res[[dist]] <- list(par = NA, nll = NA, u_hat = NA, pen = NA, 
                          U_alpha = NA, I_alpha = NA, u_hat_P_surv = NA)
    }
  }
  return(res)
}

###### MAIN SIMULATION LOOP (PARALLELIZED & CACHED)
# Safely grab the number of cores assigned by the Slurm scheduler
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores)) n_cores <- 1 # Fallback for local testing
cat(sprintf("Setting up cluster with %d cores...\n", n_cores))

# Set up an output file for the parallel workers to write logs to
log_file <- "simulation_progress.log"
file.create(log_file) # Initialize empty log file

cl <- makeCluster(n_cores, outfile = "") # outfile="" allows workers to print to standard output
registerDoParallel(cl)

results_list <- foreach(
  sim = 1:n_simulations,
  .packages = c("copula", "numDeriv", "dplyr")
) %dopar% {
  # Write progress to the log file
  cat(sprintf("Starting Simulation %d of %d at %s\n", sim, n_simulations, Sys.time()), 
      file = log_file, append = TRUE)
  # Generate true data
  U_true <- rCopula(n_obs, true_copula)
  Y_true <- cbind(
    qweibull(U_true[, 1], 1.5, 4),
    qweibull(U_true[, 2], 2, 3),
    qgamma(U_true[, 3], 2, 1),
    qgamma(U_true[, 4], 3, 1)
  )
  
  sim_results <- data.frame(
    SimID = rep(sim, nrow(model_grid)),
    ModelID = 1:nrow(model_grid),
    CIC = NA,
    AIC2ML = NA,
    P_est = NA,
    pen_CIC = NA,
    pen_AIC = NA
  )
  
  # THE CACHE: Pre-compute margins once per simulation
  cache_margin1 <- precompute_margin(Y_true[, 1], candidates_margin1, q07[1])
  cache_margin2 <- precompute_margin(Y_true[, 2], candidates_margin2, q07[2])
  cache_margin3 <- precompute_margin(Y_true[, 3], candidates_margin3, q07[3])
  cache_margin4 <- precompute_margin(Y_true[, 4], candidates_margin4, q07[4])
  
  # Inner loop: evaluate all 324 candidate models using cached values
  for (m_idx in 1:nrow(model_grid)) {
    model_def <- model_grid[m_idx, ]
    
    # Fast Lookup from Cache
    margin1_data <- cache_margin1[[model_def$margin1]]
    margin2_data <- cache_margin2[[model_def$margin2]]
    margin3_data <- cache_margin3[[model_def$margin3]]
    margin4_data <- cache_margin4[[model_def$margin4]]
    
    # Check if any of the marginal fits failed (indicated by NA in the nll)
    if (is.na(margin1_data$nll) || is.na(margin2_data$nll) || 
        is.na(margin3_data$nll) || is.na(margin4_data$nll)) {
      next
    }
    
    u_hat <- cbind(margin1_data$u_hat, margin2_data$u_hat, margin3_data$u_hat, margin4_data$u_hat)
    
    # GUARDRAIL: Clamp u_hat
    u_hat <- pmin(pmax(u_hat, 1e-8), 1 - 1e-8)
    
    # Fit Copula
    cop_obj <- switch(
      model_def$Copula,
      "gumbel" = gumbelCopula(dim = 4),
      "normal" = normalCopula(dim = 4),
      "frank" = frankCopula(dim = 4),
      "surv_clayton" = rotCopula(claytonCopula(dim = 4))
    )
    
    fit_cop <- try(fitCopula(cop_obj, u_hat, method = "ml"), silent = TRUE)
    if (inherits(fit_cop, "try-error")) next
    
    # Log-Likelihoods & AIC
    logL_margins <- -(margin1_data$nll + margin2_data$nll + margin3_data$nll + margin4_data$nll) # minus sign: negative log likelihood
    total_logL <- logL_margins + as.numeric(logLik(fit_cop))
    
    # Dynamically calculate parameter count
    p_eta <- 8 + length(fit_cop@estimate)
    sim_results$AIC2ML[m_idx] <- -2 * total_logL + 2 * p_eta
    sim_results$pen_AIC[m_idx] <- p_eta # Save it!
    
    if (inherits(margin1_data$pen, "try-error") || inherits(margin2_data$pen, "try-error") ||
        inherits(margin3_data$pen, "try-error") || inherits(margin4_data$pen, "try-error")) next
    
    # EXACT COPULA PENALTY (p_theta_star)
    p_theta_star <- calc_copula_penalty(margin1_data, margin2_data, margin3_data, margin4_data, model_df, fit_cop, Y_true)
    
    # Catch the NA returned by the function if numDeriv crashed
    if (is.na(p_theta_star)) next
   
    # Sum everything for final CIC
    p_eta_star <- margin1_data$pen + margin2_data$pen + margin3_data$pen + margin4_data$pen + p_theta_star
    sim_results$CIC[m_idx] <- -2 * total_logL + 2 * p_eta_star
    sim_results$pen_CIC[m_idx] <- p_eta_star
    
    # Evaluate P (Quantiles Looked up)
    u_hat_P_surv_all_margins <- c(margin1_data$u_hat_P_surv, margin2_data$u_hat_P_surv, margin3_data$u_hat_P_surv, margin4_data$u_hat_P_surv)
    sim_results$P_est[m_idx] <- prob(fit_cop@copula, l = u_hat_P_surv_all_margins, u = c(1, 1, 1, 1))
  }
  
  return(sim_results)
}

stopCluster(cl)
cat("Simulations complete!\n")

all_sims <- do.call(rbind, results_list)

# Update your Aggregation to include the penalties
agg_results <- all_sims %>%
  group_by(ModelID) %>%
  summarize(
    Mean_CIC = mean(CIC, na.rm = TRUE),
    Mean_AIC2ML = mean(AIC2ML, na.rm = TRUE),
    MSE_P = mean((P_est - true_P_surv)^2, na.rm = TRUE),
    Mean_pen_CIC = mean(pen_CIC, na.rm = TRUE), # Average CIC penalty
    Mean_pen_AIC = mean(pen_AIC, na.rm = TRUE)  # Average AIC penalty
  ) %>%
  filter(!is.na(Mean_CIC))

agg_results$Scaled_MSE <- agg_results$MSE_P * 1000

# Save the aggregated results to an RDS file
saveRDS(agg_results, file = "Simulation1_Results.rds")

# Also save the raw, unaggregated data just in case
saveRDS(all_sims, file = "Simulation1_Raw_Data.rds")

cat("Results successfully saved to RDS files.\n")


