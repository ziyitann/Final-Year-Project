# FYP Final R Code
# Zi Tan - 220016955

# Check whether packages are installed
if (!require("ChainLadder")) install.packages("ChainLadder")
if (!require("DCL")) install.packages("DCL")
library(ChainLadder)
library(DCL)

# Baseline simulation settings
set.seed(123)
m_obs <- 5L # only the first 5 development years are observed
m_true <- 10L # true development runs for 10 years
AY <- 5L # number of accident years

# Severity distribution setup
mu <- 208.37
phi <- 10069.1
gamma_scale <- phi / mu
gamma_shape <- mu / gamma_scale

case_factor_obs <- c(1.05, 0.85, 0.60, 0.40, 0.25)  # case reserve factor

beta_fast <- diff(c(0, 0.60, 0.76, 0.83, 0.89, 0.92, 0.96, 0.985, 0.992, 0.996, 1.00))  # fast reporting
beta_slow <- diff(c(0, 0.20, 0.45, 0.65, 0.80, 0.90, 0.95, 0.97, 0.985, 0.995, 1.00))  # slow reporting

delta_long  <- c(0.10, 0.24, 0.23, 0.16, 0.10, 0.07, 0.04, 0.03, 0.02, 0.01)  # long settlement tail
delta_short <- c(0.60, 0.25, 0.10, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)  # short settlement tail

# Keep only the observed upper triangle
mask_triangle <- function(mat, m_obs) {
  mat[row(mat) + col(mat) > (m_obs + 1)] <- NA
  mat
}

# Simulate full reported and paid triangles
simulate_portfolio <- function(alpha, beta_pattern, delta_pattern,
                               m_true, gamma_shape, gamma_scale) {
  N_reported_full <- matrix(0, nrow = length(alpha), ncol = m_true) # matrix for reported claim counts
  X_paid_full <- matrix(0, nrow = length(alpha), ncol = m_true) # matrix for paid claim amounts
  N_paid_full <- matrix(0, nrow = length(alpha), ncol = m_true) # matrix for paid claim counts
  
  for (i in seq_along(alpha)) {
    total_claims <- rpois(1, lambda = alpha[i]) # simulate total number of claims in accident year i
    if (total_claims == 0) next
    
    report_dev <- sample(1:m_true, size = total_claims,
                         replace = TRUE, prob = beta_pattern) # simulate reporting year for each claim
    
    for (k in seq_len(total_claims)) {
      r <- report_dev[k] # reporting development year for this claim
      if (r > m_true) next
      
      N_reported_full[i, r] <- N_reported_full[i, r] + 1
      
      max_delay <- m_true - r
      delay_probs <- delta_pattern[1:(max_delay + 1)]
      delay_probs <- delay_probs / sum(delay_probs)
      
      d <- sample(0:max_delay, size = 1, prob = delay_probs) # simulate payment delay after reporting
      p <- r + d
      
      X_paid_full[i, p] <- X_paid_full[i, p] +
        rgamma(1, shape = gamma_shape, scale = gamma_scale) # simulate payment amount
      N_paid_full[i, p] <- N_paid_full[i, p] + 1  # add one paid claim count
    }
  }
  
  list(
    N_reported_full = N_reported_full,
    X_paid_full = X_paid_full,
    N_paid_full = N_paid_full
  )
}

# Build incremental incurred triangle
build_incurred_triangle <- function(cum_X_pay, open_counts,
                                    case_factor_obs, mu, bias_vec = NULL) {
  case_reserves <- sweep(open_counts, 2, case_factor_obs * mu, "*") # open claims × severity mean × reserve factor
  
  if (!is.null(bias_vec)) {
    bias_mat <- matrix(bias_vec, nrow = nrow(open_counts),
                       ncol = ncol(open_counts), byrow = FALSE)
    case_reserves <- case_reserves * bias_mat
  }
  
  I_obs_cum <- t(apply(cum_X_pay + case_reserves, 1, cummax))
  I_obs_inc <- I_obs_cum
  I_obs_inc[, 2:ncol(I_obs_cum)] <-
    I_obs_cum[, 2:ncol(I_obs_cum)] - I_obs_cum[, 1:(ncol(I_obs_cum) - 1)]
  
  I_obs_inc
}

# -------------------------------------------------------------
# Part 1: Example triangles used in Methodology (Tables 2 to 4)
# -------------------------------------------------------------

# Starting claim count for the 5 accident years, with 40% annual growth
alpha_example <- 50 * (1.40)^(0:(AY - 1))

# Simulate the full portfolio using the baseline scenario assumptions
sim_example <- simulate_portfolio(
  alpha = alpha_example,
  beta_pattern = beta_fast,
  delta_pattern = delta_long,
  m_true = m_true,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale
)

# Keep only the first 5 observed development years
N_obs_inc <- sim_example$N_reported_full[, 1:m_obs]
X_obs_inc <- sim_example$X_paid_full[, 1:m_obs]

cum_N_rep <- t(apply(N_obs_inc, 1, cumsum))
cum_N_pay <- t(apply(sim_example$N_paid_full[, 1:m_obs], 1, cumsum))
cum_X_pay <- t(apply(X_obs_inc, 1, cumsum))

open_counts <- pmax(cum_N_rep - cum_N_pay, 0) # open claims = cumulative reported counts - cumulative paid counts

# Build the incremental incurred triangle
I_obs_inc <- build_incurred_triangle(
  cum_X_pay = cum_X_pay,
  open_counts = open_counts,
  case_factor_obs = case_factor_obs,
  mu = mu
)

# Mask each matrix so only the observed upper triangle remains
N_tri <- mask_triangle(N_obs_inc, m_obs)
X_tri <- mask_triangle(X_obs_inc, m_obs)
I_tri <- mask_triangle(I_obs_inc, m_obs)

# Print the three example triangles
cat("\nTable 2: Reported claim counts triangle (N_ij):\n")
print(N_tri)

cat("\nTable 3: Incremental paid triangle (X_ij):\n")
print(round(X_tri, 0))

cat("\nTable 4: Incremental incurred triangle (I_ij):\n")
print(round(I_tri, 0))

# ----------------------------------------------------------------
# Part 2: Main scenario results used in the Results chapter
# ----------------------------------------------------------------

run_scenario <- function(scenario_name, alpha_base, beta_pattern, delta_pattern,
                         result_table_no) {
  # Reset seed
  set.seed(123)
  
  alpha <- alpha_base * (1.40)^(0:(AY - 1))
  
  # Simulate the full portfolio under this scenario
  sim <- simulate_portfolio(
    alpha = alpha,
    beta_pattern = beta_pattern,
    delta_pattern = delta_pattern,
    m_true = m_true,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale
  )
  
  # Keep only the observed part
  N_obs_inc <- sim$N_reported_full[, 1:m_obs]
  X_obs_inc <- sim$X_paid_full[, 1:m_obs]
  
  cum_N_rep <- t(apply(N_obs_inc, 1, cumsum))
  cum_N_pay <- t(apply(sim$N_paid_full[, 1:m_obs], 1, cumsum))
  cum_X_pay <- t(apply(X_obs_inc, 1, cumsum))
  
  open_counts <- pmax(cum_N_rep - cum_N_pay, 0)
  
  I_obs_inc <- build_incurred_triangle(
    cum_X_pay = cum_X_pay,
    open_counts = open_counts,
    case_factor_obs = case_factor_obs,
    mu = mu
  )
  
  N_tri <- mask_triangle(N_obs_inc, m_obs)
  X_tri <- mask_triangle(X_obs_inc, m_obs)
  I_tri <- mask_triangle(I_obs_inc, m_obs)
  
  true_ult <- rowSums(sim$X_paid_full) # true ultimate = full row sums of paid amounts
  paid_to_date <- rowSums(X_tri, na.rm = TRUE) # paid to date = row sums of observed paid triangle
  
  # CLM
  cl_fit <- MackChainLadder(as.triangle(incr2cum(X_tri)), tail = FALSE)
  cl_ult <- cl_fit$FullTriangle[, m_obs]
  
  # DCL and BDCL
  suppressWarnings({
    dcl_fit <- dcl.estimation(X_tri, N_tri, adj = 1, Tables = FALSE)
    dcl_pred <- dcl.predict(dcl_fit, Ntriangle = N_tri, Tail = TRUE,
                            summ.by = "row", Tables = FALSE)
    dcl_ult <- paid_to_date + rowSums(dcl_pred$Xtotal, na.rm = TRUE)
    
    bdcl_fit <- bdcl.estimation(X_tri, N_tri, I_tri, adj = 1, Tables = FALSE)
    bdcl_pred <- dcl.predict(bdcl_fit, Ntriangle = N_tri, Tail = TRUE,
                             summ.by = "row", Tables = FALSE)
    bdcl_ult <- paid_to_date + rowSums(bdcl_pred$Xtotal, na.rm = TRUE)
  })
  
  # Table showing projected ultimates and percentage forecast errors
  results <- data.frame(
    AY = 1:AY,
    True_Ult = round(true_ult, 0),
    CL = round(cl_ult, 0),
    DCL = round(dcl_ult, 0),
    BDCL = round(bdcl_ult, 0),
    CL_Error = round(100 * (cl_ult - true_ult) / true_ult, 1),
    DCL_Error = round(100 * (dcl_ult - true_ult) / true_ult, 1),
    BDCL_Error = round(100 * (bdcl_ult - true_ult) / true_ult, 1),
    check.names = FALSE
  )
  
  cat("\n--------------------------------------------------------\n")
  cat("SCENARIO:", scenario_name, "\n")
  cat("--------------------------------------------------------\n")
  
  cat(sprintf("\nTable %d: Projected ultimate losses and percentage forecast error:\n",
              result_table_no))
  print(results)
  
  invisible(list(
    results = results
  ))
}

run_scenario(
  scenario_name = "Scenario 1: Fast Reporting with Long-Tail Settlement",
  alpha_base = 50,
  beta_pattern = beta_fast,
  delta_pattern = delta_long,
  result_table_no = 5
)

run_scenario(
  scenario_name = "Scenario 2: Short-Tail Development",
  alpha_base = 50,
  beta_pattern = beta_fast,
  delta_pattern = delta_short,
  result_table_no = 6
)

run_scenario(
  scenario_name = "Scenario 3: Delayed Reporting and Sparse Data",
  alpha_base = 35,
  beta_pattern = beta_slow,
  delta_pattern = delta_long,
  result_table_no = 7
)

# ----------------------------------------------
# Part 3: Scenario 3 sensitivity check (Table 8)
# ----------------------------------------------

run_scenario3_sensitivity <- function() {
  set.seed(123)
  
  alpha_base <- 35
  alpha <- alpha_base * (1.40)^(0:(AY - 1))
  
  # Bias is concentrated in the more recent accident years
  bias_recent_vec <- c(1.00, 1.00, 0.95, 0.80, 0.60)
  
  # Simulate the Scenario 3 portfolio
  sim <- simulate_portfolio(
    alpha = alpha,
    beta_pattern = beta_slow,
    delta_pattern = delta_long,
    m_true = m_true,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale
  )
  
  N_obs_inc <- sim$N_reported_full[, 1:m_obs]
  X_obs_inc <- sim$X_paid_full[, 1:m_obs]
  
  cum_N_rep <- t(apply(N_obs_inc, 1, cumsum))
  cum_N_pay <- t(apply(sim$N_paid_full[, 1:m_obs], 1, cumsum))
  cum_X_pay <- t(apply(X_obs_inc, 1, cumsum))
  
  open_counts <- pmax(cum_N_rep - cum_N_pay, 0)
  
  N_tri <- mask_triangle(N_obs_inc, m_obs)
  X_tri <- mask_triangle(X_obs_inc, m_obs)
  
  true_ult <- rowSums(sim$X_paid_full)
  paid_to_date <- rowSums(X_tri, na.rm = TRUE)
  
  # DCL
  suppressWarnings({
    dcl_fit <- dcl.estimation(X_tri, N_tri, adj = 1, Tables = FALSE)
    dcl_pred <- dcl.predict(dcl_fit, Ntriangle = N_tri, Tail = TRUE,
                            summ.by = "row", Tables = FALSE)
    dcl_ult <- paid_to_date + rowSums(dcl_pred$Xtotal, na.rm = TRUE)
  })
  dcl_err <- 100 * (dcl_ult - true_ult) / true_ult  # percentage forecast error for DCL
  
  # BDCL with known mean severity
  I_obs_inc_known <- build_incurred_triangle(
    cum_X_pay = cum_X_pay,
    open_counts = open_counts,
    case_factor_obs = case_factor_obs,
    mu = mu
  )
  I_tri_known <- mask_triangle(I_obs_inc_known, m_obs)
  
  suppressWarnings({
    bdcl_fit_known <- bdcl.estimation(X_tri, N_tri, I_tri_known, adj = 1, Tables = FALSE)
    bdcl_pred_known <- dcl.predict(bdcl_fit_known, Ntriangle = N_tri, Tail = TRUE,
                                   summ.by = "row", Tables = FALSE)
    bdcl_ult_known <- paid_to_date + rowSums(bdcl_pred_known$Xtotal, na.rm = TRUE)
  })
  bdcl_err_known <- 100 * (bdcl_ult_known - true_ult) / true_ult  # percentage forecast error for BDCL with known mean severity
  
  # BDCL with biased incurred data
  I_obs_inc_biased <- build_incurred_triangle(
    cum_X_pay = cum_X_pay,
    open_counts = open_counts,
    case_factor_obs = case_factor_obs,
    mu = mu,
    bias_vec = bias_recent_vec
  )
  I_tri_biased <- mask_triangle(I_obs_inc_biased, m_obs)
  
  suppressWarnings({
    bdcl_fit_biased <- bdcl.estimation(X_tri, N_tri, I_tri_biased, adj = 1, Tables = FALSE)
    bdcl_pred_biased <- dcl.predict(bdcl_fit_biased, Ntriangle = N_tri, Tail = TRUE,
                                    summ.by = "row", Tables = FALSE)
    bdcl_ult_biased <- paid_to_date + rowSums(bdcl_pred_biased$Xtotal, na.rm = TRUE)
  })
  bdcl_err_biased <- 100 * (bdcl_ult_biased - true_ult) / true_ult # percentage forecast error for BDCL with biased incurred data
  
  sensitivity_table <- data.frame(
    Method = c(
      "True Ultimate",
      "DCL",
      "BDCL with known mean severity",
      "BDCL with biased incurred data"
    ),
    Ultimate = round(c(
      true_ult[5],
      dcl_ult[5],
      bdcl_ult_known[5],
      bdcl_ult_biased[5]
    ), 0),
    `Percent Error` = c(
      NA,
      round(dcl_err[5], 1),
      round(bdcl_err_known[5], 1),
      round(bdcl_err_biased[5], 1)
    ),
    check.names = FALSE
  )
  
  cat("\n--------------------------------------------------------\n")
  cat("TABLE 8: SCENARIO 3 ALTERNATIVE COMPARISON FOR AY 5\n")
  cat("--------------------------------------------------------\n")
  print(sensitivity_table, row.names = FALSE)
  
  invisible(list(
    true_ult_ay5 = true_ult[5],
    sensitivity_table = sensitivity_table,
    dcl_ult = dcl_ult,
    bdcl_ult_known = bdcl_ult_known,
    bdcl_ult_biased = bdcl_ult_biased,
    dcl_err = dcl_err,
    bdcl_err_known = bdcl_err_known,
    bdcl_err_biased = bdcl_err_biased
  ))
}

# Run the Scenario 3 alternative comparison
scenario3_sensitivity <- run_scenario3_sensitivity()
