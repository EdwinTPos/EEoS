########################
# Code for
# A coherent signature of disturbance in successional forests
# Code written and maintained by Edwin Pos and Erica Newman
# E.T.Pos@uu.nl, ericanewman.sci@gmail.com
# Updated July 13, 2026
########################

# These functions are needed to express the sum over R and nR analytically.

# METE unnormalized structure function: n,e are microscopic variables, l are lambdas and remainder are state variables, call S, N, or E
R <- function(n, e, l, s) {
return(exp(-l[1]*n-l[2]*n*e))}

# METE unnormalized structure function summed over n: n,e are microscopic variables, l are lambdas and remainder are state variables, call S, N, or E. Define exponent for lambdas 1-2.
Rsum <- function(e, l, s) {
  l12 <- l[1] + l[2] * e
  num <- exp(-l12) - exp(-l12 * (s$N + 1))
  denom <- 1 - exp(-l12)
  return(num / denom)}

# METE unnormalized structure function multiplied by n, then summed over n: e is a microscopic variable, l are lambdas and remainder are state variables (S, N, or E). Define exponent for lambdas 1-2, with n=1 since we've done sum over n already.
nRsum <- function(e, l, s) {
  l12 <- l[1] + l[2] * e
  num <- exp(-l12) - (s$N + 1) * exp(-l12 * (s$N + 1)) + s$N * exp(-l12 * (s$N + 2))
  denom <- (1 - exp(-l12))^2
  return(num / denom)}

# Calculates partition function for DynaMETE, uses the analytic sum from Rsum and single integral over e. Integral is done with quad over log e, returns only value not error.
zfun <- function(l, s) {
integrand = 
function(loge){exp(loge)*Rsum(exp(loge), l, s)}
return(integrate(integrand,0,log(s$E)))}

# Constraint function takes in lambdas and state variables and tests the constraints numerically. The return is an array with the percent deviation from the exact constraint. First get the normalization, then take  the sum over n and ne, normalize, get exact values and return percent difference:

constraints <- function(l, s) {
z <- zfun(l, s)$value
nR <- integrate(function(loge){exp(loge)* nRsum(exp(loge), l, s)}, 0, log(s$E))$value
neR <- integrate(function(loge) {exp(loge * 2) * nRsum(exp(loge), l, s)}, 0, log(s$E))$value
  nR <- nR / z
  neR <- neR / z
  ns <- s$N / s$S
  es <- s$E / s$S
  return(c((nR - ns) / ns, (neR - es) / es))
}


#Function to call and calculate beta_constraints in METE with given state variables. 
#Inputs s as state variables, call S, N, or E. Outputs beta constraint to minimize.

beta_constraint <- function(x, s) {
narr <- seq(1, s$N)
a <- mpfr(exp(-100), precBits=64)
vec = (narr * x)/100
vec_mpfr = a^vec
return(asNumeric(sum(a^vec) / sum(a^vec/narr) - s$N / s$S))}

#Mete_lambdas returns the METE lambdas for a given set of state variables. 
#Inputs are state variables (S,N or E), outputs are arrays of lambdas necessary for the other functions. 
#It also checks to see if conditions are actually satisfied, and outputs a warning if the percent difference of the result is greater than 5% (by default). 
#This can be modified to exclude this test by setting thresh to False, or changing it to a difference percent

mete_lambdas <- function(s, thresh=0.05) {
  beta <- as.numeric(uniroot(beta_constraint, interval = c(1e-05,1e-03), lower=1e-05,upper =1 ,s=s)[1])
  l2 <- s$S / (s$E - s$N)
  ls <- c(beta - l2, l2)
  con <- constraints(ls, s)
    if (any(con > thresh)) {
      print(paste("Constraints are not satisfied within", thresh*100, "%"))
    } 
  return(ls) }

biomass_approx <- function(s, order =0){

# The approximation for the biomass equation with power=4/3 as from the original publication by Harte et al.
# order indicates whether or not to include a first order correction. The default is to not, thus order=0.
# Change order=1 if you want to include the correction.
l <- mete_lambdas(s)
b <- 4.17*s$E^(4/3)/(s$S^(1/3)*log(1/l[1]+l[2]))
if (order == 1) {
b <- b * (1-1.16*(l[1]+l[2])^(1/3))
}
return(b)
}

############################################################
### FUNCTION TO COMPARE IPD AND RAD FOR METE PREDICTION ####
############################################################

# Samples specifies the unique identifiers for each location, metric is a choicen of comparison quantification 
# and can be either "AIC", "LL" (loglikelihood), "LLZ" (z-score transformed loglikelihood) or "KL" (Kullback-Leibler)

METE_comparison <- function(Samples, metric,data){

# Create output table for all calculations of METE and the quantitative metrics for comparisons
output_table = data.frame(Samples)
split_code_year <- t(data.frame(lapply(output_table[,1], function(x) as.vector(t(str_match(x, '(.*)_(.*)') [,-1])))))
output_table$Plot = split_code_year[,1]
output_table$Year = split_code_year[,2]

metabolname = paste("metabol","_",metric,sep="")
abundname = paste("abundan","_",metric,sep="")
output_table[,metabolname] = 0
output_table[,abundname] = 0

for(i in 1:length(Samples)){

#Select plot and year
plot = output_table[i,1]
data_sub <- na.omit(data[data[,1]==plot,])

esf1 <- meteESF(spp=data_sub$SpeciesCode,
                abund=data_sub$Count,
                power=data_sub$B^.75,
                minE=min(data_sub$B^.75))

treeIPD <- ipd(esf1)
treeRAD <- sad(esf1)

if("AIC" %in% metric){
output_table[i,metabolname] <- AIC(treeIPD)
output_table[i,abundname] <- AIC(treeRAD)}

if("LL" %in% metric){
invisible(capture.output(output_table[i,metabolname] <- as.numeric(logLik(treeIPD))))
invisible(capture.output(output_table[i,abundname] <- as.numeric(logLik(treeRAD))))}

if("LLZ" %in% metric){invisible(capture.output(output_table[i,
		metabolname] <- as.numeric(logLikZ(treeIPD)$z,nrep=50, return.sim=TRUE)))
invisible(capture.output(output_table[i,abundname] <- as.numeric(logLikZ(treeRAD)$z,nrep=50, return.sim=TRUE)))}

if("KL" %in% metric){
P = rev(sort(treeIPD$r(s <-length(treeIPD$data))))
Q = as.numeric(treeIPD$data)
output_table[i,metabolname] <- kl_divergence(P,Q)

P = rev(sort(treeRAD$r(s <-length(treeRAD$data))))
Q = as.numeric(treeRAD$data)
output_table[i,abundname] <- kl_divergence(P,Q)}

if("MSE" %in% metric){
output_table[i,metabolname] <- mse(treeIPD, type='rank',
									 relative=FALSE)
output_table[i,abundname] <- mse(treeRAD, type='rank',
									 relative=FALSE)}

print((length(Samples)+1)-i)}

output_table$Foresttype = "secondary"
output_table[output_table$Plot == "SV","Foresttype"] = "primary"
output_table[output_table$Plot == "LEPviejo","Foresttype"] = "primary"

cat("Your values have been stored in the object output_table")

return(output_table)

}

#Function to calculate KL divergence needed for the METE_comparison function above
kl_divergence <- function(P, Q, eps = 1e-12) {

# Ensure both vectors have the same length
  if (length(P) != length(Q)) {
    stop("P and Q must have the same length")
  }

# Convert possible NA/NaN to 0 (they shouldn't affect the divergence)
  P[is.na(P)] <- 0
  Q[is.na(Q)] <- 0

# Check that there is at least one positive entry in each vector
  if (all(P == 0) || all(Q == 0)) {
    stop("Both distributions must contain at least one non‑zero element")
  }

# Normalise to proper probability mass functions
  sumP <- sum(P)
  sumQ <- sum(Q)

# Add eps before normalising to keep zeros from staying exactly zero
  p_norm <- (P + eps) / (sumP + eps * length(P))
  q_norm <- (Q + eps) / (sumQ + eps * length(Q))

# Compute the KL divergence
  divergence <- sum(p_norm[p_norm > 0] *
                    log(p_norm[p_norm > 0] / q_norm[p_norm > 0]))

  return(divergence)
}

# -------------------------------------------------------------
# Bootstrap function
# -------------------------------------------------------------

# -------------------------------------------------------------
# Helper 1: format a number of seconds as a readable string
# -------------------------------------------------------------
format_duration <- function(seconds) {
  seconds <- as.numeric(seconds)
  if (is.na(seconds) || !is.finite(seconds)) return("--")
  if (seconds < 60) {
    sprintf("%.0fs", seconds)
  } else if (seconds < 3600) {
    sprintf("%dm %02ds", seconds %/% 60, round(seconds %% 60))
  } else {
    sprintf("%dh %02dm", seconds %/% 3600, (seconds %% 3600) %/% 60)
  }
}

# -------------------------------------------------------------
# Helper 2: Shannon diversity index H from a vector of species
# counts. H = -sum(p_i * ln(p_i)), using only non-zero counts.
# Returns NA if there are no individuals.
# -------------------------------------------------------------
shannon_H <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  total  <- sum(counts)
  if (total <= 0) return(NA_real_)
  p <- counts / total
  -sum(p * log(p))
}

# -------------------------------------------------------------
# Helper 3: aggregate one metric column per PlotCode + Year
# -------------------------------------------------------------
aggregate_robust <- function(df, col) {
     df %>%
      group_by(PlotCode, Year) %>%
      summarise(value = sum(.data[[col]], na.rm = TRUE), .groups = "drop")
  }

# =============================================================
# Main function: bootstrap_treedata()
# -------------------------------------------------------------
# Arguments:
#   data         : data frame with columns PlotCode, Year,
#                  Quadrat, SpeciesCode, Count, B, E, Foresttype,
#                  PlotName
#   n_boot       : number of bootstrap iterations (default 1000)
#   sample_frac  : fraction of quadrats sampled per group
#                  (default 0.80); ignored if n_quadrats is set
#   n_quadrats   : absolute number of quadrats to sample per
#                  group. Overrides sample_frac. NULL = use
#                  sample_frac. Capped at the quadrats available.
#   replace      : sample quadrats with replacement? (default TRUE)
#   normalize    : min-normalize B and E first?      (default TRUE)
#   return_replicates : if TRUE, return a list with both the
#                  summary and the per-iteration replicates;
#                  if FALSE (default), return only the summary
#   output_file  : path to write the summary CSV, or NULL
#   progress_every : print a progress line every k iters (50)
#   verbose      : print progress to the console? (default TRUE)
#
# Returns:
#   return_replicates = FALSE -> a data frame of per-plot means
#                                and bootstrap SEs (the summary).
#   return_replicates = TRUE  -> list(summary = <data frame>,
#                                     replicates = <data frame>).
# =============================================================
bootstrap_treedata <- function(data,
                               n_boot            = 1000,
                               sample_frac       = 0.80,
                               n_quadrats        = NULL,
                               replace           = TRUE,
                               normalize         = TRUE,
                               return_replicates = FALSE,
                               output_file       = NULL,
                               progress_every    = 50,
                               verbose           = TRUE) {

  # verbosity-aware printer: all console output goes through this
  say <- function(...) if (isTRUE(verbose)) cat(...)

  # Input validation
  required_cols <- c("PlotCode", "Year", "Quadrat", "SpeciesCode",
                      "Count", "B", "E", "Foresttype", "PlotName")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Input data is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  if (n_boot < 1) stop("n_boot must be >= 1")
  if (sample_frac <= 0 || sample_frac > 1) stop("sample_frac must be in (0, 1]")
  if (!is.null(n_quadrats) && n_quadrats < 1) stop("n_quadrats must be >= 1")

  mode_str <- if (is.null(n_quadrats)) {
    sprintf("sample_frac: %.2f", sample_frac)
  } else {
    sprintf("n_quadrats: %d", n_quadrats)
  }
  say("=== bootstrap_treedata() ===\n")
  say(sprintf("Iterations: %d | %s | replace: %s ",
              n_boot, mode_str, replace))

  # Step 1: Preparation 
  say("[1/4] Preparing data...\n")

  if (normalize) {
    min_E <- min(data$E, na.rm = TRUE)
    min_B <- min(data$B, na.rm = TRUE)
    if (is.na(min_E) || min_E <= 0) stop("min(E) is NA or <= 0 - cannot normalize E")
    if (is.na(min_B) || min_B <= 0) stop("min(B) is NA or <= 0 - cannot normalize B")
    data$E <- data$E / min_E
    data$B <- data$B / min_B
  }

  data$PlotCode <- as.character(data$PlotCode)
  data$Year     <- as.character(data$Year)
  data$Quadrat  <- as.character(data$Quadrat)

  meta_info <- data %>%
    select(PlotCode, Foresttype, PlotName, Year) %>%
    distinct()

  unique_quads <- data %>%
    select(PlotCode, Year, Quadrat) %>%
    distinct() %>%
    group_by(PlotCode, Year) %>%
    summarise(quadrats = list(unique(Quadrat)), .groups = "drop")

  # Step 2: Bootstrap loop
  say(sprintf("[2/4] Running %d bootstrap iterations...\n", n_boot))
  boot_results_list <- vector("list", n_boot)
  start_time <- Sys.time()

  for (i in 1:n_boot) {

    # Step A: sample quadrats per PlotCode + Year.
    # The sample size is an absolute count (n_quadrats) if given,
    # otherwise a fraction of the available quadrats. quad_instance
    # gives each sampled draw a unique ID so a quadrat drawn more
    # than once (replace = TRUE) counts as multiple observations.
    sampled_quads_df <- unique_quads %>%
      rowwise() %>%
      mutate(
        sampled_quads = list(sample(
          quadrats,
          size = if (is.null(n_quadrats)) {
            ceiling(length(quadrats) * sample_frac)
          } else {
            min(n_quadrats, length(quadrats))
          },
          replace = replace
        ))
      ) %>%
      ungroup() %>%
      select(-quadrats) %>%
      unnest(sampled_quads) %>%
      rename(Quadrat = sampled_quads) %>%
      mutate(quad_instance = row_number())

    # Step B: keep only trees from sampled quadrats from the original data
    sampled_data <- data %>%
      inner_join(sampled_quads_df,
                 by = c("PlotCode", "Year", "Quadrat"),
                 relationship = "many-to-many")

    # Step C: richness (S), abundance (N), Shannon diversity (H)
    # and Hill evenness (exp(H)/S) per PlotCode + Year.
    SN_boot <- sampled_data %>%
      group_by(PlotCode, Year, SpeciesCode) %>%
      summarise(sp_count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
      group_by(PlotCode, Year) %>%
      summarise(
        S = sum(sp_count > 0),
        N = sum(sp_count),
        H = shannon_H(sp_count),
        .groups = "drop"
      ) %>%
      mutate(
        Evenness = ifelse(S > 0, exp(H) / S, NA_real_)
      )

    # Step D: robust aggregation of B and E per PlotCode + Year
    B_boot <- aggregate_robust(sampled_data, "B") %>%
      rename(B_sum = value)
    E_boot <- aggregate_robust(sampled_data, "E") %>%
      rename(E_sum = value)

    # Step E: merge metrics with metadata
    boot_results_list[[i]] <- meta_info %>%
      left_join(SN_boot, by = c("PlotCode", "Year")) %>%
      left_join(B_boot,  by = c("PlotCode", "Year")) %>%
      left_join(E_boot,  by = c("PlotCode", "Year")) %>%
      rename(B = B_sum, E = E_sum)

    # Progress + ETA
    if (verbose && (i %% progress_every == 0 || i == n_boot)) {
      elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      per_iter  <- elapsed / i
      remaining <- per_iter * (n_boot - i)
      say(sprintf("  iter %d/%d (%.0f%%) | elapsed %s | ETA %s\n",
                  i, n_boot, 100 * i / n_boot,
                  format_duration(elapsed), format_duration(remaining)))
    }
  }

  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  say(sprintf("  Loop finished in %s.\n", format_duration(total_time)))

  # Step 3: Aggregate
  say("[3/4] Aggregating means and bootstrap standard errors...\n")
  all_boot_data <- bind_rows(boot_results_list, .id = "iteration")

  final_aggdata <- all_boot_data %>%
    group_by(PlotCode, Foresttype, PlotName, Year) %>%
    summarise(
      S_mean = mean(S, na.rm = TRUE),
      S_se   = sd(S,   na.rm = TRUE),
      N_mean = mean(N, na.rm = TRUE),
      N_se   = sd(N,   na.rm = TRUE),
      H_mean = mean(H, na.rm = TRUE),
      H_se   = sd(H,   na.rm = TRUE),
      Evenness_mean = mean(Evenness, na.rm = TRUE),
      Evenness_se   = sd(Evenness,   na.rm = TRUE),
      B_mean = mean(B, na.rm = TRUE),
      B_se   = sd(B,   na.rm = TRUE),
      E_mean = mean(E, na.rm = TRUE),
      E_se   = sd(E,   na.rm = TRUE),
      .groups = "drop"
    )

  # Step 4: Save
  if (!is.null(output_file)) {
    write.csv(final_aggdata, file = output_file, row.names = FALSE)
    say(sprintf("[4/4] Results written to: %s\n", output_file))
  } else {
    say("[4/4] No output_file given - results returned only.\n")
  }

  say("=== Done. ===\n")

  if (isTRUE(return_replicates)) {
    return(list(summary = final_aggdata, replicates = all_boot_data))
  }
  return(final_aggdata)
}


# =============================================================
# Wrapper: accumulation_curve_B()
# -------------------------------------------------------------
# Runs bootstrap_treedata() across a gradient of sampled-quadrat
# counts for each Plot x Year combination to show how the 
# plot-level basal area estimate and its
# uncertainty behave as sampling intensity increases.
#
# For each quadrat count k the sampled B (a sum over k quadrats)
# is rescaled to a FULL-PLOT estimate by multiplying by
# (full quadrat count / k). The rescaled estimate is unbiased at
# every k - so its bootstrap mean stays roughly flat - while the
# 95% CI narrows as k grows. The narrowing CI is the
# representativeness signal: it shows how many quadrats are
# needed for a precise plot-level B.
#
# Subsampling is done WITHOUT replacement (in contrast to the
# bootstrapping function above so the x-axis is
# literally "number of distinct quadrats sampled".
#
# Arguments:
#   data      : tree data (see bootstrap_treedata)
#   quad_grid : vector of quadrat counts to evaluate
#   n_boot    : bootstrap iterations per quadrat level (default 500)
#   robust    : aggregation mode passed through (default "none")
#   normalize : min-normalize B / E? (default TRUE)
#
# Returns: a long data frame with one row per PlotCode x Year x
# quadrat level, including B_estimate (full-plot scale), B_est_se,
# and the 95% CI (B_lower, B_upper).
# =============================================================
accumulation_curve_B <- function(data,
                                 quad_grid = c(5, 10, 20, 30, 40, 50,
                                               60, 70, 80, 90, 100),
                                 n_boot    = 500,
                                 normalize = TRUE,
                                 ...) {

  # full quadrat count per PlotCode + Year (for the rescaling)
  plot_n <- data %>%
    mutate(PlotCode = as.character(PlotCode),
           Year     = as.character(Year),
           Quadrat  = as.character(Quadrat)) %>%
    distinct(PlotCode, Year, Quadrat) %>%
    count(PlotCode, Year, name = "full_n")

  cat(sprintf("=== accumulation_curve_B(): %d quadrat levels x %d iterations ===\n",
              length(quad_grid), n_boot))
  start_time <- Sys.time()

  out <- vector("list", length(quad_grid))
  for (j in seq_along(quad_grid)) {
    k <- quad_grid[j]
    cat(sprintf("[%d/%d] Bootstrapping with %d quadrats...\n",
                j, length(quad_grid), k))
    res <- bootstrap_treedata(
      data       = data,
      n_boot     = n_boot,
      n_quadrats = k,
      replace    = FALSE,        # subsample distinct quadrats
      normalize  = normalize,
      verbose    = FALSE,
      ...
    )
    res$n_quadrats <- k
    out[[j]] <- res
  }
  cat(sprintf("Accumulation analysis finished in %s.\n",
              format_duration(as.numeric(difftime(Sys.time(), start_time,
                                                   units = "secs")))))

  # combine; rescale summed B to a full-plot estimate
  bind_rows(out) %>%
    left_join(plot_n, by = c("PlotCode", "Year")) %>%
    mutate(
      k_eff      = pmin(n_quadrats, full_n),
      scale      = full_n / k_eff,
      B_estimate = B_mean * scale,
      B_est_se   = B_se   * scale,
      B_lower    = B_estimate - 1.96 * B_est_se,
      B_upper    = B_estimate + 1.96 * B_est_se
    )
}
