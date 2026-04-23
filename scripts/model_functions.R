# model_functions.R — config loading + model function definitions
# Source this file before any execution. Sets globals: paramsData, n_t, n_c,
# discount, config_dist, config_dist_p1, config_dist_p2, strategies (via create_*).
# Requires working directory = project root (thesis/).

# -----------------------------------------------------------------------------
# Load all model parameters from config.csv
# Each row in config.csv specifies a parameter value, its probability
# distribution (for PSA), and its source reference.
# -----------------------------------------------------------------------------
config <- read.csv("input/config.csv", stringsAsFactors = FALSE)
config_list    <- setNames(config$value,                    config$parameter)
config_dist    <- setNames(config$distribution,             config$parameter)
config_dist_p1 <- setNames(as.numeric(config$dist_param1), config$parameter)
config_dist_p2 <- setNames(as.numeric(config$dist_param2), config$parameter)

# -----------------------------------------------------------------------------
# Model structural parameters
# Time horizon: 55 years (660 monthly cycles), entry age ~25 to age ~80 (lifetime; NICE reference case)
# Cohort size: 100,000 hypothetical migrants
# Discount rate: 3.5% per annum, applied to both costs and effects (NICE TA)
# -----------------------------------------------------------------------------
n_t      <- as.integer(config_list["n_t"])
n_c      <- as.integer(config_list["n_c"])
discount <- as.numeric(config_list["discount"])
discount_m   <- (1 + discount)^(1/12) - 1          # monthly compound discount rate
disc_factors <- 1 / (1 + discount_m)^seq_len(n_t)  # pre-computed 660-element discount vector

# ONS National Life Tables 2020-2022 (England & Wales, persons, both sexes).
# All-cause mortality rates per 100,000/year, by 5-year age band.
# Cohort entry age ~25; band 1 = 25-29, ..., band 11 = 75-79.
# 11 bands x 60 cycles (5yr) = 660 cycles = 55yr lifetime horizon.
# Converted to monthly probability: rate / 100,000 / 12.
# Source: ONS National Life Tables 2020-2022, Table 3 (England & Wales, persons).
# Values in config.csv (parameters ons_mort_25_29 through ons_mort_75_79).
ons_bg_mortality_monthly <- as.numeric(config_list[
  c("ons_mort_25_29","ons_mort_30_34","ons_mort_35_39","ons_mort_40_44",
    "ons_mort_45_49","ons_mort_50_54","ons_mort_55_59","ons_mort_60_64",
    "ons_mort_65_69","ons_mort_70_74","ons_mort_75_79")
]) / 1e5 / 12

# 16 health states representing the TB disease and care cascade.
# States are grouped into: uninfected, latent TB (LTBI), active TB, and dead.
# Within LTBI and active TB, states reflect diagnostic and treatment status.
v_state_names <- c(
  "Uninfected",
  "LatentUndiagnosed", "LatentDiagnosed",
  "ActiveUndiagnosed", "ActiveDiagnosed",
  "LatentTreated", "LatentNotreated",
  "ActiveTreated", "ActiveNotreated",
  "LatentCompleted", "LatentDiscontinued", "LatentLtfu",
  "ActiveCompleted", "ActiveDiscontinued", "ActiveLtfu",
  "Dead"
)

# Extract the subset of parameters used directly inside the model() function
# (state costs, QALY weights, transition probabilities). Other parameters
# such as initial conditions and effective detection rates are read separately.
param_names  <- names(config_list)
model_params <- param_names[grepl("^(Cost_|qaly_|p_)", param_names)]
paramsData   <- as.list(as.numeric(config_list[model_params]))
names(paramsData) <- model_params

# -----------------------------------------------------------------------------
# Core Markov model function
#
# Implements a monthly-cycle Markov cohort model with 16 health states.
# The model accepts an optional initial state distribution (init_dist) so
# that strategy-specific starting conditions (from the diagnostic decision
# tree) can be passed in. If no init_dist is provided, the model uses the
# background prevalence values from config.csv (no-screening baseline).
#
# Returns: state membership trace, discounted costs, discounted QALYs,
#          transition matrix, and per-state payoff vectors.
# -----------------------------------------------------------------------------
model <- function(.params, init_dist = NULL, n_t_override = NULL) {
  # Allow extended-horizon scenario analyses to shadow the global n_t without
  # modifying global state. If n_t_override is supplied, it is used for all
  # cycle counts in this call; otherwise the global n_t is used unchanged.
  if (!is.null(n_t_override)) n_t <- n_t_override
  with(.params, {
    n_s <- length(v_state_names)

    # Construct the 16x16 monthly transition probability matrix.
    # Rows = origin state, columns = destination state.
    # Off-diagonal entries are filled from config.csv parameters;
    # diagonal entries are set so each row sums to exactly 1.
    m_p <- matrix(0, nrow = n_s, ncol = n_s, dimnames = list(from = v_state_names, to = v_state_names))

    # -------------------------------------------------------------------------
    # Complete list of allowed state transitions (42 connections total).
    # All other cell entries remain zero (transitions not listed are forbidden).
    #
    # ENTRY (initial allocation at t = 0, handled outside model() via init_dist):
    #   U  -> U          uninfected remain uninfected
    #   L  -> Lu         latent TB missed by screening (proportion = 1 - eff_ltbi)
    #   L  -> Ld         latent TB detected by screening (proportion = eff_ltbi)
    #   I  -> Iu         active TB missed by screening (proportion = 1 - eff_active)
    #   I  -> Id         active TB detected by screening (proportion = eff_active)
    #
    # UNINFECTED (U):
    #   U  -> Lu         new LTBI acquisition (background transmission rate)
    #   U  -> Dead       background mortality (ONS)
    #
    # LATENT UNDIAGNOSED (Lu):
    #   Lu -> Ld         opportunistic/background screening detection
    #   Lu -> Iu         progression from latent to active TB (undiagnosed; observed incidence)
    #   Lu -> Dead       background mortality
    #
    # LATENT DIAGNOSED (Ld):
    #   Ld -> LTt        starts preventive treatment
    #   Ld -> LTn        declines / not offered treatment
    #   Ld -> Iu         reactivates despite diagnosis (before treatment starts)
    #   Ld -> Dead       background mortality
    #
    # LATENT ON TREATMENT (LTt):
    #   LTt -> LTOc      completes treatment course
    #   LTt -> LTOd      discontinues (adverse effects, ~9% UK)
    #   LTt -> LTOltfu   lost to follow-up during treatment
    #   LTt -> Dead      background + hepatotoxicity mortality
    #
    # LATENT NOT TREATED (LTn):
    #   LTn -> Iu        progression from latent to active TB (elevated risk, ~1%/year)
    #   LTn -> Dead      background mortality
    #
    # LATENT COMPLETED (LTOc):
    #   LTOc -> U        protected / effectively cured (treatment 60-90% effective)
    #   LTOc -> Dead     background mortality
    #
    # LATENT DISCONTINUED (LTOd):
    #   LTOd -> Iu       progression from latent to active TB (partial protection lost)
    #   LTOd -> Ld       re-engages with care
    #   LTOd -> Dead     background mortality
    #
    # LATENT LTFU (LTOltfu):
    #   LTOltfu -> Iu    progression from latent to active TB (no treatment protection)
    #   LTOltfu -> Ld    re-engages with care
    #   LTOltfu -> Dead  background mortality
    #
    # ACTIVE UNDIAGNOSED (Iu):
    #   Iu -> Id         diagnosed (spontaneous presentation / contact tracing)
    #   Iu -> Dead       elevated TB mortality (untreated ~50% 5-year mortality)
    #
    # ACTIVE DIAGNOSED (Id):
    #   Id -> ITt        starts treatment (2HRZE/4HR)
    #   Id -> ITn        declines / delays treatment
    #   Id -> Dead       pre-treatment mortality
    #
    # ACTIVE ON TREATMENT (ITt):
    #   ITt -> ITOc      completes treatment (84.4% UK 2023)
    #   ITt -> ITOd      discontinues treatment (~2% UK)
    #   ITt -> ITOltfu   lost to follow-up during treatment
    #   ITt -> Dead      treatment-phase mortality
    #
    # ACTIVE NOT TREATED (ITn):
    #   ITn -> Id        returns to care
    #   ITn -> Dead      high mortality without treatment
    #
    # ACTIVE COMPLETED (ITOc):
    #   ITOc -> U        cured (88% treatment success UK; post-TB sequelae captured via utility)
    #   ITOc -> Dead     post-treatment mortality (near background)
    #
    # ACTIVE DISCONTINUED (ITOd):
    #   ITOd -> Iu       returns to infectious / undiagnosed pool
    #   ITOd -> Dead     elevated mortality
    #
    # ACTIVE LTFU (ITOltfu):
    #   ITOltfu -> Iu    returns to infectious pool
    #   ITOltfu -> Dead  elevated mortality
    #
    # DEAD:
    #   Dead -> Dead     absorbing state
    # -------------------------------------------------------------------------

    # --- Disease progression and care cascade transitions ---
    m_p["Uninfected", "LatentUndiagnosed"] <- p_Uninfected_LatentUndiagnosed
    m_p["LatentUndiagnosed", "LatentDiagnosed"] <- p_LatentUndiagnosed_LatentDiagnosed
    m_p["LatentUndiagnosed", "ActiveUndiagnosed"] <- p_LatentUndiagnosed_ActiveUndiagnosed
    m_p["LatentDiagnosed", "ActiveUndiagnosed"] <- p_LatentDiagnosed_ActiveUndiagnosed
    m_p["ActiveUndiagnosed", "ActiveDiagnosed"] <- p_ActiveUndiagnosed_ActiveDiagnosed
    m_p["LatentDiagnosed", "LatentTreated"] <- p_LatentDiagnosed_LatentTreated
    m_p["LatentDiagnosed", "LatentNotreated"] <- p_LatentDiagnosed_LatentNotreated
    m_p["ActiveDiagnosed", "ActiveTreated"] <- p_ActiveDiagnosed_ActiveTreated
    m_p["ActiveDiagnosed", "ActiveNotreated"] <- p_ActiveDiagnosed_ActiveNotreated
    m_p["LatentTreated", "LatentCompleted"] <- p_LatentTreated_LatentCompleted
    m_p["LatentTreated", "LatentDiscontinued"] <- p_LatentTreated_LatentDiscontinued
    m_p["LatentTreated", "LatentLtfu"] <- p_LatentTreated_LatentLtfu
    m_p["ActiveTreated", "ActiveCompleted"] <- p_ActiveTreated_ActiveCompleted
    m_p["ActiveTreated", "ActiveDiscontinued"] <- p_ActiveTreated_ActiveDiscontinued
    m_p["ActiveTreated", "ActiveLtfu"] <- p_ActiveTreated_ActiveLtfu
    m_p["LatentNotreated", "ActiveUndiagnosed"] <- p_LatentNotreated_ActiveUndiagnosed
    m_p["ActiveDiscontinued", "ActiveUndiagnosed"] <- p_ActiveDiscontinued_ActiveUndiagnosed
    m_p["ActiveLtfu", "ActiveUndiagnosed"] <- p_ActiveLtfu_ActiveUndiagnosed
    m_p["LatentCompleted", "Uninfected"]        <- p_LatentCompleted_Uninfected
    # Residual reactivation after LTBI treatment: 14% of untreated rate (Berrocal-Almanza 2022: HR 0.14)
    m_p["LatentCompleted", "ActiveUndiagnosed"] <- p_LatentCompleted_ActiveUndiagnosed
    m_p["ActiveCompleted", "Uninfected"] <- p_ActiveCompleted_Uninfected

    # --- Progression from latent to active TB (observed incidence) and return-to-care transitions ---
    m_p["LatentDiscontinued", "ActiveUndiagnosed"] <- p_LatentDiscontinued_ActiveUndiagnosed
    m_p["LatentLtfu", "ActiveUndiagnosed"]         <- p_LatentLtfu_ActiveUndiagnosed
    m_p["ActiveNotreated", "ActiveDiagnosed"]      <- p_ActiveNotreated_ActiveDiagnosed

    # Re-engagement with care after LTBI treatment discontinuation or LTFU
    if ("p_LatentDiscontinued_LatentDiagnosed" %in% names(.params))
      m_p["LatentDiscontinued", "LatentDiagnosed"] <- .params[["p_LatentDiscontinued_LatentDiagnosed"]]
    if ("p_LatentLtfu_LatentDiagnosed" %in% names(.params))
      m_p["LatentLtfu", "LatentDiagnosed"] <- .params[["p_LatentLtfu_LatentDiagnosed"]]

    # --- Mortality transitions ---
    # State-specific monthly probabilities of death are read from config.csv.
    # Active TB states have elevated mortality; LTBI states use background rate.
    for(st in v_state_names[v_state_names != "Dead"]) {
      dead_par_name <- paste0("p_", st, "_Dead")
      if (dead_par_name %in% names(.params)) {
        m_p[st, "Dead"] <- .params[[dead_par_name]]
      }
    }
    # Dead is an absorbing state — once entered, no further transitions occur.
    m_p["Dead", "Dead"] <- 1

    # Diagonal = 1 − row sum (each row must sum to 1).
    off <- rowSums(m_p)
    for(i in seq_len(n_s)) {
      st <- v_state_names[i]
      if(st != "Dead") diag(m_p)[i] <- 1 - off[i]
    }

    # Numerical safeguard: during PSA, sampled parameter combinations can
    # occasionally cause row sums to exceed 1 (diagonal goes negative).
    # In this case, outgoing probabilities are rescaled proportionally so
    # the row still sums to 1. This is standard practice in PSA for Markov models.
    for(i in seq_len(n_s)) {
      st <- v_state_names[i]
      if(st == "Dead") next
      if(diag(m_p)[i] < 0) {
        off_diag_sum <- off[i]
        for(j in seq_len(n_s)) {
          if(i != j) m_p[i, j] <- m_p[i, j] / off_diag_sum
        }
        diag(m_p)[i] <- 0
      }
    }

    # --- Run the Markov cohort through n_t monthly cycles ---
    state_membership <- matrix(0, nrow = n_t, ncol = n_s, dimnames = list(1:n_t, v_state_names))

    # Populate cycle 1 from the strategy-specific initial distribution (if
    # provided) or from background prevalence values in config.csv (used for
    # the no-screening base case analysis).
    if (!is.null(init_dist)) {
      state_membership[1, ] <- init_dist
    } else {
      init_Uninfected        <- as.numeric(config_list["init_Uninfected"])
      init_LatentUndiagnosed <- as.numeric(config_list["init_LatentUndiagnosed"])
      init_ActiveUndiagnosed <- as.numeric(config_list["init_ActiveUndiagnosed"])
      state_membership[1, ]                        <- rep(0, n_s)
      state_membership[1, "Uninfected"]            <- n_c * init_Uninfected
      state_membership[1, "LatentUndiagnosed"]     <- n_c * init_LatentUndiagnosed
      state_membership[1, "ActiveUndiagnosed"]     <- n_c * init_ActiveUndiagnosed
    }

    # --- Age-varying background mortality (ONS 5-year age bands) ---
    # The cohort enters at age ~25. Background mortality increases with age;
    # the transition probability for background-mortality states is updated
    # every 60 cycles (5 years) to the relevant ONS age-band rate.
    # Active TB states keep their fixed TB-specific mortality (dominates background).
    # LatentTreated uses background + hepatotoxicity excess from NICE NG33.
    bg_mort_base     <- .params[["p_Uninfected_Dead"]]
    hepatotox_excess <- max(0, .params[["p_LatentTreated_Dead"]] - bg_mort_base)
    # ActiveCompleted is included so post-treatment mortality rises with age (ONS bands).
    # Omitting it caused a fixed 0.0002/mo rate that falls BELOW background at ages >65.
    bg_states <- c("Uninfected", "LatentUndiagnosed", "LatentDiagnosed",
                   "LatentNotreated", "LatentCompleted", "LatentDiscontinued", "LatentLtfu",
                   "ActiveCompleted")

    # Initialise matrix to age-band 1 (ages 25-29) before entering the loop.
    current_band <- 1L
    new_bg <- ons_bg_mortality_monthly[current_band]
    for (st in bg_states) {
      delta <- m_p[st, "Dead"] - new_bg
      m_p[st, "Dead"] <- new_bg
      m_p[st, st]     <- m_p[st, st] + delta
    }
    lt_delta <- m_p["LatentTreated", "Dead"] - (new_bg + hepatotox_excess)
    m_p["LatentTreated", "Dead"]          <- new_bg + hepatotox_excess
    m_p["LatentTreated", "LatentTreated"] <- m_p["LatentTreated", "LatentTreated"] + lt_delta

    # Two-phase reactivation: all latent→active transitions scale by
    # (p_react_phase2 / phase1_rate) at cycle 61 (year 5).
    # To disable: set p_react_phase2 = p_LatentUndiagnosed_ActiveUndiagnosed in config.csv.
    react_phase1         <- .params[["p_LatentUndiagnosed_ActiveUndiagnosed"]]
    p_react_phase2       <- .params[["p_react_phase2"]]
    react_phase_switched <- FALSE
    react_latent_states  <- c("LatentUndiagnosed", "LatentDiagnosed", "LatentNotreated",
                               "LatentCompleted",   "LatentDiscontinued", "LatentLtfu")

    # Matrix multiplication advances the cohort one cycle at a time.
    # At each 5-year band boundary (cycles 61, 121, ..., 601) the background
    # mortality rows are updated to the next ONS age-band rate (11 bands total,
    # covering ages 25-79; band 11 applies for ages 75-80 at end of horizon).
    for(t in 2:n_t) {
      new_band <- min(11L, as.integer(floor((t - 1) / 60)) + 1L)
      if (new_band != current_band) {
        current_band <- new_band
        new_bg <- ons_bg_mortality_monthly[current_band]
        for (st in bg_states) {
          delta <- m_p[st, "Dead"] - new_bg
          m_p[st, "Dead"] <- new_bg
          m_p[st, st]     <- m_p[st, st] + delta
        }
        lt_delta <- m_p["LatentTreated", "Dead"] - (new_bg + hepatotox_excess)
        m_p["LatentTreated", "Dead"]          <- new_bg + hepatotox_excess
        m_p["LatentTreated", "LatentTreated"] <- m_p["LatentTreated", "LatentTreated"] + lt_delta
      }
      # Switch reactivation to phase 2 rate at year 5 (cycle 61)
      if (!react_phase_switched && t == 61L) {
        react_phase_switched <- TRUE
        react_scale <- p_react_phase2 / react_phase1
        for (st in react_latent_states) {
          old_r            <- m_p[st, "ActiveUndiagnosed"]
          new_r            <- old_r * react_scale
          m_p[st, "ActiveUndiagnosed"] <- new_r
          m_p[st, st]      <- m_p[st, st] + (old_r - new_r)
        }
      }
      state_membership[t, ] <- as.numeric(state_membership[t-1, ] %*% m_p)
      if (abs(sum(state_membership[t, ]) - n_c) > 1.0)
        warning(sprintf("Cycle %d: cohort sum = %.1f (expected %.0f; check transition matrix row sums)", t, sum(state_membership[t, ]), n_c))
    }

    # --- Per-state monthly cost and QALY payoff vectors ---
    # Costs are in 2024 GBP. QALYs are divided by 12 to
    # convert annual utility weights to monthly values.
    m_payoffsCost <- matrix(0, nrow = n_s, ncol = 1, dimnames = list(v_state_names, "Costs"))
    m_payoffsQALY <- matrix(0, nrow = n_s, ncol = 1, dimnames = list(v_state_names, "QALYs"))
    m_payoffsCost[,"Costs"] <- c(
      Cost_Uninfected,
      Cost_LatentUndiagnosed,Cost_LatentDiagnosed,
      Cost_ActiveUndiagnosed,Cost_ActiveDiagnosed,
      Cost_LatentTreated,Cost_LatentNotreated,
      Cost_ActiveTreated,Cost_ActiveNotreated,
      Cost_LatentCompleted,Cost_LatentDiscontinued,Cost_LatentLtfu,
      Cost_ActiveCompleted,Cost_ActiveDiscontinued,Cost_ActiveLtfu,
      Cost_Dead
    )
    m_payoffsQALY[,"QALYs"] <- c(
      qaly_Uninfected,
      qaly_LatentUndiagnosed,qaly_LatentDiagnosed,
      qaly_ActiveUndiagnosed,qaly_ActiveDiagnosed,
      qaly_LatentTreated,qaly_LatentNotreated,
      qaly_ActiveTreated,qaly_ActiveNotreated,
      qaly_LatentCompleted,qaly_LatentDiscontinued,qaly_LatentLtfu,
      qaly_ActiveCompleted,qaly_ActiveDiscontinued,qaly_ActiveLtfu,
      qaly_Dead
    ) / 12

    # Apply discounting at 3.5% per annum (converted to monthly compound rate)
    # per NICE reference case. Future costs and QALYs are worth less than
    # immediate ones; dividing by (1 + r)^t adjusts for this.
    payoff_traceCost <- state_membership %*% m_payoffsCost
    payoff_traceQaly <- state_membership %*% m_payoffsQALY

    # Vectorised discounting: multiply by pre-computed discount vector.
    # If n_t was overridden (scenario analysis), compute locally; otherwise use
    # the global disc_factors (660-element vector, computed once at startup).
    dv <- if (n_t == length(disc_factors)) disc_factors else 1/(1+discount_m)^seq_len(n_t)
    payoff_traceCost_d <- payoff_traceCost * dv
    payoff_traceQaly_d <- payoff_traceQaly * dv

    # Half-cycle correction: assumes transitions occur on average mid-cycle
    # rather than at the start. Subtracts half the first and last cycle values
    # to correct for the overcount introduced by the standard Markov assumption.
    total_cost_hcc <- sum(payoff_traceCost_d) -
      0.5 * payoff_traceCost_d[1,1] - 0.5 * payoff_traceCost_d[n_t,1]
    total_qaly_hcc <- sum(payoff_traceQaly_d) -
      0.5 * payoff_traceQaly_d[1,1] - 0.5 * payoff_traceQaly_d[n_t,1]

    list(
      state_membership = state_membership,
      payoff_trace_perCycle_cost = payoff_traceCost_d,
      payoff_trace_perCycle_qaly = payoff_traceQaly_d,
      total_cost = total_cost_hcc,
      total_qaly = total_qaly_hcc,
      m_p = m_p,
      m_payoffsCost = m_payoffsCost,
      m_payoffsQALY = m_payoffsQALY
    )
  })
}

# -----------------------------------------------------------------------------
# PSA parameter sampling
#
# Draws one random sample of all model parameters from the probability
# distributions specified in config.csv. Beta distributions are used for
# probabilities and utilities (bounded 0-1); gamma distributions are used
# for costs (bounded >0, right-skewed). Fixed parameters are not varied.
# This function is called once per PSA simulation inside run_psa().
# -----------------------------------------------------------------------------
sample_params <- function() {
  sampled <- paramsData  # start from base case values

  # Vectorised sampling: draw all beta and gamma parameters in two batch calls.
  # R's rbeta/rgamma recycle shape vectors, sampling one value per parameter set.
  beta_names  <- names(sampled)[config_dist[names(sampled)] == "beta"  & !is.na(config_dist[names(sampled)])]
  gamma_names <- names(sampled)[config_dist[names(sampled)] == "gamma" & !is.na(config_dist[names(sampled)])]

  if (length(beta_names) > 0) {
    sampled[beta_names] <- as.list(rbeta(length(beta_names),
                                         shape1 = config_dist_p1[beta_names],
                                         shape2 = config_dist_p2[beta_names]))
  }
  if (length(gamma_names) > 0) {
    sampled[gamma_names] <- as.list(rgamma(length(gamma_names),
                                            shape = config_dist_p1[gamma_names],
                                            scale = config_dist_p2[gamma_names]))
  }

  # Correlated asymptomatic-state QALY weights.
  # All latent TB states and Uninfected share the same underlying utility (0.92);
  # independent draws create artificial within-sim variation. One shared draw.
  asym_names <- c("qaly_Uninfected", "qaly_LatentUndiagnosed", "qaly_LatentDiagnosed",
                  "qaly_LatentTreated", "qaly_LatentNotreated", "qaly_LatentCompleted",
                  "qaly_LatentDiscontinued", "qaly_LatentLtfu")
  sampled[asym_names] <- rbeta(1,
                                shape1 = config_dist_p1["qaly_Uninfected"],
                                shape2 = config_dist_p2["qaly_Uninfected"])

  # Correlated active-TB state QALY weights (Group A).
  # Five active-TB states share the same source utility (0.68, Jit 2011; Beta(14.12,6.64));
  # independent draws create spurious variation in incremental active-TB detection benefit.
  # qaly_ActiveTreated (0.81) and qaly_ActiveCompleted (0.92) have distinct distributions
  # and are intentionally excluded.
  active_tb_names <- c("qaly_ActiveUndiagnosed", "qaly_ActiveDiagnosed",
                       "qaly_ActiveNotreated", "qaly_ActiveDiscontinued",
                       "qaly_ActiveLtfu")
  active_tb_names <- active_tb_names[active_tb_names %in% beta_names]
  if (length(active_tb_names) > 0) {
    sampled[active_tb_names] <- rbeta(1,
                                      shape1 = config_dist_p1["qaly_ActiveUndiagnosed"],
                                      shape2 = config_dist_p2["qaly_ActiveUndiagnosed"])
  }

  # Correlated untreated-LTBI reactivation rates — early phase (Group B).
  # Four states share the same Berrocal-Almanza 2022 early-phase rate (0.000912/mo;
  # Beta(9.23,10114.77)); independent draws create spurious variation in the QALY
  # gain from preventing reactivation — the central mechanism of LTBI treatment benefit.
  # p_LatentDiscontinued_ActiveUndiagnosed has a distinct distribution and is excluded.
  react_early_names <- c("p_LatentUndiagnosed_ActiveUndiagnosed",
                          "p_LatentDiagnosed_ActiveUndiagnosed",
                          "p_LatentNotreated_ActiveUndiagnosed",
                          "p_LatentLtfu_ActiveUndiagnosed")
  react_early_names <- react_early_names[react_early_names %in% beta_names]
  if (length(react_early_names) > 0) {
    sampled[react_early_names] <- rbeta(1,
                                        shape1 = config_dist_p1["p_LatentUndiagnosed_ActiveUndiagnosed"],
                                        shape2 = config_dist_p2["p_LatentUndiagnosed_ActiveUndiagnosed"])
  }

  # Correlated background mortality rates — non-TB (Group C).
  # Seven states share the same ONS age-standardised population mortality
  # (Beta(17.39,416982)); independent draws add unnecessary noise across identical-
  # source states. p_LatentTreated_Dead has a distinct distribution and is excluded.
  bg_mort_names <- c("p_Uninfected_Dead", "p_LatentUndiagnosed_Dead",
                     "p_LatentDiagnosed_Dead", "p_LatentNotreated_Dead",
                     "p_LatentCompleted_Dead", "p_LatentDiscontinued_Dead",
                     "p_LatentLtfu_Dead")
  bg_mort_names <- bg_mort_names[bg_mort_names %in% beta_names]
  if (length(bg_mort_names) > 0) {
    sampled[bg_mort_names] <- rbeta(1,
                                    shape1 = config_dist_p1["p_Uninfected_Dead"],
                                    shape2 = config_dist_p2["p_Uninfected_Dead"])
  }

  return(sampled)
}

#------------------------------------------------------------------------------#
# Diagnostic strategy initial conditions
#
# Active TB allocation: TP and FN taken directly from the Zenner et al. 2025 ERJ
#   decision tree (Initial_conditions.xlsx, 100k cohort). These account for the
#   full pathway (symptom screen → test → culture confirmation).
#   TP → ActiveDiagnosed | FN → ActiveUndiagnosed
#   FP + TN → Uninfected (UK assumption: culture confirmation before treatment)
#
# LTBI allocation: taken directly from Zenner et al. 2025 ERJ decision tree (Excel col 5).
#   LatentUndiagnosed = 17,800 (Berrocal-Almanza 2022: 17.8% per 100k) − Excel detected.
#   Excel LTBI column used directly; consistent with decision tree design.
#   The decision tree embeds its own lower LTBI prevalence assumption (~0.1–0.15%),
#   which differs from the Berrocal-Almanza 2022 17.8% background prevalence used elsewhere.
#
# Upfront diagnostic costs: Tot_costs from Excel, per-person = Tot_costs / 100,000.
#------------------------------------------------------------------------------#

# Background prevalence at entry — read from config.csv (UKHSA 2024)
prev_uninfected <- as.numeric(config_list["init_Uninfected"])
prev_ltbi       <- as.numeric(config_list["init_LatentUndiagnosed"])
prev_active     <- as.numeric(config_list["init_ActiveUndiagnosed"])

# IGRA LTBI sensitivities — read from config.csv (Zenner et al. 2025 Eur Respir J Table 1)
igra_sens_qft  <- as.numeric(config_list["igra_ltbi_sens_qft"])   # QFT-TB Gold Plus: 0.83
igra_sens_tspt <- as.numeric(config_list["igra_ltbi_sens_tspt"])  # T-SPOT.TB: 0.88

# IGRA LTBI specificities — read from config.csv (Pai et al. 2008, BCG-vaccinated populations)
# QFT: 0.96 (95% CI 0.94-0.98); T-SPOT: 0.93 (95% CI 0.86-1.00)
igra_spec_qft  <- as.numeric(config_list["igra_spec_qft"])        # QFT specificity:   0.96
igra_spec_tspt <- as.numeric(config_list["igra_spec_tspt"])       # T-SPOT specificity: 0.93

# IGRA uptake — base case proportion of LTBI-positive migrants accepting IGRA test.
# Real-world: 60-80% (UKHSA 2024). Base case: 75% (midpoint). SA: 50-100% in run_dsa.R.
igra_uptake_base <- as.numeric(config_list["igra_uptake_base"])   # 0.75

# eff_active parameter names (for PSA sampling — all have beta distributions in config.csv)
eff_active_param_names <- c("eff_active_CXR_only", "eff_active_QFT",
                              "eff_active_QFT_allsx_ultra", "eff_active_QFT_allsx_xpert",
                              "eff_active_QFT_cough_ultra", "eff_active_QFT_cough_xpert",
                              "eff_active_TSPOT",   "eff_active_TST",
                              "eff_active_CXR_QFT", "eff_active_TST_QFT",
                              "eff_active_CXR_Xpert")
eff_active_base <- setNames(
  as.numeric(config_list[eff_active_param_names]),
  eff_active_param_names
)
eff_active_dist_p1 <- setNames(
  as.numeric(config_dist_p1[eff_active_param_names]),
  eff_active_param_names
)
eff_active_dist_p2 <- setNames(
  as.numeric(config_dist_p2[eff_active_param_names]),
  eff_active_param_names
)

# Read decision tree outputs from Excel (active TB TP/FN and diagnostic costs)
# Columns (after skipping 2 header rows): 1=strategy, 4=TP, 5=LTBI, 6=Tot_costs, 9=FN
excel_ic <- read_excel("input/Screening_Decision_Tree.xlsx", col_names = FALSE, skip = 2)
excel_denom <- 100000  # Excel cohort denominator (TP+FN+TN+FP = 100,000)

# Column index validation — guards against silent errors if Excel columns shift
# Expected layout (after skip=2): col1=strategy, col4=TP, col5=LTBI, col6=Tot_costs, col9=FN
# Col 4 (TP): active TB true positives per 100k — expect 100–1,000
# Col 6 (Tot_costs): programme costs per 100k — expect >100,000 (typically 400k–19M)
# Col 9 (FN): active TB false negatives per 100k — expect 0–900
stopifnot(
  "Excel col 1 (strategy names): first entry should be a string" =
    is.character(excel_ic[[1]][1]),
  "Excel col 4 (TP): values should be 100-1000 (active TB true positives per 100k)" =
    all(na.omit(as.numeric(excel_ic[[4]])) > 50 &
        na.omit(as.numeric(excel_ic[[4]])) < 1500),
  "Excel col 6 (Tot_costs): values should be >100,000 (programme costs per 100k)" =
    median(na.omit(as.numeric(excel_ic[[6]]))) > 100000,
  "Excel col 9 (FN): values should be 0-900 (active TB false negatives per 100k)" =
    all(na.omit(as.numeric(excel_ic[[9]])) >= 0 &
        na.omit(as.numeric(excel_ic[[9]])) < 1000)
)
cat("Excel column validation passed.\n")

get_ic <- function(row_name) {
  r <- excel_ic[excel_ic[[1]] == row_name, ]
  list(
    TP         = as.numeric(r[[4]]),
    LTBI       = as.numeric(r[[5]]),
    Tot_costs  = as.numeric(r[[6]]),
    FN         = as.numeric(r[[9]])
  )
}


# Project-wide colour palette — Lancet scheme throughout
project_pal <- colorRampPalette(c(
  "#00468b", "#0099b4", "#42b540", "#925e9f",
  "#ed0000", "#ad002a", "#fdaf91", "#adb6b6"
))

# Build named colour vector for any set of strategy names
build_strategy_colours <- function(nms) {
  n <- length(nms)
  cols <- project_pal(n)
  setNames(cols, sort(nms))
}

# Clean Excel row names -> human-readable labels for plots
clean_excel_name <- function(nm) {
  known <- c(
    "anycough_cxr(TB)"           = "Cough+CXR (TB sx)",
    "anysx_cxr(any)"             = "Symptom screen+CXR",
    "cough/CXR (parallel)_xpert" = "CXR+Xpert",
    "cough/CXR (parallel)_ultra" = "CXR+Ultra",
    "qft_xpert"                  = "QFT-GIT+Xpert",
    "tspt_xpert"                 = "T-SPOT.TB+Xpert",
    "tst_cxr(TB)"                = "CXR+TST",
    "qft_cxr(TB)"                = "CXR+QFT-GIT",
    "parallel allsx qft_xpert"   = "Parallel Sx+QFT (Xpert)",
    "parallel allsx qft_ultra"   = "Parallel Sx+QFT (Ultra)",
    "parallel cough qft_xpert"   = "Parallel Cough+QFT (Xpert)",
    "parallel cough qft_ultra"   = "Parallel Cough+QFT (Ultra)",
    "parallel allsx tspt_xpert"  = "Parallel Sx+T-SPOT (Xpert)",
    "parallel allsx tspt_ultra"  = "Parallel Sx+T-SPOT (Ultra)",
    "parallel cough tspt_xpert"  = "Parallel Cough+T-SPOT (Xpert)",
    "parallel cough tspt_ultra"  = "Parallel Cough+T-SPOT (Ultra)"
  )
  if (nm %in% names(known)) return(known[[nm]])
  out <- nm
  out <- gsub("cough/CXR \\(parallel\\)", "CXR", out)
  out <- gsub("^anycough", "Cough", out)
  out <- gsub("^anysx",    "Sx screen", out)
  out <- gsub("_", " + ", out)
  out <- gsub("\\btspt\\b|\\btspot\\b", "T-SPOT.TB", out, perl = TRUE)
  out <- gsub("\\bqft\\b",   "QFT-GIT",  out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\btst\\b",   "TST",      out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bcxr\\b",   "CXR",      out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bxpert\\b", "Xpert",    out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\bultra\\b", "Ultra",    out, perl = TRUE, ignore.case = TRUE)
  out <- gsub("\\(TB\\)",    "(TB sx)",  out)
  out <- gsub("([[:alpha:]])\\(", "\\1 (", out)  # ensure space before (
  out <- gsub("\\s*\\+\\s*", "+", out)
  trimws(out)
}

# Assign test family by primary (first-line) test, not confirmatory test.
# Order matters: check IGRA/TST families before Xpert/Ultra so that
# e.g. "QFT-GIT+Xpert" is classified as QFT-GIT, not Molecular.
assign_family <- function(nm) {
  nm_lc <- tolower(nm)
  if (grepl("qft",         nm_lc))  return("QFT-GIT")
  if (grepl("tspot|tspt",  nm_lc))  return("T-SPOT.TB")
  if (grepl("^tst|_tst",   nm_lc))  return("TST Mantoux")
  if (grepl("xpert|ultra", nm_lc))  return("Molecular (Xpert/Ultra)")
  return("CXR / Symptom screen")
}

# eff_active key map: maps Excel row names to the eff_active_* config parameter
# that governs active TB detection uncertainty for that test pathway.
eff_active_key_map <- c(
  "anycough_cxr(any)"                       = "eff_active_CXR_only",
  "2 wk cough_cxr (any)"                    = "eff_active_CXR_only",
  "2 wk cough_cxr (TB)"                     = "eff_active_CXR_only",
  "anycough_cxr(TB)"                        = "eff_active_CXR_only",
  "anysx_cxr(any)"                          = "eff_active_CXR_only",
  "anysx_cxr(TB)"                           = "eff_active_CXR_only",
  "cough/CXR (parallel)_xpert"              = "eff_active_CXR_Xpert",
  "cough/CXR (parallel)_ultra"              = "eff_active_CXR_Xpert",
  "anycough_cxr(any)_xpert"                 = "eff_active_CXR_Xpert",
  "anycough_cxr(any)_ultra"                 = "eff_active_CXR_Xpert",
  "qft_cxr(any)"                            = "eff_active_QFT",
  "qft_cxr(TB)"                             = "eff_active_QFT",
  "qft_xpert"                               = "eff_active_QFT",
  "qft_ultra"                               = "eff_active_QFT",
  "qft_cxr(TB)_xpert"                       = "eff_active_CXR_QFT",
  "qft_cxr(TB)_ultra"                       = "eff_active_CXR_QFT",
  "qft_cough/CXR (parallel)_xpert"          = "eff_active_CXR_QFT",
  "qft_cough/CXR (parallel)_ultra"          = "eff_active_CXR_QFT",
  "tspt_cxr(any)"                           = "eff_active_TSPOT",
  "tspt_cxr(TB)"                            = "eff_active_TSPOT",
  "tspt_xpert"                              = "eff_active_TSPOT",
  "tspt_ultra"                              = "eff_active_TSPOT",
  "tspot_cxr(TB)_xpert"                     = "eff_active_TSPOT",
  "tspot_cxr(TB)_ultra"                     = "eff_active_TSPOT",
  "tspot_cough/CXR (parallel)_xpert"        = "eff_active_TSPOT",
  "tspot_cough/CXR (parallel)_ultra"        = "eff_active_TSPOT",
  "tst_cxr(any)"                            = "eff_active_TST",
  "tst_cxr(TB)"                             = "eff_active_TST",
  "tst_xpert"                               = "eff_active_TST",
  "tst_ultra"                               = "eff_active_TST",
  "tst_cxr(TB)_xpert"                       = "eff_active_TST_QFT",
  "tst_cxr(TB)_ultra"                       = "eff_active_TST_QFT",
  "tst_cough/CXR (parallel)_xpert"          = "eff_active_TST_QFT",
  "tst_cough/CXR (parallel)_ultra"          = "eff_active_TST_QFT",
  "parallel allsx qft_xpert"                = "eff_active_QFT_allsx_xpert",
  "parallel allsx qft_ultra"                = "eff_active_QFT_allsx_ultra",
  "parallel cough qft_xpert"                = "eff_active_QFT_cough_xpert",
  "parallel cough qft_ultra"                = "eff_active_QFT_cough_ultra",
  "parallel allsx tspt_xpert"               = "eff_active_TSPOT",
  "parallel allsx tspt_ultra"               = "eff_active_TSPOT",
  "parallel cough tspt_xpert"               = "eff_active_TSPOT",
  "parallel cough tspt_ultra"               = "eff_active_TSPOT"
)

create_diagnostic_strategies <- function() {
  make_init <- function() setNames(rep(0, length(v_state_names)), v_state_names)

  # Total LTBI and active TB in model cohort (from prevalence assumptions)
  n_ltbi   <- n_c * prev_ltbi    # Berrocal-Almanza 2022: 17.8% → 17,800 per 100k
  n_active <- n_c * prev_active  # 1% → 1,000 per 100k

  # Helper: build init vector from Excel TP/FN/LTBI columns directly
  make_from_excel <- function(ic) {
    init <- make_init()
    ltbi_detected               <- ic$LTBI * (n_c / excel_denom)
    init["ActiveDiagnosed"]     <- ic$TP * (n_c / excel_denom)
    init["ActiveUndiagnosed"]   <- ic$FN * (n_c / excel_denom)
    init["LatentDiagnosed"]     <- ltbi_detected
    init["LatentUndiagnosed"]   <- n_ltbi - ltbi_detected
    init["Uninfected"]          <- n_c * prev_uninfected
    init
  }

  # ---------------------------------------------------------------------------
  # All 34 feasible pathways from Zenner et al. 2025 ERJ decision tree
  # (rows with non-NA Tot_costs) plus Passive_case_finding reference = 35 total.
  # 8 parallel IGRA strategies with Option B LTBI override added below = 43 total.
  # ---------------------------------------------------------------------------

  # ---- Reference: Passive case finding ----
  s <- list(name = "Passive case finding", excel_name = NA_character_,
            init = make_init(), test_cost = 0)
  s$init["Uninfected"]        <- n_c * prev_uninfected
  s$init["LatentUndiagnosed"] <- n_c * prev_ltbi
  s$init["ActiveUndiagnosed"] <- n_c * prev_active
  strategies <- list(Passive_case_finding = s)
  strategies[["Passive_case_finding"]]$eff_active_key <- NA_character_
  strategies[["Passive_case_finding"]]$base_tp        <- 0
  strategies[["Passive_case_finding"]]$base_fn        <- n_active  # all active TB undetected

  # ---- All Excel rows with non-NA Tot_costs ----
  for (i in seq_len(nrow(excel_ic))) {
    row      <- excel_ic[i, ]
    excel_nm <- as.character(row[[1]])
    if (is.na(excel_nm) || excel_nm == "") next
    tc <- suppressWarnings(as.numeric(row[[6]]))
    if (is.na(tc)) next
    ic <- list(
      TP        = as.numeric(row[[4]]),
      LTBI      = as.numeric(row[[5]]),
      Tot_costs = tc,
      FN        = as.numeric(row[[9]])
    )
    key <- gsub("[^A-Za-z0-9]", "_", excel_nm)
    strategies[[key]] <- list(
      name       = clean_excel_name(excel_nm),
      excel_name = excel_nm,
      init       = make_from_excel(ic),
      test_cost  = ic$Tot_costs / excel_denom
    )
    strategies[[key]]$eff_active_key <- unname(eff_active_key_map[excel_nm])
    strategies[[key]]$base_tp        <- ic$TP * (n_c / excel_denom)
    strategies[[key]]$base_fn        <- ic$FN * (n_c / excel_denom)
  }
  # ---------------------------------------------------------------------------
  # OPTION B: 8 parallel IGRA strategies — LTBI recalculated from
  # population prevalence × pooled IGRA sensitivity (Zenner 2025 Table 1).
  # Excel LTBI values for these rows reflect sequential pathway yield and are
  # overridden here with universal-IGRA estimates.
  #
  # LTBI override: n_ltbi × igra_sensitivity
  #   QFT pooled sensitivity:     0.83 (Zenner 2025 Table 1)
  #   T-SPOT.TB sensitivity:      0.88 (Zenner 2025 Table 1)
  #
  # Cost estimation: existing IGRA+confirmatory base cost + IGRA testing programme cost
  #   allsx (any TB symptom) screen: cost_igra_programme_allsx per 100k (base £500,000; £5/pp)
  #   cough screen only:             cost_igra_programme_cough per 100k (base £300,000; £3/pp)
  #
  # Programme delivery cost values are read from paramsData (config.csv rows
  # cost_igra_programme_cough / cost_igra_programme_allsx) so they can be varied in PSA and DSA.
  # ---------------------------------------------------------------------------

  prog_cost_cough <- as.numeric(config_list["cost_igra_programme_cough"])  # £ per cohort of 100k
  prog_cost_allsx <- as.numeric(config_list["cost_igra_programme_allsx"])

  # Net LTBI detected at 100% uptake = TP − FP
  # FP = (1 - specificity) × n_uninfected; source: Pai et al. 2008 BCG-vaccinated meta-analysis
  # QFT: 17,800 × 0.83 − 0.04 × 81,200 = 14,774 − 3,248 = 11,526/100k
  # T-SPOT: 17,800 × 0.88 − 0.07 × 81,200 = 15,664 − 5,684 = 9,980/100k
  igra_ltbi_qft_100pct  <- max(0, n_ltbi * igra_sens_qft  - (1 - igra_spec_qft)  * n_c * prev_uninfected)
  igra_ltbi_tspt_100pct <- max(0, n_ltbi * igra_sens_tspt - (1 - igra_spec_tspt) * n_c * prev_uninfected)
  # Apply base-case uptake (default 75%; real-world 60-80%, UKHSA 2024)
  igra_ltbi_qft  <- igra_ltbi_qft_100pct  * igra_uptake_base
  igra_ltbi_tspt <- igra_ltbi_tspt_100pct * igra_uptake_base

  # Read TP and FN for parallel IGRA strategies from Excel rows.
  # These rows (parallel allsx/cough qft/tspt xpert/ultra) have TP/FN populated
  # in the decision tree but Tot_costs = NA (hence excluded from the main loop).
  # If a row is missing or has NA values, falls back to Zenner 2025 Table 1 values.
  lookup_parallel_tp_fn <- function(nm, tp_fallback, fn_fallback) {
    r <- excel_ic[!is.na(excel_ic[[1]]) & excel_ic[[1]] == nm, ]
    if (nrow(r) == 0 || is.na(as.numeric(r[[4]])))  {
      warning(paste0("Excel row '", nm, "' not found or TP is NA — using Zenner 2025 Table 1 fallback: TP=", tp_fallback))
      return(list(tp = tp_fallback, fn = fn_fallback))
    }
    list(tp = as.numeric(r[[4]]), fn = as.numeric(r[[9]]))
  }

  make_parallel_igra_init <- function(tp, fn, ltbi_detected) {
    init <- make_init()
    init["ActiveDiagnosed"]   <- tp * (n_c / excel_denom)
    init["ActiveUndiagnosed"] <- fn * (n_c / excel_denom)
    init["LatentDiagnosed"]   <- ltbi_detected
    init["LatentUndiagnosed"] <- n_ltbi - ltbi_detected
    init["Uninfected"]        <- n_c * prev_uninfected
    init
  }

  qft_xpert_base  <- as.numeric(excel_ic[excel_ic[[1]] == "qft_xpert",  ][[6]])
  qft_ultra_base  <- as.numeric(excel_ic[excel_ic[[1]] == "qft_ultra",  ][[6]])
  tspt_xpert_base <- as.numeric(excel_ic[excel_ic[[1]] == "tspt_xpert", ][[6]])
  tspt_ultra_base <- as.numeric(excel_ic[excel_ic[[1]] == "tspt_ultra", ][[6]])

  parallel_igra_specs <- list(
    list(nm = "parallel allsx qft_xpert",  pf = lookup_parallel_tp_fn("parallel allsx qft_xpert",  836.6, 163.4), ltbi = igra_ltbi_qft,  ltbi_100pct = igra_ltbi_qft_100pct,  cost = qft_xpert_base  + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel allsx qft_ultra",  pf = lookup_parallel_tp_fn("parallel allsx qft_ultra",  874.6, 125.4), ltbi = igra_ltbi_qft,  ltbi_100pct = igra_ltbi_qft_100pct,  cost = qft_ultra_base  + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel cough qft_xpert",  pf = lookup_parallel_tp_fn("parallel cough qft_xpert",  791.7, 208.3), ltbi = igra_ltbi_qft,  ltbi_100pct = igra_ltbi_qft_100pct,  cost = qft_xpert_base  + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel cough qft_ultra",  pf = lookup_parallel_tp_fn("parallel cough qft_ultra",  827.7, 172.3), ltbi = igra_ltbi_qft,  ltbi_100pct = igra_ltbi_qft_100pct,  cost = qft_ultra_base  + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel allsx tspt_xpert", pf = lookup_parallel_tp_fn("parallel allsx tspt_xpert", 849.4, 150.6), ltbi = igra_ltbi_tspt, ltbi_100pct = igra_ltbi_tspt_100pct, cost = tspt_xpert_base + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel allsx tspt_ultra", pf = lookup_parallel_tp_fn("parallel allsx tspt_ultra", 888.0, 112.0), ltbi = igra_ltbi_tspt, ltbi_100pct = igra_ltbi_tspt_100pct, cost = tspt_ultra_base + prog_cost_allsx, prog_cost = prog_cost_allsx, prog_cost_type = "allsx"),
    list(nm = "parallel cough tspt_xpert", pf = lookup_parallel_tp_fn("parallel cough tspt_xpert", 817.7, 182.3), ltbi = igra_ltbi_tspt, ltbi_100pct = igra_ltbi_tspt_100pct, cost = tspt_xpert_base + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough"),
    list(nm = "parallel cough tspt_ultra", pf = lookup_parallel_tp_fn("parallel cough tspt_ultra", 854.9, 145.1), ltbi = igra_ltbi_tspt, ltbi_100pct = igra_ltbi_tspt_100pct, cost = tspt_ultra_base + prog_cost_cough, prog_cost = prog_cost_cough, prog_cost_type = "cough")
  )

  for (ps in parallel_igra_specs) {
    key <- gsub("[^A-Za-z0-9]", "_", ps$nm)
    ps_eff_key <- unname(eff_active_key_map[ps$nm])
    if (is.na(ps_eff_key)) ps_eff_key <- if (grepl("tspt|tspot", ps$nm, ignore.case = TRUE)) "eff_active_TSPOT" else "eff_active_QFT"
    strategies[[key]] <- list(
      name            = clean_excel_name(ps$nm),
      excel_name      = ps$nm,
      init            = make_parallel_igra_init(ps$pf$tp, ps$pf$fn, ps$ltbi),
      test_cost       = ps$cost / excel_denom,
      prog_cost_total = ps$prog_cost,
      prog_cost_type  = ps$prog_cost_type,
      eff_active_key  = ps_eff_key,
      base_tp         = ps$pf$tp * (n_c / excel_denom),
      base_fn         = ps$pf$fn * (n_c / excel_denom),
      ltbi_100pct     = ps$ltbi_100pct  # LTBI detected at 100% uptake — used by uptake SA
    )
  }

  cat(sprintf("Parallel IGRA strategies added (Option B): LTBI detection = %.0f (QFT) / %.0f (T-SPOT) per 100k at %.0f%% uptake\n",
              igra_ltbi_qft, igra_ltbi_tspt, igra_uptake_base * 100))

  return(strategies)
}

# -----------------------------------------------------------------------------
# Run the Markov model for a single diagnostic strategy
#
# Passes the strategy-specific initial state distribution to model() and
# adds the one-time screening test cost to the total (test costs are incurred
# at entry and are not part of the recurring state-based cost structure).
# -----------------------------------------------------------------------------
run_strategy <- function(strategy, params = paramsData, n_t_override = NULL) {
  result <- model(params, init_dist = strategy$init, n_t_override = n_t_override)
  result$total_cost    <- result$total_cost + (strategy$test_cost * n_c)
  result$test_cost     <- strategy$test_cost * n_c
  result$strategy_name <- strategy$name
  return(result)
}

# -----------------------------------------------------------------------------
# Incremental cost-effectiveness ratio (ICER) calculation
#
# Performs a full sequential incremental analysis following NICE methods:
#   1. Pairwise ICERs vs No Screening (retained for reference).
#   2. Simple dominance: any strategy with higher cost AND lower QALYs than a
#      cheaper alternative is flagged "simply dominated".
#   3. Extended dominance: strategies whose sequential ICER exceeds that of the
#      next comparator on the efficient frontier are iteratively removed.
#   4. Sequential ICERs on the efficient frontier are reported.
#
# The NICE willingness-to-pay threshold is £25,000–£35,000 per QALY (updated Dec 2025, effective Apr 2026).
# -----------------------------------------------------------------------------
calculate_icer <- function(strategy_results) {
  df <- tibble(
    strategy = sapply(strategy_results, function(x) x$strategy_name),
    cost     = sapply(strategy_results, function(x) x$total_cost),
    qaly     = sapply(strategy_results, function(x) x$total_qaly),
    deaths   = sapply(strategy_results, function(x) x$state_membership[nrow(x$state_membership), "Dead"])
  )

  df <- df %>% arrange(cost)

  ref_idx <- which(df$strategy == "Passive case finding")
  if (length(ref_idx) == 0) ref_idx <- 1

  # Pairwise ICERs vs No Screening (retained for backward compatibility)
  df <- df %>%
    mutate(
      inc_cost        = cost - cost[ref_idx],
      inc_qaly        = qaly - qaly[ref_idx],
      icer            = ifelse(inc_qaly != 0, inc_cost / inc_qaly, NA),
      cost_per_person = cost / n_c,
      qaly_per_person = qaly / n_c
    )

  # -------------------------------------------------------------------
  # Sequential incremental analysis with extended dominance detection
  # -------------------------------------------------------------------
  df$dominance      <- "non-dominated"
  df$sequential_icer <- NA_real_
  df$dominance[ref_idx] <- "ref"

  # Step 1: simple dominance — flag any strategy where a cheaper alternative
  #         achieves at least as many QALYs.
  for (i in seq_len(nrow(df))) {
    if (df$dominance[i] == "ref") next
    cheaper <- df[df$cost < df$cost[i], ]
    if (nrow(cheaper) > 0 && any(cheaper$qaly >= df$qaly[i])) {
      df$dominance[i] <- "simply dominated"
    }
  }

  # Step 2: extended dominance — iteratively remove strategies whose sequential
  #         ICER exceeds that of the next comparator (frontier not convex).
  changed <- TRUE
  while (changed) {
    changed  <- FALSE
    frontier <- df[df$dominance %in% c("ref", "non-dominated"), ] %>% arrange(cost)
    n_f      <- nrow(frontier)
    if (n_f < 3) break  # need at least 3 points to identify extended dominance
    for (i in 2:(n_f - 1)) {
      d_cost_i    <- frontier$cost[i]   - frontier$cost[i - 1]
      d_qaly_i    <- frontier$qaly[i]   - frontier$qaly[i - 1]
      d_cost_next <- frontier$cost[i + 1] - frontier$cost[i]
      d_qaly_next <- frontier$qaly[i + 1] - frontier$qaly[i]
      if (d_qaly_i <= 0 || d_qaly_next <= 0) next
      seq_icer_i    <- d_cost_i    / d_qaly_i
      seq_icer_next <- d_cost_next / d_qaly_next
      if (seq_icer_i > seq_icer_next) {
        idx <- which(df$strategy == frontier$strategy[i])
        df$dominance[idx] <- "extendedly dominated"
        changed <- TRUE
        break  # restart after each removal
      }
    }
  }

  # Step 3: assign sequential ICERs on the final efficient frontier.
  frontier_final <- df[df$dominance %in% c("ref", "non-dominated"), ] %>% arrange(cost)
  for (i in seq_len(nrow(frontier_final))) {
    idx <- which(df$strategy == frontier_final$strategy[i])
    if (i == 1) {
      df$sequential_icer[idx] <- NA_real_
    } else {
      d_cost <- frontier_final$cost[i] - frontier_final$cost[i - 1]
      d_qaly <- frontier_final$qaly[i] - frontier_final$qaly[i - 1]
      df$sequential_icer[idx] <- if (d_qaly != 0) d_cost / d_qaly else NA_real_
    }
  }

  return(df)
}

# -----------------------------------------------------------------------------
# Probabilistic sensitivity analysis (PSA)
#
# Runs n_sim Monte Carlo simulations. In each simulation:
#   1. Markov parameters (costs, QALYs, transition probabilities) are drawn
#      from their distributions via sample_params().
#   2. All 43 strategies are run with fixed initial conditions and their
#      costs/QALYs recorded.
#
# Note: initial state distributions are fixed from the Excel decision tree outputs.
# eff_active_* parameters are re-sampled per simulation via beta distributions in config.csv;
# background LTBI prevalence is held fixed (addressed in scenario analyses).
# -----------------------------------------------------------------------------
run_psa <- function(strategies_base, n_sim = 1000, n_cores = NULL) {
  # Samples Markov parameters and eff_active_* detection rates from config.csv distributions.
  # Initial LTBI state distributions are fixed; their uncertainty is addressed in DSA.
  #
  # Parallelisation: uses parallel::mclapply (fork-based; works on macOS/Linux).
  # Each simulation gets a pre-assigned seed derived from the master seed so results
  # are reproducible regardless of core count. n_cores defaults to (detected - 1).

  if (is.null(n_cores)) n_cores <- max(1L, detectCores() - 1L)
  cat(sprintf("\nRunning PSA with %d simulations (%d cores)...\n", n_sim, n_cores))

  # Pre-generate per-simulation seeds from the current RNG state so parallelism
  # does not break reproducibility — same seeds regardless of core count.
  sim_seeds <- sample.int(.Machine$integer.max, n_sim)

  one_sim <- function(i) {
    set.seed(sim_seeds[i])
    sampled <- sample_params()

    # Sample active TB detection efficiency parameters for this PSA simulation.
    # Each eff_active_* parameter has a beta distribution in config.csv.
    eff_sampled <- setNames(
      rbeta(length(eff_active_param_names),
            shape1 = eff_active_dist_p1[eff_active_param_names],
            shape2 = eff_active_dist_p2[eff_active_param_names]),
      eff_active_param_names
    )

    rows    <- vector("list", length(strategies_base))
    for (k in seq_along(strategies_base)) {
      s  <- strategies_base[[k]]
      ek <- s$eff_active_key
      if (!is.null(ek) && !is.na(ek)) {
        # Rescale ActiveDiagnosed/ActiveUndiagnosed using the sampled detection efficiency.
        # Total active TB at entry = base_tp + base_fn (fixed from Excel/prevalence).
        new_eff <- eff_sampled[[ek]]
        s$init["ActiveDiagnosed"]   <- n_c * prev_active * new_eff
        s$init["ActiveUndiagnosed"] <- n_c * prev_active * (1 - new_eff)
      }
      result <- run_strategy(s, params = sampled)
      rows[[k]] <- tibble(
        sim      = i,
        strategy = strategies_base[[k]]$name,
        cost     = result$total_cost,
        qaly     = result$total_qaly
      )
    }
    rows
  }

  raw <- mclapply(seq_len(n_sim), one_sim, mc.cores = n_cores)

  psa_df <- bind_rows(lapply(raw, bind_rows))
  return(psa_df)
}

#------------------------------------------------------------------------------#
# PSA Visualisation: Cost-Effectiveness Plane
#------------------------------------------------------------------------------#
plot_ce_plane <- function(psa_df, reference = "Passive case finding", wtp = 20000) {
  ref_data <- psa_df %>% filter(strategy == reference)

  ce_df <- psa_df %>%
    filter(strategy != reference) %>%
    left_join(ref_data %>% select(sim, ref_cost = cost, ref_qaly = qaly), by = "sim") %>%
    mutate(
      inc_cost = cost - ref_cost,
      inc_qaly = qaly - ref_qaly
    )

  p <- ggplot(ce_df, aes(x = inc_qaly, y = inc_cost, color = strategy)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(intercept = 0, slope = wtp, linetype = "dotted", color = "#ad002a", linewidth = 1) +
    annotate("text", x = max(ce_df$inc_qaly) * 0.8, y = wtp * max(ce_df$inc_qaly) * 0.8,
             label = paste0("WTP = £", format(wtp, big.mark = ",")), color = "#ad002a", size = 3.5) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_y_continuous(labels = function(x) paste0("\u00a3", format(round(x), big.mark = ",", scientific = FALSE))) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Incremental QALYs",
      y = "Incremental Cost (\u00a3)",
      title = NULL,
      subtitle = paste0("vs ", reference, " | WTP threshold = £", format(wtp, big.mark = ","), "/QALY"),
      color = "Strategy"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  return(p)
}

#------------------------------------------------------------------------------#
# PSA Visualisation: Cost-Effectiveness Acceptability Curve (CEAC)
#------------------------------------------------------------------------------#
plot_ceac <- function(psa_df, reference = "Passive case finding",
                      wtp_range = seq(0, 50000, by = 1000)) {
  ref_data <- psa_df %>% filter(strategy == reference)
  n_sim <- max(psa_df$sim)

  strategies_to_compare <- unique(psa_df$strategy)

  ceac_data <- expand.grid(wtp = wtp_range, strategy = strategies_to_compare,
                           stringsAsFactors = FALSE)
  ceac_data$prob_ce <- NA

  for (w in wtp_range) {
    # NMB = QALY × WTP − cost; strategy with highest NMB is optimal
    nmb_df <- psa_df %>%
      mutate(nmb = qaly * w - cost) %>%
      group_by(sim) %>%
      mutate(is_best = nmb == max(nmb)) %>%
      ungroup()

    # Probability of being cost-effective at this WTP
    prob_df <- nmb_df %>%
      group_by(strategy) %>%
      summarise(prob_ce = mean(is_best), .groups = "drop")

    for (s in strategies_to_compare) {
      idx <- which(ceac_data$wtp == w & ceac_data$strategy == s)
      prob_val <- prob_df$prob_ce[prob_df$strategy == s]
      if (length(prob_val) > 0) ceac_data$prob_ce[idx] <- prob_val
    }
  }

  # Distinct colours per frontier strategy for CEAC readability.
  # Canonical Lancet palette: gray / teal / dark red (Symptom screen+CXR excluded — 0% PSA-optimal).
  ceac_strat_cols <- c(
    "Passive case finding"    = "#adb6b6",
    "Cough+CXR (TB sx)"       = "#0099b4",
    "Parallel Sx+QFT (Ultra)" = "#ad002a"
  )
  # Fall back to auto palette for any strategies not in the named set
  all_ceac_strats <- unique(ceac_data$strategy)
  extra_strats    <- setdiff(all_ceac_strats, names(ceac_strat_cols))
  if (length(extra_strats) > 0) {
    extra_cols <- setNames(project_pal(length(extra_strats)), extra_strats)
    ceac_strat_cols <- c(ceac_strat_cols, extra_cols)
  }

  p <- ggplot(ceac_data, aes(x = wtp, y = prob_ce, color = strategy)) +
    geom_line(linewidth = 1.4) +
    geom_vline(xintercept = 25000, linetype = "dashed", color = "grey40") +
    annotate("text", x = 25000, y = 0.05, label = "NICE\n£25k", color = "grey40", size = 3, hjust = -0.1) +
    geom_vline(xintercept = 35000, linetype = "dashed", color = "grey60") +
    annotate("text", x = 35000, y = 0.05, label = "NICE\n£35k", color = "grey60", size = 3, hjust = -0.1) +
    scale_x_continuous(labels = function(x) paste0("\u00a3", format(x, big.mark = ",", scientific = FALSE))) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_colour_manual(values = ceac_strat_cols) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Willingness-to-Pay Threshold (£/QALY)",
      y = "Probability Cost-Effective",
      title = NULL,
      subtitle = "Probability each strategy is optimal at given WTP thresholds",
      color = "Strategy"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  return(p)
}
