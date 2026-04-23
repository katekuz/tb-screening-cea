#!/usr/bin/env Rscript
# =============================================================================
# verify_project.R — TB Screening CEA Project Verification
# Run: Rscript verify_project.R
# =============================================================================

suppressMessages({ library(dplyr) })

# Counters and helpers --------------------------------------------------------
n_fails <- 0L; n_warns <- 0L

PASS <- function(msg) cat(sprintf("  \x1b[32m[PASS]\x1b[0m  %s\n", msg))
FAIL <- function(msg) { cat(sprintf("  \x1b[31m[FAIL]\x1b[0m  %s\n", msg)); n_fails <<- n_fails + 1L }
WARN <- function(msg) { cat(sprintf("  \x1b[33m[WARN]\x1b[0m  %s\n", msg)); n_warns <<- n_warns + 1L }
HEAD <- function(msg) cat(sprintf("\n\x1b[1m%s\x1b[0m\n%s\n", msg, strrep("-", nchar(msg))))

ok   <- function(cond, pass_msg, fail_msg) if (cond) PASS(pass_msg) else FAIL(fail_msg)
okw  <- function(cond, pass_msg, warn_msg) if (cond) PASS(pass_msg) else WARN(warn_msg)

# Load data -------------------------------------------------------------------
bc      <- read.csv("output/csv/icer_nhs_perspective.csv",             stringsAsFactors = FALSE)

tx      <- read.csv("output/csv/icer_societal_perspective.csv",         stringsAsFactors = FALSE)
cfg     <- read.csv("input/config.csv",                                 stringsAsFactors = FALSE)
psa     <- read.csv("output/csv/psa_results.csv",                     stringsAsFactors = FALSE)
ltbi_sa <- read.csv("output/csv/ltbi_completion_sensitivity.csv",     stringsAsFactors = FALSE)
prev_sa <- read.csv("output/csv/ltbi_prevalence_sensitivity.csv",     stringsAsFactors = FALSE)

tx_counts <- read.csv("output/csv/transmission_prevention_counts.csv",  stringsAsFactors = FALSE)

cfg_val <- function(p) { v <- cfg$value[cfg$parameter == p]; if (length(v)) suppressWarnings(as.numeric(v[1])) else NA_real_ }

# =============================================================================
HEAD("1. STRUCTURAL INTEGRITY")
# =============================================================================

n_strat <- nrow(bc)
ok(n_strat == 43,
   sprintf("43 strategies in base-case table (%d rows)", n_strat),
   sprintf("Expected 43 strategies, found %d", n_strat))

n_c_cfg     <- cfg_val("n_c")
n_c_implied <- round(bc$cost[1] / bc$cost_per_person[1])
ok(abs(n_c_implied - n_c_cfg) < 1,
   sprintf("Cohort size consistent: n_c = %s", format(n_c_cfg, big.mark = ",")),
   sprintf("Cohort size mismatch: config=%s vs implied=%s", n_c_cfg, n_c_implied))

n_t <- cfg_val("n_t")
okw(n_t == 660,
    sprintf("Time horizon: %d cycles (55 yr lifetime)", n_t),
    sprintf("Time horizon is %d cycles (expected 660)", n_t))

disc <- cfg_val("discount")
ok(abs(disc - 0.035) < 1e-6,
   "Discount rate: 3.5% (NICE reference case)",
   sprintf("Discount rate: %.4f (expected 0.035)", disc))

ref_row <- bc[bc$dominance == "ref", ]
ok(nrow(ref_row) == 1 && ref_row$strategy == "Passive case finding",
   "Reference strategy: 'Passive case finding' (single ref row)",
   "Reference strategy not 'Passive case finding' or multiple ref rows")

ok(ref_row$ltbi_detected == 0 && ref_row$active_detected == 0,
   "Passive case finding: 0 LTBI detected, 0 active detected",
   sprintf("Passive case finding LTBI=%s active=%s (expected 0/0)", ref_row$ltbi_detected, ref_row$active_detected))

no_scr_qaly <- bc[bc$strategy == "Passive case finding", "qaly_per_person"]
okw(no_scr_qaly > 21.5 && no_scr_qaly < 21.7,
    sprintf("Passive case finding QALY/person: %.4f (expected 21.5\u201321.7)", no_scr_qaly),
    sprintf("Passive case finding QALY = %.4f \u2014 outside expected range 21.5\u201321.7", no_scr_qaly))

ok(all(bc$qaly_per_person >= no_scr_qaly),
   "All strategies have QALYs \u2265 Passive case finding (no negative-QALY interventions)",
   "Some strategies have lower QALYs than Passive case finding \u2014 check transition matrix")

# =============================================================================
HEAD("2. EFFICIENT FRONTIER (BASE CASE)")
# =============================================================================

frontier <- bc[bc$dominance %in% c("ref", "non-dominated"), ]
frontier <- frontier[order(frontier$cost), ]
f_names  <- frontier$strategy
expected_frontier <- c("Passive case finding", "Cough+CXR (TB sx)", "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)")

ok(length(f_names) == 4,
   sprintf("Frontier has 4 strategies: %s", paste(f_names, collapse = " \u2192 ")),
   sprintf("Expected 4 frontier strategies, found %d: %s", length(f_names), paste(f_names, collapse = ", ")))

ok(identical(f_names, expected_frontier),
   "Frontier strategies: Passive case finding \u2192 Cough+CXR \u2192 Symptom screen+CXR \u2192 Parallel Sx+QFT (Ultra)",
   sprintf("Frontier mismatch. Expected: %s | Got: %s",
           paste(expected_frontier, collapse=", "), paste(f_names, collapse=", ")))

ok(all(diff(frontier$cost) > 0) && all(diff(frontier$qaly) > 0),
   "Frontier is strictly cost- and QALY-increasing (no reversals)",
   "Frontier cost or QALY ordering is not strictly increasing")

cqft <- bc[bc$strategy == "Parallel Cough+QFT (Ultra)", ]
ok(nrow(cqft) == 1 && cqft$dominance == "extendedly dominated",
   "Parallel Cough+QFT (Ultra): extendedly dominated (correct post-ISSUE-18)",
   sprintf("Parallel Cough+QFT dominance = '%s' (expected 'extendedly dominated')",
           if (nrow(cqft) == 1) cqft$dominance else "not found"))

tspt_dom <- bc[grepl("T-SPOT", bc$strategy), "dominance"]
ok(length(tspt_dom) > 0 && all(tspt_dom %in% c("simply dominated", "extendedly dominated")),
   sprintf("All %d T-SPOT strategies dominated (simply or extendedly; lower spec \u2192 more FP)", length(tspt_dom)),
   sprintf("%d T-SPOT strategies on frontier (unexpected)", sum(!tspt_dom %in% c("simply dominated","extendedly dominated"))))

seq_icers <- frontier$sequential_icer[!is.na(frontier$sequential_icer)]
ok(all(seq_icers < 35000),
   sprintf("All frontier sequential ICERs < \u00a335,000 (NICE WTP range \u00a325k\u2013\u00a335k; range: \u00a3%s\u2013\u00a3%s)",
           format(round(min(seq_icers)), big.mark=","), format(round(max(seq_icers)), big.mark=",")),
   sprintf("Frontier ICER(s) exceed \u00a335k (NICE upper WTP): %s",
           paste(format(round(seq_icers[seq_icers >= 35000]), big.mark=","), collapse=", ")))

par_min  <- min(bc$ltbi_detected[grepl("^Parallel", bc$strategy)])
seq_max  <- max(bc$ltbi_detected[!grepl("^Parallel|^No", bc$strategy)])
okw(par_min > seq_max,
    sprintf("Parallel IGRA LTBI detection (%s min) > sequential max (%s)",
            format(round(par_min), big.mark=","), format(round(seq_max), big.mark=",")),
    sprintf("Parallel IGRA min LTBI (%s) not > sequential max (%s)",
            format(round(par_min), big.mark=","), format(round(seq_max), big.mark=",")))

qft_ltbi  <- bc$ltbi_detected[bc$strategy == "Parallel Cough+QFT (Ultra)"]
tspt_ltbi <- bc$ltbi_detected[bc$strategy == "Parallel Cough+T-SPOT (Ultra)"]
ok(length(qft_ltbi)==1 && length(tspt_ltbi)==1 && qft_ltbi > tspt_ltbi,
   sprintf("QFT detects more LTBI than T-SPOT (%s vs %s; spec 0.96 vs 0.93)",
           format(round(qft_ltbi),big.mark=","), format(round(tspt_ltbi),big.mark=",")),
   "QFT does not detect more LTBI than T-SPOT (specificity correction may be broken)")

# =============================================================================
HEAD("3. PRIMARY ANALYSIS (TRANSMISSION) — LOGICAL CONSISTENCY")
# =============================================================================

# Transmission: Parallel Sx+QFT cost must fall vs base case
cost_bc <- bc[bc$strategy  == "Parallel Sx+QFT (Ultra)", "cost_per_person"]
cost_tx <- tx[tx$strategy  == "Parallel Sx+QFT (Ultra)", "cost_per_person"]
ok(length(cost_tx)==1 && cost_tx < cost_bc,
   sprintf("Primary: Parallel Sx+QFT cost falls \u00a3%s \u2192 \u00a3%s (secondary case savings)",
           round(cost_bc), round(cost_tx)),
   "Primary: Parallel Sx+QFT cost did not fall")

# Transmission frontier: same 4 strategies, Parallel Cough+QFT remains ext. dominated
tx_front <- tx[tx$dominance %in% c("ref", "non-dominated"), "strategy"]
ok(!"Parallel Cough+QFT (Ultra)" %in% tx_front,
   "Primary: Parallel Cough+QFT remains extendedly dominated (confirmed CSV)",
   "Primary: Parallel Cough+QFT unexpectedly entered transmission frontier")

# Transmission frontier: same 4 strategies
tx_front_names <- tx[tx$dominance %in% c("ref", "non-dominated"), "strategy"]
ok(identical(sort(tx_front_names), sort(expected_frontier)),
   sprintf("Primary frontier: same 4 strategies as robustness check"),
   sprintf("Primary frontier differs: %s", paste(tx_front_names, collapse=", ")))

# Societal ICER
tx_qft_icer <- tx[tx$strategy == "Parallel Sx+QFT (Ultra)", "sequential_icer"]
ok(length(tx_qft_icer)==1 && abs(tx_qft_icer - 12738) <= 300,
   sprintf("Parallel Sx+QFT societal ICER: \u00a3%s (expected ~\\u00a312,738)",
           format(round(tx_qft_icer), big.mark=",")),
   sprintf("Parallel Sx+QFT societal ICER = \u00a3%s (expected ~\\u00a312,738)",
           format(round(tx_qft_icer), big.mark=",")))

# Transmission counts
qft_atb <- tx_counts[tx_counts$strategy == "Parallel Sx+QFT (Ultra)", "atb_cases_prevented"]
ok(length(qft_atb)==1 && qft_atb > 0,
   sprintf("ATB cases prevented by Parallel Sx+QFT: %s per 100k", format(round(qft_atb), big.mark=",")),
   "Parallel Sx+QFT ATB prevented is 0 or missing")

# SA checks — unchanged
base_icer_qft <- bc[bc$strategy == "Parallel Sx+QFT (Ultra)", "sequential_icer"]
rw_row <- ltbi_sa[ltbi_sa$scenario == "realworld_55pct" & ltbi_sa$strategy == "Parallel Sx+QFT (Ultra)", ]
okw(nrow(rw_row)==1 && !is.na(rw_row$sequential_icer) && rw_row$sequential_icer > base_icer_qft,
    sprintf("LTBI completion SA: ICER increases at real-world completion (\u00a3%s \u2192 \u00a3%s)",
            format(round(base_icer_qft),big.mark=","),
            format(round(rw_row$sequential_icer),big.mark=",")),
    "LTBI completion SA: ICER direction unexpected or strategy not found")

prev_row <- prev_sa[prev_sa$scenario == "ukhsa2025_15.1pct" & prev_sa$strategy == "Parallel Sx+QFT (Ultra)", ]
okw(nrow(prev_row)==1 && !is.na(prev_row$sequential_icer) && prev_row$sequential_icer < 25000,
    sprintf("LTBI prevalence SA: Parallel Sx+QFT CE at 15.1%% LTBI prevalence (ICER \u00a3%s)",
            format(round(prev_row$sequential_icer),big.mark=",")),
    "LTBI prevalence SA: strategy not CE or not found at 15.1%% LTBI prevalence")

# =============================================================================
HEAD("4. KEY PARAMETER VALUES (config.csv)")
# =============================================================================

param_checks <- list(
  c("p_LatentTreated_LatentCompleted", 0.253,   0.001, "Surey 2021 76%/3mo"),
  c("p_LatentTreated_LatentLtfu",      0.057,   0.001, "Surey 2021 LTFU fix"),
  c("p_ActiveTreated_ActiveLtfu",      0.005,   0.001, "UKHSA 2025 ~3%/6mo (ISSUE-18)"),
  c("p_LatentDiagnosed_ActiveUndiagnosed", 0.000912,1e-5,"Berrocal-Almanza 2022 1.09%/yr reactivation"),
  c("igra_ltbi_sens_qft",              0.83,    0.01,  "Zenner 2025 pooled"),
  c("igra_ltbi_sens_tspt",             0.88,    0.01,  "Zenner 2025 pooled"),
  c("igra_spec_qft",                   0.96,    0.01,  "Pai 2008 BCG"),
  c("igra_spec_tspt",                  0.93,    0.01,  "Pai 2008 BCG"),
  c("p_LatentDiagnosed_LatentTreated", 0.151,   0.001, "UKHSA 2024 74.8% initiation"),
  c("discount",                        0.035,   1e-6,  "NICE 3.5%")
)

for (chk in param_checks) {
  pname <- chk[1]; expected <- as.numeric(chk[2]); tol <- as.numeric(chk[3]); src <- chk[4]
  actual <- cfg_val(pname)
  if (is.na(actual)) {
    FAIL(sprintf("Parameter '%s' not found in config.csv", pname))
  } else {
    ok(abs(actual - expected) <= tol,
       sprintf("%-45s = %.5f  (%s)", pname, actual, src),
       sprintf("%-45s = %.5f (expected %.5f) \u2014 %s", pname, actual, expected, src))
  }
}

beta_params   <- cfg[!is.na(cfg$distribution) & cfg$distribution == "beta", ]
bad_beta      <- beta_params[is.na(suppressWarnings(as.numeric(beta_params$dist_param1))) |
                              suppressWarnings(as.numeric(beta_params$dist_param1)) <= 0 |
                              is.na(suppressWarnings(as.numeric(beta_params$dist_param2))) |
                              suppressWarnings(as.numeric(beta_params$dist_param2)) <= 0, ]
ok(nrow(bad_beta) == 0,
   sprintf("All %d beta-distributed parameters have valid shape params (>0)", nrow(beta_params)),
   sprintf("%d beta parameters have invalid shape params: %s",
           nrow(bad_beta), paste(bad_beta$parameter, collapse=", ")))

gamma_params  <- cfg[!is.na(cfg$distribution) & cfg$distribution == "gamma", ]
bad_gamma     <- gamma_params[is.na(suppressWarnings(as.numeric(gamma_params$dist_param1))) |
                               suppressWarnings(as.numeric(gamma_params$dist_param1)) <= 0, ]
ok(nrow(bad_gamma) == 0,
   sprintf("All %d gamma-distributed cost parameters valid", nrow(gamma_params)),
   sprintf("%d gamma parameters invalid: %s", nrow(bad_gamma), paste(bad_gamma$parameter, collapse=", ")))

# =============================================================================
HEAD("5. PSA CONSISTENCY")
# =============================================================================

n_psa_sims <- length(unique(psa$sim))
ok(nrow(psa) == 43000 && n_psa_sims == 1000,
   sprintf("PSA: 43,000 rows (43 strategies \u00d7 1,000 simulations)"),
   sprintf("PSA rows: %d (expected 43,000); sims: %d (expected 1,000)", nrow(psa), n_psa_sims))

psa_ref <- psa[psa$strategy == "Passive case finding", c("sim","cost","qaly")]
names(psa_ref)[2:3] <- c("ref_cost","ref_qaly")

psa_cxr  <- psa[psa$strategy == "Cough+CXR (TB sx)", ] %>%
  left_join(psa_ref, by="sim") %>%
  mutate(icer = ifelse((qaly - ref_qaly) > 0, (cost - ref_cost) / (qaly - ref_qaly), NA_real_))
pce_cxr_25  <- mean(psa_cxr$icer < 25000, na.rm=TRUE) * 100
pce_cxr_35  <- mean(psa_cxr$icer < 35000, na.rm=TRUE) * 100
okw(pce_cxr_25 >= 99,
    sprintf("Cough+CXR P(CE at \u00a325k/\u00a335k) = %.1f%%/%.1f%% (expected ~100%%)", pce_cxr_25, pce_cxr_35),
    sprintf("Cough+CXR P(CE at \u00a325k) = %.1f%% (expected ~100%%)", pce_cxr_25))

psa_qft  <- psa[psa$strategy == "Parallel Sx+QFT (Ultra)", ] %>%
  left_join(psa_ref, by="sim") %>%
  mutate(icer = ifelse((qaly - ref_qaly) > 0, (cost - ref_cost) / (qaly - ref_qaly), NA_real_))
pce_qft_25  <- mean(psa_qft$icer < 25000, na.rm=TRUE) * 100
pce_qft_35  <- mean(psa_qft$icer < 35000, na.rm=TRUE) * 100
okw(pce_qft_25 >= 30 && pce_qft_25 <= 100,
    sprintf("Parallel Sx+QFT P(CE) pairwise NHS: %.1f%% (\u00a325k) / %.1f%% (\u00a335k)", pce_qft_25, pce_qft_35),
    sprintf("Parallel Sx+QFT P(CE at \u00a325k) = %.1f%% (unexpectedly low)", pce_qft_25))

# PSA — societal perspective (transmission savings applied)
tx_psa_sum <- read.csv("output/csv/transmission_psa_summary.csv", stringsAsFactors = FALSE)
if ("prob_ce_25k" %in% names(tx_psa_sum)) {
  pce_tx <- tx_psa_sum[tx_psa_sum$strategy == "Parallel Sx+QFT (Ultra)", "prob_ce_25k"]
  okw(length(pce_tx)==1 && pce_tx > 0.70,
      sprintf("Parallel Sx+QFT P(CE) societal @ \u00a325k = %.1f%%", pce_tx * 100),
      sprintf("Parallel Sx+QFT P(CE) societal = %.1f%% (expected >70%%)", pce_tx * 100))
} else {
  WARN(paste0("transmission_psa_summary.csv cols: ",
              paste(names(tx_psa_sum), collapse=", "), " — update column name in verify_project.R"))
}

# =============================================================================
HEAD("6. OUTPUT FILES PRESENT")
# =============================================================================

required <- c(
  # CSVs
  "output/csv/icer_nhs_perspective.csv",
  "output/csv/icer_societal_perspective.csv",
  "output/csv/ltbi_completion_sensitivity.csv",
  "output/csv/ltbi_prevalence_sensitivity.csv",
  "output/csv/ltbi_efficacy_sensitivity.csv",
  "output/csv/active_tb_prev_sensitivity.csv",
  "output/csv/igra_programme_cost_sensitivity.csv",
  "output/csv/igra_uptake_sensitivity.csv",
  "output/csv/psa_results.csv",
  "output/csv/icer_confidence_intervals.csv",
  "output/csv/dsa_results.csv",
  "output/csv/transmission_prevention_counts.csv",
  "output/csv/post_tb_sequelae_sensitivity.csv",
  "output/csv/sequential_ltbi_sa.csv",
  # Body figures (TIFF)
  "output/body/ce_scatter.tiff",
  "output/body/ceac_panel.tiff",
  "output/body/tornado_panel.tiff",
  "output/body/transmission_prevention_counts.tiff",
  # Supplementary figures — patient flow (manual drawio exports)
  "output/supplementary/model_scheme_supplementary.png",
  "output/supplementary/analytical_framework_supplementary.png",
  "output/supplementary/patientflow_CXRonly_supplementary.png",
  "output/supplementary/patientflow_CXRandQFT_supplementary.png",
  "output/supplementary/patientflow_CXRUltra_supplementary.png",
  "output/supplementary/patientflow_CoughCXR_supplementary.png",
  # Supplementary figures — model output (generated by R scripts)
  "output/supplementary/ltbi_detection_gap_supplementary.tiff",
  "output/supplementary/igra_specificity_supplementary.tiff",
  "output/supplementary/igra_uptake_supplementary.tiff",
  "output/supplementary/ltbi_completion_supplementary.tiff",
  "output/supplementary/ltbi_prevalence_supplementary.tiff",
  "output/supplementary/ltbi_sa_combined_supplementary.tiff",
  "output/supplementary/active_tb_prevalence_supplementary.tiff",
  "output/supplementary/psa_ce_plane_supplementary.tiff",
  # Paper supplementary tables
  "output/supplementary/CHEERS_supplementary.xlsx",
  "output/supplementary/Model_Parameters_supplementary.xlsx",
  "output/supplementary/CE_Results_supplementary.xlsx",
  "output/supplementary/Sensitivity_Analyses_supplementary.xlsx",
  "output/supplementary/PSA_Results_supplementary.xlsx",
  "input/Screening_Decision_Tree.xlsx"
)

# ResultsTable — date-stamped; check any exists
results_tables <- Sys.glob("output/body/ResultsTable_paper_*.xlsx")
ok(length(results_tables) > 0,
   sprintf("ResultsTable_paper found: %s", basename(results_tables[length(results_tables)])),
   "No output/body/ResultsTable_paper_*.xlsx found — run Rscript scripts/run_plots_and_tables.R")

missing <- required[!file.exists(required)]
ok(length(missing) == 0,
   sprintf("All %d required output files present", length(required)),
   paste(c(sprintf("%d missing files:", length(missing)), paste(" ", missing)), collapse="\n"))

# =============================================================================
HEAD("7. KNOWN RESULTS SPOT-CHECK")
# =============================================================================

spot_checks <- list(
  list("Passive case finding",      "cost_per_person", 203,    5),
  list("Cough+CXR (TB sx)",         "cost_per_person", 209,    5),
  list("Parallel Sx+QFT (Ultra)",   "cost_per_person", 363,    5),
  list("Cough+CXR (TB sx)",         "sequential_icer", 2205,   100),
  list("Symptom screen+CXR",        "sequential_icer", 7963,   200),
  list("Parallel Sx+QFT (Ultra)",   "sequential_icer", 14344,  200),
  list("Parallel Cough+QFT (Ultra)","ltbi_detected",    8644,  20)
)

for (chk in spot_checks) {
  strat <- chk[[1]]; col <- chk[[2]]; exp <- chk[[3]]; tol <- chk[[4]]
  row   <- bc[bc$strategy == strat, ]
  if (nrow(row) == 0) { FAIL(sprintf("Strategy '%s' not found", strat)); next }
  actual <- row[[col]]
  ok(!is.na(actual) && abs(actual - exp) <= tol,
     sprintf("%-40s %-20s = \u00a3%s", strat, col, format(round(actual), big.mark=",")),
     sprintf("%-40s %-20s = \u00a3%s (expected ~\u00a3%s, tol %s)",
             strat, col, format(round(actual), big.mark=","), format(exp, big.mark=","), tol))
}

# =============================================================================
HEAD("8. DOMINANCE COUNTS")
# =============================================================================

n_ref  <- sum(bc$dominance == "ref")
n_nd   <- sum(bc$dominance == "non-dominated")
n_dom  <- sum(bc$dominance == "simply dominated")
n_ext  <- sum(bc$dominance == "extendedly dominated")
cat(sprintf("  Breakdown: %d ref | %d non-dominated | %d simply dominated | %d extendedly dominated\n",
            n_ref, n_nd, n_dom, n_ext))
ok(n_ref==1 && n_nd==3 && n_dom+n_ext==39,
   sprintf("Dominance counts correct: 1 ref | 3 non-dominated | %d simply | %d extendedly dominated", n_dom, n_ext),
   sprintf("Unexpected counts: ref=%d non-dom=%d simply=%d ext=%d (expected 1+3+39=43)",
           n_ref, n_nd, n_dom, n_ext))

# =============================================================================
HEAD("9. ADDITIONAL ROBUSTNESS CHECKS")
# =============================================================================

# 9a. DSA parameter bounds validity (catch alpha↔beta swaps and mis-specified dists)
dsa_df <- read.csv("output/csv/dsa_results.csv", stringsAsFactors = FALSE)
bad_bounds <- dsa_df[!(dsa_df$lo_value < dsa_df$base_value &
                        dsa_df$base_value < dsa_df$hi_value), ]
ok(nrow(bad_bounds) == 0,
   sprintf("DSA bounds valid: lo < base < hi for all %d parameters", nrow(dsa_df)),
   sprintf("%d DSA parameters have lo >= base or base >= hi (possible alpha/beta swap): %s",
           nrow(bad_bounds), paste(head(bad_bounds$parameter, 5), collapse=", ")))

# 9b. Frontier sequential ICERs are strictly monotonically increasing
f_seq <- frontier[!is.na(frontier$sequential_icer), "sequential_icer"]
ok(length(f_seq) >= 2 && all(diff(f_seq) > 0),
   sprintf("Sequential frontier ICERs strictly increasing: %s",
           paste(paste0("\u00a3", format(round(f_seq), big.mark=",")), collapse=" < ")),
   sprintf("Sequential frontier ICERs NOT monotonically increasing: %s",
           paste(paste0("\u00a3", format(round(f_seq), big.mark=",")), collapse=", ")))

# 9c. CEAF: P(optimal) sums to 1.0 at each WTP point (from PSA results)
# At each WTP, exactly one strategy maximises NMB per sim → sum of P(optimal) = 1
psa_ref_cost <- psa[psa$strategy == "Passive case finding", c("sim","cost","qaly")]
names(psa_ref_cost)[2:3] <- c("ref_cost","ref_qaly")
for (wtp in c(25000, 35000)) {
  psa_nmb <- psa %>%
    left_join(psa_ref_cost, by = "sim") %>%
    mutate(inc_nmb = (qaly - ref_qaly) * wtp - (cost - ref_cost))
  psa_opt <- psa_nmb %>%
    group_by(sim) %>%
    slice_max(inc_nmb, n = 1, with_ties = FALSE) %>%
    ungroup()
  prob_opt <- psa_opt %>%
    group_by(strategy) %>%
    summarise(p = n() / n_psa_sims, .groups = "drop")
  p_sum <- sum(prob_opt$p)
  ok(abs(p_sum - 1.0) < 0.001,
     sprintf("CEAF probabilities sum to 1.00 at \u00a3%s/QALY (sum = %.4f)",
             format(wtp, big.mark=","), p_sum),
     sprintf("CEAF probabilities sum to %.4f at \u00a3%s/QALY (expected 1.000)",
             p_sum, format(wtp, big.mark=",")))
}

# 9d. Active TB prevalence SA: Parallel Sx+QFT remains on frontier at 0.44% prevalence
atb_sa <- read.csv("output/csv/active_tb_prev_sensitivity.csv", stringsAsFactors = FALSE)
atb_highrisk <- atb_sa[atb_sa$scenario == "highrisk_044pct", ]
qft_front_044 <- atb_highrisk[atb_highrisk$strategy == "Parallel Sx+QFT (Ultra)" &
                               atb_highrisk$dominance %in% c("ref","non-dominated"), ]
okw(nrow(qft_front_044) > 0,
    "SA (prev=0.44%): Parallel Sx+QFT remains on efficient frontier",
    "SA (prev=0.44%): Parallel Sx+QFT NOT on frontier at 0.44% prevalence — check SA output")

# =============================================================================
HEAD("10. NARRATIVE CLAIMS RECONCILIATION")
# =============================================================================
# Verify every numeric claim in paper/narrative.md against live CSV outputs.
# Any FAIL here means the narrative is outdated and must be updated before submission.

# 10a. Frontier costs and ICERs
narr_icer <- list(
  list("nhs",      "Passive case finding",    "cost_per_person",  203,    5,   "Para 1 / Table 2"),
  list("nhs",      "Cough+CXR (TB sx)",       "cost_per_person",  209,    5,   "Para 1 / Table 2"),
  list("nhs",      "Symptom screen+CXR",      "cost_per_person",  220,    5,   "Para 1 / Table 2"),
  list("nhs",      "Parallel Sx+QFT (Ultra)", "cost_per_person",  361,    5,   "Para 1 (NHS)"),
  list("societal", "Parallel Sx+QFT (Ultra)", "cost_per_person",  345,    5,   "Para 1/2 (societal)"),
  list("nhs",      "Cough+CXR (TB sx)",       "sequential_icer",  2205,   100, "Para 1 / abstract"),
  list("nhs",      "Symptom screen+CXR",      "sequential_icer",  7963,   200, "Para 1"),
  list("nhs",      "Parallel Sx+QFT (Ultra)", "sequential_icer",  14344,  300, "Para 1 / abstract (NHS)"),
  list("societal", "Parallel Sx+QFT (Ultra)", "sequential_icer",  12738,  300, "Para 1 / abstract (societal)"),
  list("nhs",      "Passive case finding",    "qaly_per_person",  21.606, 0.01,"Para 1")
)

for (chk in narr_icer) {
  persp <- chk[[1]]; strat <- chk[[2]]; col <- chk[[3]]
  exp <- chk[[4]]; tol <- chk[[5]]; src <- chk[[6]]
  tbl <- if (persp == "societal") tx else bc
  row <- tbl[tbl$strategy == strat, ]
  if (nrow(row) == 0) { FAIL(sprintf("Narrative claim: '%s' not found in %s table", strat, persp)); next }
  actual <- row[[col]]
  ok(!is.na(actual) && abs(actual - exp) <= tol,
     sprintf("Narrative (%s): %-38s %-20s = %s",
             persp, strat, col, format(round(actual, 2), big.mark = ",")),
     sprintf("NARRATIVE STALE [%s]: %s %s %s = %s (narrative claims %s)",
             src, persp, strat, col,
             format(round(actual, 2), big.mark = ","),
             format(round(exp,    2), big.mark = ",")))
}

# 10b. Transmission counts (Para 2)
tx_cnt <- read.csv("output/csv/transmission_prevention_counts.csv", stringsAsFactors = FALSE)
psq_tx <- tx_cnt[tx_cnt$strategy == "Parallel Sx+QFT (Ultra)", ]

narr_tx <- list(
  list("atb_cases_prevented",    889,  50,  "Para 2: reactivations prevented"),
  list("secondary_cases_avoided", 182, 20,  "Para 2: secondary active TB avoided"),
  list("secondary_ltbi_avoided", 1023, 50,  "Para 2: secondary LTBI avoided"),
  list("transmission_saving_pp",  15.73, 1, "Para 2: £15.73/pp savings")
)

for (chk in narr_tx) {
  col <- chk[[1]]; exp <- chk[[2]]; tol <- chk[[3]]; src <- chk[[4]]
  if (nrow(psq_tx) == 0 || !col %in% names(psq_tx)) { WARN(sprintf("Transmission column '%s' missing", col)); next }
  actual <- psq_tx[[col]]
  ok(!is.na(actual) && abs(actual - exp) <= tol,
     sprintf("Narrative (%s): %s = %s", src, col, format(round(actual, 2), big.mark = ",")),
     sprintf("NARRATIVE STALE [%s]: %s = %s (narrative claims %s)",
             src, col, format(round(actual, 2), big.mark = ","), exp))
}

# 10c. PSA P(optimal) exact values vs narrative
psa_ref_cost2 <- psa[psa$strategy == "Passive case finding", c("sim","cost","qaly")]
names(psa_ref_cost2)[2:3] <- c("ref_cost", "ref_qaly")

compute_p_optimal <- function(wtp) {
  # Filter to frontier strategies only — consistent with CEAC in paper (standard HTA practice).
  # All-strategy P(optimal) gives ~41% for Parallel Sx+QFT because some dominated strategies
  # can have highest NMB in individual PSA draws; this is not what is reported in the paper.
  psa %>%
    filter(strategy %in% expected_frontier) %>%
    left_join(psa_ref_cost2, by = "sim") %>%
    mutate(inc_nmb = (qaly - ref_qaly) * wtp - (cost - ref_cost)) %>%
    group_by(sim) %>%
    slice_max(inc_nmb, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    count(strategy) %>%
    mutate(p = n / n_psa_sims)
}

popt_25 <- compute_p_optimal(25000)
popt_35 <- compute_p_optimal(35000)

get_p <- function(df, strat) { v <- df$p[df$strategy == strat]; if (length(v)) v else 0 }

narr_psa <- list(
  list(popt_25, "Parallel Sx+QFT (Ultra)", 0.848, 0.025, "Para 4 (NHS £25k)"),
  list(popt_35, "Parallel Sx+QFT (Ultra)", 0.963, 0.025, "Para 4 (NHS £35k)"),
  list(popt_25, "Cough+CXR (TB sx)",       0.152, 0.025, "Para 4 (Cough+CXR £25k)"),
  list(popt_25, "Symptom screen+CXR",      0.000, 0.010, "Para 4 (Symptom screen 0%)")
)

for (chk in narr_psa) {
  df <- chk[[1]]; strat <- chk[[2]]; exp <- chk[[3]]; tol <- chk[[4]]; src <- chk[[5]]
  actual <- get_p(df, strat)
  ci_lo  <- actual - 1.96 * sqrt(actual * (1 - actual) / n_psa_sims)
  ci_hi  <- actual + 1.96 * sqrt(actual * (1 - actual) / n_psa_sims)
  ok(abs(actual - exp) <= tol,
     sprintf("Narrative PSA [%s]: P(optimal) = %.1f%% [95%% CI %.1f%%\u2013%.1f%%]",
             src, actual * 100, ci_lo * 100, ci_hi * 100),
     sprintf("NARRATIVE STALE [%s]: P(optimal) = %.1f%% (narrative claims %.1f%%)",
             src, actual * 100, exp * 100))
}

# 10d. LTBI treatment efficacy SA — Parallel Sx+QFT dominated at NICE NG33 range
eff_csv_path <- "output/csv/ltbi_efficacy_sensitivity.csv"
if (file.exists(eff_csv_path)) {
  ltbi_eff <- read.csv(eff_csv_path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ltbi_eff))) {
    is_dom <- ltbi_eff$parallel_dominance[i] %in%
                c("simply dominated", "extendedly dominated")
    ok(is_dom,
       sprintf("LTBI efficacy SA [%s]: Parallel Sx+QFT = %s (expected: dominated)",
               ltbi_eff$scenario[i], ltbi_eff$parallel_dominance[i]),
       sprintf("LTBI EFFICACY SA [%s]: Parallel Sx+QFT = %s but narrative states dominated",
               ltbi_eff$scenario[i], ltbi_eff$parallel_dominance[i]))
  }
} else {
  WARN("output/csv/ltbi_efficacy_sensitivity.csv not found — run run_sa.R first")
}

# =============================================================================
HEAD("11. TRANSITION PROBABILITY CONSISTENCY")
# =============================================================================
# Verify outflow probabilities from key states sum to ≤ 1.0.
# remainder = diagonal (stay in state). Sum > 1 = impossible transition matrix.

tp_checks <- list(
  list("LatentTreated",
       c("p_LatentTreated_LatentCompleted", "p_LatentTreated_LatentDiscontinued",
         "p_LatentTreated_LatentLtfu"),
       "LTBI treatment outflows"),
  list("ActiveTreated",
       c("p_ActiveTreated_ActiveCompleted", "p_ActiveTreated_ActiveDiscontinued",
         "p_ActiveTreated_ActiveLtfu"),
       "Active TB treatment outflows"),
  list("LatentDiagnosed",
       c("p_LatentDiagnosed_LatentTreated", "p_LatentDiagnosed_LatentNotreated",
         "p_LatentDiagnosed_ActiveUndiagnosed"),
       "LTBI diagnosed outflows"),
  list("ActiveUndiagnosed",
       c("p_ActiveUndiagnosed_ActiveDiagnosed", "p_ActiveUndiagnosed_Dead"),
       "Active undiagnosed outflows (excl. treatment-seeking)")
)

for (chk in tp_checks) {
  state <- chk[[1]]; params <- chk[[2]]; desc <- chk[[3]]
  vals <- sapply(params, cfg_val)
  if (any(is.na(vals))) {
    WARN(sprintf("TP check: missing params for %s state (%s)", state, paste(params[is.na(vals)], collapse=",")))
    next
  }
  total_out <- sum(vals)
  ok(total_out > 0 && total_out <= 1.0,
     sprintf("TP outflows %-20s sum = %.5f (diagonal = %.5f) \u2014 %s", state, total_out, 1-total_out, desc),
     sprintf("IMPOSSIBLE TP: %-20s outflows sum = %.5f > 1.0 \u2014 %s", state, total_out, desc))
}

# Phase 2 reactivation rate < phase 1 (switch at cycle 61 = 5 years)
react1 <- cfg_val("p_LatentUndiagnosed_ActiveUndiagnosed")   # phase 1: 0.000912/mo
react2 <- cfg_val("p_react_phase2")                           # phase 2: 0.000333/mo
ok(!is.na(react1) && !is.na(react2) && react2 < react1,
   sprintf("Phase 2 reactivation (%.6f/mo) < phase 1 (%.6f/mo) \u2014 two-phase model correct",
           react2, react1),
   sprintf("Phase 2 reactivation (%.6f) not < phase 1 (%.6f) \u2014 check model code",
           react2, react1))

ok(!is.na(react1) && abs(react1 - 0.000912) < 1e-5,
   sprintf("Phase 1 reactivation: %.6f/mo (Berrocal-Almanza 2022 1.09%%/yr confirmed)", react1),
   sprintf("Phase 1 reactivation: %.6f/mo (expected 0.000912)", react1))

ok(!is.na(react2) && abs(react2 - 0.000333) < 1e-5,
   sprintf("Phase 2 reactivation: %.6f/mo (Horsburgh 2004 ~0.4%%/yr confirmed)", react2),
   sprintf("Phase 2 reactivation: %.6f/mo (expected 0.000333)", react2))

# =============================================================================
HEAD("12. CONFIG.CSV COMPLETENESS")
# =============================================================================

# No duplicate parameter names
dup_params <- cfg$parameter[duplicated(cfg$parameter)]
ok(length(dup_params) == 0,
   sprintf("No duplicate parameter names (%d parameters in config.csv)", nrow(cfg)),
   sprintf("Duplicate parameters: %s", paste(dup_params, collapse = ", ")))

# No NA or empty values in rows that look like real parameters (valid R identifiers)
real_param_rows <- cfg[grepl("^[a-zA-Z_][a-zA-Z0-9_.]*$", cfg$parameter), ]
na_vals <- real_param_rows[is.na(real_param_rows$value) | trimws(as.character(real_param_rows$value)) == "", "parameter"]
ok(length(na_vals) == 0,
   sprintf("No NA/empty values in config.csv (%d parameters checked)", nrow(real_param_rows)),
   sprintf("%d parameters with NA/empty value: %s", length(na_vals), paste(na_vals, collapse = ", ")))

# Initial cohort distribution: use init_* parameters (base-case proportions)
init_ltbi   <- cfg_val("init_LatentUndiagnosed")   # 0.178 = 17.8%
init_active <- cfg_val("init_ActiveUndiagnosed")    # 0.010 = 1.0%
init_uninf  <- cfg_val("init_Uninfected")           # 0.812 = 81.2%

if (!is.na(init_ltbi) && !is.na(init_active) && !is.na(init_uninf)) {
  init_sum <- init_ltbi + init_active + init_uninf
  ok(abs(init_sum - 1.0) < 1e-6,
     sprintf("Initial cohort sums to 1.0: uninfected %.1f%% + LTBI %.1f%% + active %.1f%%",
             init_uninf*100, init_ltbi*100, init_active*100),
     sprintf("Initial cohort sums to %.6f (expected 1.000000)", init_sum))
  ok(abs(init_ltbi - 0.178) < 0.005,
     sprintf("LTBI prevalence at entry: %.1f%% (narrative: 17.8%%)", init_ltbi*100),
     sprintf("LTBI prevalence: %.1f%% (narrative claims 17.8%%)", init_ltbi*100))
  ok(abs(init_active - 0.01) < 0.002,
     sprintf("Active TB prevalence at entry: %.1f%% (base case, narrative: 1.0%%)", init_active*100),
     sprintf("Active TB prevalence: %.1f%% (expected 1.0%%)", init_active*100))
  ok(abs(init_uninf - 0.812) < 0.005,
     sprintf("Uninfected at entry: %.1f%% (= 1 \u2212 17.8%% \u2212 1.0%%)", init_uninf*100),
     sprintf("Uninfected at entry: %.1f%% (expected 81.2%%)", init_uninf*100))
}

# Discount rate consistent with n_t and NICE guidance
ok(!is.na(cfg_val("n_t")) && cfg_val("n_t") == 660,
   "n_t = 660 cycles (55-year horizon \u00d7 12 months, NICE reference case)",
   sprintf("n_t = %s (expected 660)", cfg_val("n_t")))

# =============================================================================
HEAD("13. ICER ARITHMETIC VERIFICATION")
# =============================================================================
# Manually recompute sequential ICERs from raw cost/QALY columns.
# Catches discounting or cycle-length bugs independently of model internals.

frontier_arith <- bc[bc$dominance %in% c("ref","non-dominated"), ]
frontier_arith  <- frontier_arith[order(frontier_arith$cost_per_person), ]

for (i in seq(2, nrow(frontier_arith))) {
  d_cost <- frontier_arith$cost_per_person[i] - frontier_arith$cost_per_person[i - 1]
  d_qaly <- frontier_arith$qaly_per_person[i] - frontier_arith$qaly_per_person[i - 1]
  manual  <- if (!is.na(d_qaly) && d_qaly > 0) d_cost / d_qaly else NA_real_
  reported <- frontier_arith$sequential_icer[i]
  if (is.na(manual) || is.na(reported)) next
  ok(abs(manual - reported) < 200,
     sprintf("ICER arithmetic: %-32s manual \u00a3%s \u2248 reported \u00a3%s",
             frontier_arith$strategy[i],
             format(round(manual), big.mark = ","),
             format(round(reported), big.mark = ",")),
     sprintf("ICER mismatch: %-32s manual=\u00a3%s vs reported=\u00a3%s",
             frontier_arith$strategy[i],
             format(round(manual),   big.mark = ","),
             format(round(reported), big.mark = ",")))
}

# =============================================================================
HEAD("14. DSA NMB MONOTONICITY")
# =============================================================================
# For parameters with a known expected direction, verify nmb_hi vs nmb_lo.
# A reversed relationship signals a sign error in the sensitivity range.
# nmb_lo = NMB when parameter at lo value; nmb_hi = NMB when at hi value.
# "positive": higher parameter value → higher NMB (more CE)
# "negative": higher parameter value → lower NMB

dsa_mono <- list(
  list("p_LatentUndiagnosed_ActiveUndiagnosed", "positive",
       "Reactivation rate: more disease \u2192 more screening benefit"),
  list("p_LatentDiagnosed_LatentTreated",       "positive",
       "LTBI treatment initiation: more treated \u2192 more benefit"),
  list("p_LatentTreated_LatentCompleted",        "positive",
       "LTBI treatment completion: more cures \u2192 more benefit"),
  list("qaly_LatentUndiagnosed",                "negative",
       "Latent QALY: higher utility \u2192 less disease burden \u2192 less CE"),
  list("p_ActiveUndiagnosed_Dead",              "positive",
       "Active TB mortality: higher \u2192 more benefit from detection")
)

for (chk in dsa_mono) {
  param <- chk[[1]]; dir <- chk[[2]]; desc <- chk[[3]]
  row   <- dsa_df[dsa_df$parameter == param, ]
  if (nrow(row) == 0) { WARN(sprintf("DSA monotonicity: '%s' not found in dsa_results.csv", param)); next }
  row <- row[1, ]
  if (dir == "positive") {
    ok(row$nmb_hi > row$nmb_lo,
       sprintf("DSA monotonic (+): %s", desc),
       sprintf("DSA NOT monotonic (+): %s \u2014 nmb_hi (%.0f) \u2264 nmb_lo (%.0f); check parameter direction",
               param, row$nmb_hi, row$nmb_lo))
  } else {
    ok(row$nmb_hi < row$nmb_lo,
       sprintf("DSA monotonic (\u2212): %s", desc),
       sprintf("DSA NOT monotonic (\u2212): %s \u2014 nmb_hi (%.0f) \u2265 nmb_lo (%.0f); check parameter direction",
               param, row$nmb_hi, row$nmb_lo))
  }
}

# =============================================================================
HEAD("15. PSA STATISTICS & PRECISION")
# =============================================================================
# Compute 95% binomial CI for P(optimal) and flag if narrative claims fall outside.
# With 1,000 sims: SE ≈ 1% for p=0.85; 95% CI ≈ ±2.2pp.

pce_ci_checks <- list(
  list(wtp = 25000, strat = "Parallel Sx+QFT (Ultra)", claimed = 0.848, label = "NHS £25k"),
  list(wtp = 35000, strat = "Parallel Sx+QFT (Ultra)", claimed = 0.963, label = "NHS £35k"),
  list(wtp = 25000, strat = "Cough+CXR (TB sx)",       claimed = 0.152, label = "NHS £25k"),
  list(wtp = 25000, strat = "Symptom screen+CXR",      claimed = 0.000, label = "NHS £25k")
)

popt_dfs <- list("25000" = popt_25, "35000" = popt_35)

for (chk in pce_ci_checks) {
  wtp_key <- as.character(chk$wtp)
  actual  <- get_p(popt_dfs[[wtp_key]], chk$strat)
  se      <- sqrt(actual * (1 - actual) / n_psa_sims)
  ci_lo   <- max(0, actual - 1.96 * se)
  ci_hi   <- min(1, actual + 1.96 * se)
  in_ci   <- chk$claimed >= ci_lo & chk$claimed <= ci_hi
  okw(in_ci,
      sprintf("PSA CI: %-30s @ %s: P=%.1f%% [%.1f%%-%.1f%%]; narrative %.1f%% within CI",
              chk$strat, chk$label, actual*100, ci_lo*100, ci_hi*100, chk$claimed*100),
      sprintf("PSA CI: %-30s @ %s: P=%.1f%% [%.1f%%-%.1f%%]; narrative %.1f%% OUTSIDE CI",
              chk$strat, chk$label, actual*100, ci_lo*100, ci_hi*100, chk$claimed*100))
}

# PSA: Symptom screen+CXR must never be optimal (0% at any frontier WTP)
symscr_opt_25 <- get_p(popt_25, "Symptom screen+CXR")
symscr_opt_35 <- get_p(popt_35, "Symptom screen+CXR")
ok(symscr_opt_25 == 0 && symscr_opt_35 == 0,
   "PSA: Symptom screen+CXR P(optimal) = 0% at both £25k and £35k (never NMB-maximising)",
   sprintf("PSA: Symptom screen+CXR P(optimal) = %.1f%% (£25k) / %.1f%% (£35k) — unexpected",
           symscr_opt_25 * 100, symscr_opt_35 * 100))

# PSA: set.seed() check for reproducibility
master_lines <- tryCatch(readLines("MasterTBModel.R"), error = function(e) character(0))
has_seed <- any(grepl("set\\.seed", master_lines))
okw(has_seed,
    "PSA: set.seed() found in MasterTBModel.R \u2014 results reproducible",
    "PSA: no set.seed() in MasterTBModel.R \u2014 PSA results vary between runs; document in narrative")

# =============================================================================
HEAD("16. FIGURE-TABLE CONSISTENCY")
# =============================================================================
suppressPackageStartupMessages({ library(readxl) })

# icer_table.xlsx: verify row count = 43 strategies
icer_xlsx_file <- "output/body/icer_table.xlsx"
if (file.exists(icer_xlsx_file)) {
  icer_xl <- tryCatch(suppressMessages(read_excel(icer_xlsx_file, col_names = FALSE)),
                      error = function(e) NULL)
  if (!is.null(icer_xl)) {
    # Rows that are strategy data (not headers, not footnotes, not NA)
    col1 <- as.character(icer_xl[[1]])
    strat_rows <- col1[!is.na(col1) &
                         !grepl("^(Strategy|ICER table|NHS direct|Highlighted)", col1, ignore.case = TRUE) &
                         trimws(col1) != ""]
    ok(length(strat_rows) == 43,
       sprintf("icer_table.xlsx: 43 strategy rows confirmed"),
       sprintf("icer_table.xlsx: %d strategy rows (expected 43)", length(strat_rows)))

    # Verify Parallel Sx+QFT ICER in xlsx matches CSV
    psq_xlsx_row <- icer_xl[col1 == "Parallel Sx+QFT (Ultra)" & !is.na(col1), ]
    if (nrow(psq_xlsx_row) > 0) {
      icer_xlsx_val <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(psq_xlsx_row[[7]]))))
      icer_csv_val  <- bc[bc$strategy == "Parallel Sx+QFT (Ultra)", "sequential_icer"]
      if (!is.na(icer_xlsx_val) && length(icer_csv_val) == 1) {
        ok(abs(icer_xlsx_val - icer_csv_val) < 200,
           sprintf("icer_table.xlsx: Parallel Sx+QFT ICER \u00a3%s matches CSV \u00a3%s",
                   format(round(icer_xlsx_val), big.mark=","),
                   format(round(icer_csv_val),  big.mark=",")),
           sprintf("icer_table.xlsx ICER mismatch: xlsx=\u00a3%s csv=\u00a3%s",
                   format(round(icer_xlsx_val), big.mark=","),
                   format(round(icer_csv_val),  big.mark=",")))
      }
    }
  }
} else {
  WARN("icer_table.xlsx not found \u2014 run scripts/run_plots_and_tables.R")
}

# ResultsTable: Body Table 2 frontier rows match icer_societal_perspective.csv
rt_file <- tail(Sys.glob("output/body/ResultsTable_paper_*.xlsx"), 1)
if (length(rt_file) == 1) {
  rt_xl <- tryCatch(suppressMessages(read_excel(rt_file, sheet = 1, col_names = FALSE)),
                    error = function(e) NULL)
  if (!is.null(rt_xl)) {
    rt_col1 <- as.character(rt_xl[[1]])
    rt_strats <- rt_col1[!is.na(rt_col1) &
                           !grepl("^(Strategy|Table|55-year|\\+IGRA)", rt_col1, ignore.case=TRUE) &
                           trimws(rt_col1) != ""]
    ok(length(rt_strats) == 4,
       sprintf("ResultsTable Body Table 2: 4 frontier strategy rows"),
       sprintf("ResultsTable Body Table 2: %d rows (expected 4)", length(rt_strats)))
  }
}

# =============================================================================
HEAD("17. LIMITATIONS DISCLOSURE CHECK")
# =============================================================================
narr_lines <- tryCatch(readLines("paper/narrative.md"), error = function(e) character(0))
narr_str   <- paste(narr_lines, collapse = " ")

lim_checks <- list(
  list("SAGER|sex.disaggregat",          "SAGER / sex disaggregation limitation"),
  list("emigrat|return.migrat",           "Emigration / return migration limitation"),
  list("100.*efficac|residual.*reactivat|NICE NG33.*10", "LTBI 100% efficacy limitation (residual reactivation)"),
  list("post.TB sequelae|bronchiectasis", "Post-TB sequelae limitation"),
  list("passive case finding calibrat",   "Passive case finding calibration assumption")
)

for (chk in lim_checks) {
  pattern <- chk[[1]]; desc <- chk[[2]]
  okw(grepl(pattern, narr_str, ignore.case = TRUE, perl = TRUE),
      sprintf("Limitation in narrative: %s", desc),
      sprintf("Limitation may be missing from narrative Discussion: %s", desc))
}

# =============================================================================
HEAD("18. PLOT QUALITY AUDIT")
# =============================================================================

# 18a — File size floor checks (catches blank/failed renders)
# Floors are set for LZW-compressed TIFFs; a blank/failed render is < 5 KB
plot_size_floors <- list(
  "output/body/ce_scatter.tiff"                                    = 100e3,
  "output/body/ceac_panel.tiff"                                    = 100e3,
  "output/body/tornado_panel.tiff"                                 = 100e3,
  "output/body/transmission_prevention_counts.tiff"                = 50e3,
  "output/supplementary/igra_specificity_supplementary.tiff"       = 100e3,
  "output/supplementary/igra_uptake_supplementary.tiff"            = 50e3,
  "output/supplementary/psa_ce_plane_supplementary.tiff"           = 100e3,
  "output/supplementary/active_tb_prevalence_supplementary.tiff"   = 100e3,
  "output/supplementary/ltbi_sa_combined_supplementary.tiff"       = 200e3,
  "output/supplementary/ltbi_detection_gap_supplementary.tiff"     = 100e3
)
for (f in names(plot_size_floors)) {
  sz  <- file.info(f)$size
  thr <- plot_size_floors[[f]]
  ok(!is.na(sz) && sz >= thr,
     sprintf("Plot size \u2265 %.0f KB: %s (%.0f KB)", thr / 1024,
             basename(f), if (is.na(sz)) 0 else sz / 1024),
     sprintf("Plot too small / blank: %s (%.0f KB < %.0f KB floor)",
             basename(f), if (is.na(sz)) 0 else sz / 1024, thr / 1024))
}

# 18b — CSV–plot consistency: active TB prevalence SA must have >= 3 scenarios
if (file.exists("output/csv/active_tb_prev_sensitivity.csv")) {
  atb_csv <- read.csv("output/csv/active_tb_prev_sensitivity.csv",
                      stringsAsFactors = FALSE)
  n_prev <- length(unique(atb_csv$active_tb_prev))
  ok(n_prev >= 3,
     sprintf("active_tb_prev_sensitivity.csv contains %d distinct prevalence values", n_prev),
     sprintf("active_tb_prev_sensitivity.csv has only %d prevalence value(s); expected >= 3", n_prev))
} else {
  FAIL("output/csv/active_tb_prev_sensitivity.csv missing — cannot check prevalence scenario count")
}

# =============================================================================
HEAD("SUMMARY")
# =============================================================================

total <- n_fails + n_warns
cat(sprintf("\n  \x1b[31mFailures:\x1b[0m  %d\n", n_fails))
cat(sprintf("  \x1b[33mWarnings:\x1b[0m  %d\n", n_warns))
if (n_fails == 0 && n_warns == 0) {
  cat("  \x1b[32mAll checks passed", "\u2014", "project verified.\x1b[0m\n\n")
  quit(status = 0)
} else if (n_fails == 0) {
  cat("  \x1b[33mNo failures", "\u2014", "review warnings above.\x1b[0m\n\n")
  quit(status = 0)
} else {
  cat("  \x1b[31mFix failures before submission.\x1b[0m\n\n")
  quit(status = 1)
}
