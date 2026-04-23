# run_sa.R — all sensitivity analyses: DSA, IGRA, LTBI, ATB prevalence, sequelae,
#             LTBI-prev, transmission analysis, transmission DSA, PSA (1,000 simulations)
# Can be run standalone: Rscript scripts/run_sa.R
# Or sourced from MasterTBModel.R (inherits strategy_results from global env).

# =============================================================================
# DETERMINISTIC SENSITIVITY ANALYSES
# =============================================================================

if (!exists("strategy_results")) {
  if (!file.exists("input/config.csv"))
    stop("Run from project root: Rscript scripts/run_sa.R")
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(parallel)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(ggsci)
    source("scripts/model_functions.R")
  })
  if (!exists("LANCET"))
    LANCET <- ggsci::pal_lancet("lanonc", alpha = 0.85)(8)
  strategies       <- create_diagnostic_strategies()
  strategy_results <- mclapply(strategies, run_strategy,
                               mc.cores = max(1L, parallel::detectCores() - 1L))
}
if (!exists("icer_table")) {
  icer_table <- calculate_icer(strategy_results)
}

# =============================================================================
# ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS (DSA)
#
# Each param varied to 2.5th/97.5th percentile (others fixed); NMB at £25k WTP.
# Tornado diagrams for 3 frontier strategies:
# (1) Cough+CXR (TB sx)       seq. ICER £2,205/QALY (NHS)
# (2) Symptom screen+CXR      seq. ICER £7,963/QALY (NHS)
# (3) Parallel Sx+QFT (Ultra) seq. ICER £11,324/QALY (NHS)
#
# Only Markov model parameters (costs, QALYs, transition probabilities) are varied
# here. Effective active TB detection rates (eff_active_*) and background prevalence
# are addressed in separate scenario analyses.
# =============================================================================
# =============================================================================
# IGRA SPECIFICITY SENSITIVITY ANALYSIS
#
# Base case uses test-specific IGRA specificities from Pai et al. 2008 meta-analysis
# (Ann Intern Med 149:177-184; BCG-vaccinated populations):
#   QFT-TB Gold Plus: 0.96 (95% CI 0.94-0.98)   [config: igra_spec_qft]
#   T-SPOT.TB:        0.93 (95% CI 0.86-1.00)    [config: igra_spec_tspt]
#
# True LTBI detection = Sensitivity × LTBI_pool − (1 − Specificity) × Uninfected
# At base case: QFT = 11,526/100k; T-SPOT = 9,980/100k
#
# This SA applies a uniform specificity sweep to both tests to show sensitivity
# across the plausible range:
#   0.90 — conservative lower bound
#   0.93 — T-SPOT.TB base case (Pai 2008 BCG-vaccinated)
#   0.96 — QFT-TB base case    (Pai 2008 BCG-vaccinated)
#   1.00 — naive upper bound (no false positives; shown for comparison only)
#
# For each spec value, 8 parallel IGRA strategies are re-initialised with adjusted
# LTBI detection; non-parallel strategy results held fixed from base case run.
# =============================================================================
cat("\n================================================================================\n")
cat("         IGRA SPECIFICITY SENSITIVITY — PARALLEL LTBI DETECTION                \n")
cat("================================================================================\n\n")

igra_spec_values  <- c(0.90, 0.93, 0.96, 1.00)
n_ltbi_global     <- n_c * prev_ltbi        # 17,800 per 100k
n_uninfected_100k <- n_c * prev_uninfected  # 81,200 per 100k

cat(sprintf("  %-8s  %-20s  %-20s  %-20s\n",
    "Spec.", "QFT LTBI det./100k", "T-SPOT LTBI det./100k", "Note"))
cat(strrep("-", 75), "\n")
for (spec in igra_spec_values) {
  ltbi_qft_s  <- (n_ltbi_global * igra_sens_qft  - (1 - spec) * n_uninfected_100k) * igra_uptake_base
  ltbi_tspt_s <- (n_ltbi_global * igra_sens_tspt - (1 - spec) * n_uninfected_100k) * igra_uptake_base
  note <- if (spec == 0.96) "QFT base case (Pai 2008 BCG-vaccinated)" else
          if (spec == 0.93) "T-SPOT base case (Pai 2008 BCG-vaccinated)" else
          if (spec == 1.00) "naive upper bound (no FP; shown for comparison)" else
          "conservative lower bound"
  cat(sprintf("  %-6.2f  %-20.0f  %-20.0f  %s\n",
      spec, ltbi_qft_s, ltbi_tspt_s, note))
}
cat("\n")

# Run frontier ICERs for each specificity level
# (Works directly from strategies list — no dependency on build_strategies internals)
igra_spec_results <- lapply(igra_spec_values, function(spec) {
  # Apply igra_uptake_base so SA is consistent with base-case (75% uptake)
  ltbi_qft_s  <- max(0, n_ltbi_global * igra_sens_qft  - (1 - spec) * n_uninfected_100k) * igra_uptake_base
  ltbi_tspt_s <- max(0, n_ltbi_global * igra_sens_tspt - (1 - spec) * n_uninfected_100k) * igra_uptake_base

  # Re-create only the 8 parallel IGRA strategies with adjusted LTBI detection
  # Adjust LatentDiagnosed + LatentUndiagnosed in their init vectors directly
  strats_spec <- strategies
  parallel_keys <- names(strats_spec)[sapply(names(strats_spec), function(nm)
    !is.null(strats_spec[[nm]]$prog_cost_type))]
  for (k in parallel_keys) {
    ltbi_s <- if (grepl("tspt", k, ignore.case = TRUE)) ltbi_tspt_s else ltbi_qft_s
    strats_spec[[k]]$init["LatentDiagnosed"]   <- ltbi_s
    strats_spec[[k]]$init["LatentUndiagnosed"] <- n_ltbi_global - ltbi_s
  }

  # Run only affected strategies; copy base results for non-parallel
  spec_results <- strategy_results
  for (k in parallel_keys) {
    r <- run_strategy(strats_spec[[k]])
    idx <- which(sapply(spec_results, function(x) x$strategy_name == r$strategy_name))
    if (length(idx) > 0) spec_results[[idx]] <- r else spec_results <- c(spec_results, list(r))
  }

  icer_s <- calculate_icer(spec_results)
  icer_s$igra_spec <- spec
  icer_s
})

igra_spec_df <- bind_rows(igra_spec_results)
write.csv(igra_spec_df %>% mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
          "output/csv/igra_specificity_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_specificity_sensitivity.csv\n\n")

# Print frontier comparison across specificity values
frontier_strats_spec <- c("Passive case finding", "Cough+CXR (TB sx)",
                           "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)",
                           "Parallel Sx+T-SPOT (Ultra)")
cat(sprintf("Efficient frontier ICERs by IGRA specificity (QFT base=0.96; T-SPOT base=0.93; Pai 2008 BCG):\n"))
cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n",
    "Strategy", "Spec=0.90", "Spec=0.93", "Spec=0.96", "Spec=1.00"))
cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n",
    "", "(conserv.)", "(T-SPOT)", "(QFT base)", "(naive UB)"))
cat(strrep("-", 80), "\n")
for (st in frontier_strats_spec) {
  icers <- sapply(igra_spec_values, function(spec) {
    r <- igra_spec_df[igra_spec_df$strategy == st & abs(igra_spec_df$igra_spec - spec) < 1e-9, ]
    if (nrow(r) == 0 || is.na(r$sequential_icer)) return("—")
    paste0("\u00a3", format(round(r$sequential_icer), big.mark = ","))
  })
  cat(sprintf("  %-32s  %10s  %10s  %10s  %10s\n", st, icers[1], icers[2], icers[3], icers[4]))
}
cat("\n")

# IGRA specificity sensitivity plot
# Caption note: specificity values sourced from Pai M et al. Meta-analysis of IGRAs.
# Ann Intern Med. 2008;149:177-184. QFT base 0.96; T-SPOT base 0.93 (BCG-vaccinated subgroup).
# Cough+CXR (TB sx) and Symptom screen+CXR have no IGRA component; their ICERs are
# invariant to IGRA specificity and are excluded from this plot. Only parallel
# IGRA strategies (where specificity affects LTBI detection) are shown.
spec_plot_strats <- c("Parallel Sx+QFT (Ultra)", "Parallel Cough+QFT (Ultra)")
spec_pal <- c(
  "Parallel Sx+QFT (Ultra)"    = LANCET[7],   # dark red
  "Parallel Cough+QFT (Ultra)" = "#7a8a99"    # steel gray
)
spec_df_plot <- igra_spec_df %>%
  filter(strategy %in% spec_plot_strats, !is.na(sequential_icer))

# y-axis range: parallel strategy ICERs sit well below £25k; pad to 20k
spec_y_max <- max(spec_df_plot$sequential_icer, na.rm = TRUE) * 1.15

p_igra_spec <- ggplot(spec_df_plot,
                      aes(x = igra_spec, y = sequential_icer, colour = strategy)) +
  geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 0.905, y = 25800,
           label = "\u00a325,000/QALY (NICE threshold)", hjust = 0, size = 3.0,
           colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = spec_pal, name = "Strategy") +
  scale_x_continuous(
    breaks = igra_spec_values,
    labels = c("0.90\n(conservative)", "0.93\n(T-SPOT base)", "0.96\n(QFT base)", "1.00\n(naive UB)")
  ) +
  scale_y_continuous(
    labels = function(x) paste0("\u00a3", format(round(x), big.mark = ",")),
    limits = c(0, max(spec_y_max, 28000)),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "IGRA specificity",
    y = "Sequential ICER (\u00a3/QALY)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9)
  )

ggsave("output/supplementary/igra_specificity_supplementary.tiff", p_igra_spec,
       width = 11, height = 7, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/igra_specificity_supplementary.tiff\n\n")

cat("\n")
cat("================================================================================\n")
cat("             ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS                         \n")
cat("================================================================================\n\n")

# Find a strategy object by its Excel row name
find_strategy <- function(strats, excel_nm) {
  idx <- which(sapply(strats, function(s) identical(s$excel_name, excel_nm)))
  if (length(idx) == 0) stop(paste("Strategy not found:", excel_nm))
  strats[[idx[1]]]
}

run_owsa <- function(strat_test, strat_ref = NULL) {

  # Restrict DSA to parameters with specified distributions (costs, QALYs,
  # transition probabilities). Fixed structural parameters (n_t, n_c, discount)
  # and detection rates are excluded from one-way variation here.
  dsa_params <- config %>%
    filter(distribution != "fixed", distribution != "", !is.na(distribution)) %>%
    filter(grepl("^(Cost_|qaly_|p_)", parameter))

  if (is.null(strat_ref)) strat_ref <- create_diagnostic_strategies()[["Passive_case_finding"]]
  base_result_test <- run_strategy(strat_test)
  base_result_ref  <- run_strategy(strat_ref)
  base_inc_cost    <- base_result_test$total_cost - base_result_ref$total_cost
  base_inc_qaly    <- base_result_test$total_qaly - base_result_ref$total_qaly
  base_icer        <- base_inc_cost / base_inc_qaly
  base_nmb         <- base_inc_qaly * 25000 - base_inc_cost  # NMB at £25,000/QALY WTP (NICE updated Dec 2025)

  # Parallelised one-way SA: each parameter's lo/hi pair is independent,
  # so all rows are computed in parallel with mclapply.
  n_cores_dsa <- max(1L, detectCores() - 1L)

  raw_results <- mclapply(seq_len(nrow(dsa_params)), function(r) {
    pname <- dsa_params$parameter[r]
    pval  <- as.numeric(dsa_params$value[r])
    pdist <- dsa_params$distribution[r]
    p1    <- as.numeric(dsa_params$dist_param1[r])
    p2    <- as.numeric(dsa_params$dist_param2[r])

    if (pdist == "beta") {
      lo <- qbeta(0.025, shape1 = p1, shape2 = p2)
      hi <- qbeta(0.975, shape1 = p1, shape2 = p2)
    } else if (pdist == "gamma") {
      lo <- qgamma(0.025, shape = p1, scale = p2)
      hi <- qgamma(0.975, shape = p1, scale = p2)
    } else {
      return(NULL)  # skip fixed / unsupported distributions
    }

    params_lo          <- paramsData
    params_lo[[pname]] <- lo
    res_test_lo  <- run_strategy(strat_test, params = params_lo)
    res_ref_lo   <- run_strategy(strat_ref,  params = params_lo)
    nmb_lo       <- (res_test_lo$total_qaly - res_ref_lo$total_qaly) * 25000 -
                    (res_test_lo$total_cost  - res_ref_lo$total_cost)

    params_hi          <- paramsData
    params_hi[[pname]] <- hi
    res_test_hi  <- run_strategy(strat_test, params = params_hi)
    res_ref_hi   <- run_strategy(strat_ref,  params = params_hi)
    nmb_hi       <- (res_test_hi$total_qaly - res_ref_hi$total_qaly) * 25000 -
                    (res_test_hi$total_cost  - res_ref_hi$total_cost)

    tibble(
      parameter   = pname,
      param_label = dsa_params$parameter_name[r],
      base_value  = pval,
      lo_value    = lo,
      hi_value    = hi,
      nmb_lo      = nmb_lo,
      nmb_hi      = nmb_hi,
      nmb_base    = base_nmb,
      nmb_range   = abs(nmb_hi - nmb_lo)
    )
  }, mc.cores = n_cores_dsa)

  dsa_df <- bind_rows(Filter(Negate(is.null), raw_results)) %>% arrange(desc(nmb_range))
  return(list(dsa_df = dsa_df, base_nmb = base_nmb, base_icer = base_icer))
}

# Helper: build and save a tornado diagram from a run_owsa() result.
# bar_col  = fill colour for bars; cap_col = endpoint marker colour.
# Short human-readable labels for tornado diagram y-axis
dsa_short_labels <- c(
  # LTBI states
  "QALY weight - Latent TB undiagnosed"              = "Utility: LTBI, undetected",
  "QALY weight - Uninfected"                         = "Utility: uninfected",
  "QALY weight - Latent TB not treated"              = "Utility: LTBI, untreated",
  "QALY weight - Latent TB lost to follow-up"        = "Utility: LTBI, lost to follow-up",
  "QALY weight - Latent TB treatment completed"      = "Utility: LTBI, treatment complete",
  "QALY weight - Latent TB diagnosed"                = "Utility: LTBI, diagnosed (awaiting tx)",
  "QALY weight - Latent TB treatment discontinued"   = "Utility: LTBI, treatment discontinued",
  # Active TB states
  "QALY weight - Active TB undiagnosed"              = "Utility: active TB, undetected",
  "QALY weight - Active TB under treatment"          = "Utility: active TB, on treatment",
  "QALY weight - Active TB treatment completed"      = "Utility: active TB, treatment complete",
  "QALY weight - Active TB treatment discontinued"   = "Utility: active TB, treatment discontinued",
  "QALY weight - Active TB diagnosed"                = "Utility: active TB, diagnosed (no tx)",
  "QALY weight - Active TB not treated"              = "Utility: active TB, untreated",
  "QALY weight - Active TB lost to follow-up"        = "Utility: active TB, lost to follow-up",
  # Reactivation
  "Prob: Latent undiagnosed to Active undiagnosed"   = "LTBI reactivation rate",
  "Prob: Latent not treated to Active undiagnosed"   = "Reactivation: untreated LTBI",
  "Prob: Latent LTFU to Active undiagnosed"          = "Reactivation: LTBI lost to follow-up",
  # LTBI treatment pathway
  "Prob: Latent diagnosed to treatment"              = "LTBI treatment initiation rate",
  "Prob: Latent diagnosed to not treated"            = "LTBI treatment non-initiation",
  "Prob: Latent treatment to discontinued"           = "LTBI treatment discontinuation",
  "Prob: Latent treatment to lost to follow-up"      = "LTBI treatment: lost to follow-up",
  # Active TB diagnosis / progression
  "Prob: Active undiagnosed to Dead"                 = "Active TB mortality (undetected)",
  "Prob: Active undiagnosed to diagnosed"            = "Active TB spontaneous detection rate",
  "Prob: Active diagnosed to treatment"              = "Active TB treatment initiation",
  "Prob: Active diagnosed to not treated"            = "Active TB: detected, no treatment",
  "Prob: Active diagnosed to Dead"                   = "Active TB mortality (detected, no tx)",
  "Prob: Active not treated to diagnosed"            = "Active TB: untreated → detected",
  "Prob: Active not treated to Dead"                 = "Active TB mortality (untreated)",
  # Active TB treatment
  "Prob: Active treatment to completed"              = "Active TB treatment completion",
  "Prob: Active treatment to discontinued"           = "Active TB treatment discontinuation",
  "Prob: Active treatment to lost to follow-up"      = "Active TB treatment: lost to follow-up",
  "Prob: Active treated to Dead"                     = "Active TB mortality (on treatment)",
  "Prob: Active discontinued to undiagnosed"         = "Active TB: discontinued → undetected",
  "Prob: Active discontinued to Dead"                = "Active TB mortality (discontinued)",
  # Cost
  "Cost - Active TB under treatment"                 = "Cost: active TB treatment (monthly)"
)

plot_tornado_diagram <- function(owsa_result, subtitle_str, bar_col, cap_col, out_file,
                                 top_n_override = NULL, base_size = 13, save = TRUE,
                                 label_width = 48, show_caption = FALSE, show_title = TRUE) {
  top_n  <- if (!is.null(top_n_override)) min(top_n_override, nrow(owsa_result$dsa_df))
             else min(15, nrow(owsa_result$dsa_df))
  t_df   <- owsa_result$dsa_df %>%
    slice_head(n = top_n) %>%
    mutate(
      param_short = dplyr::coalesce(dsa_short_labels[param_label], param_label),
      param_short = ifelse(nchar(param_short) > label_width,
                           paste0(substr(param_short, 1, label_width - 3), "..."),
                           param_short),
      param_short = make.unique(param_short, sep = " "),
      param_short = factor(param_short, levels = rev(param_short)),
      nmb_lo_dev  = nmb_lo - owsa_result$base_nmb,
      nmb_hi_dev  = nmb_hi - owsa_result$base_nmb,
      bar_left    = pmin(nmb_lo_dev, nmb_hi_dev),
      bar_right   = pmax(nmb_lo_dev, nmb_hi_dev)
    )
  p <- ggplot(t_df) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_segment(aes(y = param_short, yend = param_short,
                     x = bar_left, xend = bar_right),
                 linewidth = 8, color = bar_col, alpha = 0.8) +
    geom_point(aes(y = param_short, x = bar_left),  shape = "|", size = 3, color = cap_col) +
    geom_point(aes(y = param_short, x = bar_right), shape = "|", size = 3, color = cap_col) +
    scale_x_continuous(
      labels = function(x) {
        sign_str <- ifelse(x < 0, "\u2212\u00a3", "\u00a3")
        ax <- abs(x)
        ifelse(ax >= 1e6, paste0(sign_str, round(ax / 1e6, 1), "M"),
        ifelse(ax >= 1e3, paste0(sign_str, round(ax / 1e3, 0), "k"),
                          paste0(sign_str, round(ax))))
      },
      n.breaks = 7
    ) +
    theme_minimal(base_size = base_size) +
    labs(
      x        = "Change in NMB from base case (2023/24 \u00a3)",
      y        = "",
      title    = if (show_title) "One-way deterministic sensitivity analysis" else NULL,
      subtitle = subtitle_str,
      caption  = if (show_caption) "NMB = net monetary benefit; WTP = willingness-to-pay threshold (NICE \u00a325,000/QALY); DSA = deterministic sensitivity analysis.\nAll costs in 2023/24 GBP. Parameters varied individually across 95% plausible range; all others held at base-case values." else NULL
    ) +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = element_text(color = "grey40", size = base_size - 1),
      plot.caption       = element_text(color = "grey40", size = max(base_size - 3, 7),
                                        hjust = 0),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = base_size)
    )
  if (save) {
    ggsave(out_file, p, width = 14, height = 8, dpi = 300)
    cat(sprintf("Tornado diagram saved to %s\n", out_file))
  }
  invisible(p)
}

cat("Running one-way DSA (Parallel Sx+QFT (Ultra) vs No Screening)...\n")
owsa <- run_owsa(
  strat_test = find_strategy(strategies, "parallel allsx qft_ultra"),
  strat_ref  = strategies[["Passive_case_finding"]]
)
dsa_df <- owsa$dsa_df

# Report the 15 parameters with the largest impact on NMB
cat("\nTop 15 most influential parameters on NMB (£25k WTP):\n")
cat(sprintf("Base case NMB: £%s | Base case ICER: £%s\n\n",
            format(round(owsa$base_nmb), big.mark = ","),
            format(round(owsa$base_icer), big.mark = ",")))
cat(sprintf("%-45s %10s %12s %12s %12s\n",
            "Parameter", "Base", "NMB (Low)", "NMB (High)", "Range"))
cat(sprintf("%-45s %10s %12s %12s %12s\n",
            "---------------------------------------------", "----------",
            "------------", "------------", "------------"))
for (i in 1:min(15, nrow(dsa_df))) {
  cat(sprintf("%-45s %10.4f %11s %11s %11s\n",
              substr(dsa_df$param_label[i], 1, 45),
              dsa_df$base_value[i],
              paste0("£", format(round(dsa_df$nmb_lo[i]), big.mark = ",")),
              paste0("£", format(round(dsa_df$nmb_hi[i]), big.mark = ",")),
              paste0("£", format(round(dsa_df$nmb_range[i]), big.mark = ","))))
}

# Save full DSA results table for supplementary materials
write.csv(dsa_df, "output/csv/dsa_results.csv", row.names = FALSE)

# -------------------- Tornado diagram -----------------------------------------
# Displays the top 15 parameters ranked by their impact on NMB. Each bar shows
# how NMB changes when the parameter moves from its lower to upper bound.
# The width of each bar reflects the parameter's contribution to overall
# model uncertainty. Bars are expressed as deviation from the base case NMB
# so that the central reference line (x = 0) represents no change.
p_tornado_sx_qft <- plot_tornado_diagram(owsa,
  subtitle_str = "Parallel Sx+QFT (Ultra) vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#0099b4", cap_col = "#ad002a",
  out_file = NULL, save = FALSE)

# -------------------- Second tornado: Symptom screen+CXR vs No Screening
# Symptom screen+CXR is the second strategy on the efficient frontier
# (seq. ICER £3,822/QALY). This tornado characterises its structural uncertainty.
cat("Running one-way DSA (Symptom screen+CXR vs No Screening)...\n")
owsa_cxr <- run_owsa(
  strat_test = find_strategy(strategies, "anysx_cxr(any)"),
  strat_ref  = strategies[["Passive_case_finding"]]
)

cat(sprintf("Base case NMB (Symptom screen+CXR): £%s | Base case ICER: £%s\n\n",
            format(round(owsa_cxr$base_nmb), big.mark = ","),
            format(round(owsa_cxr$base_icer), big.mark = ",")))

p_tornado_cough_qft <- plot_tornado_diagram(owsa_cxr,
  subtitle_str = "Symptom screen+CXR vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#ed0000", cap_col = "#ad002a",
  out_file = NULL, save = FALSE)

# -------------------- Third tornado: Cough+CXR (TB sx) vs No Screening --------
# Cough+CXR (TB sx) is the lowest-cost strategy on the efficient frontier
# (seq. ICER ~£1,210/QALY vs No Screening; well within NICE £25,000 threshold).
# This tornado characterises the structural uncertainty around the cheapest
# CE strategy — the non-IGRA entry-level recommendation.
cat("Running one-way DSA (Cough+CXR (TB sx) vs No Screening)...\n")
owsa_cough_cxr <- run_owsa(
  strat_test = find_strategy(strategies, "anycough_cxr(TB)"),
  strat_ref  = strategies[["Passive_case_finding"]]
)

cat(sprintf("Base case NMB (Cough+CXR TB sx): £%s | Base case ICER: £%s\n\n",
            format(round(owsa_cough_cxr$base_nmb), big.mark = ","),
            format(round(owsa_cough_cxr$base_icer), big.mark = ",")))

p_tornado_cough_cxr <- plot_tornado_diagram(owsa_cough_cxr,
  subtitle_str = "Cough+CXR (TB sx) vs No Screening | NMB at \u00a325,000/QALY WTP | Top 15 parameters",
  bar_col = "#ed6060", cap_col = "#ad002a",
  out_file = NULL, save = FALSE)

# Save DSA results for Cough+CXR and Symptom screen+CXR (for panel regeneration)
write.csv(owsa_cough_cxr$dsa_df, "output/csv/dsa_results_coughcxr.csv",  row.names = FALSE)
write.csv(owsa_cxr$dsa_df,       "output/csv/dsa_results_symscrCXR.csv",  row.names = FALSE)

# -------------------- Tornado panel: all 3 frontier strategies ----------------
# Three panels, one per non-reference frontier strategy, ordered by ascending cost:
#   A. Cough+CXR (TB sx)        — seq. ICER £2,205/QALY (NHS)
#   B. Symptom screen+CXR       — seq. ICER £7,963/QALY (NHS)
#   C. Parallel Sx+QFT (Ultra)  — seq. ICER £11,324/QALY (NHS)
# Top 7 parameters each (top 8 too cramped at 3-panel width); shared caption.
# Panel construction delegated to run_plots_and_tables.R (reads CSVs written above).
source("scripts/run_plots_and_tables.R")

# =============================================================================
# IGRA TESTING PROGRAMME COST SENSITIVITY ANALYSIS — PARALLEL IGRA STRATEGIES
#
# Analytical one-way sensitivity on IGRA testing programme delivery cost.
# Programme delivery cost is a one-time lump-sum test cost that does not affect
# Markov transitions, so ICER changes linearly with programme cost:
#
#   ICER(pc) = [base_inc_cost + (pc − pc_base)] / base_inc_qaly
#
# where base_inc_cost and base_inc_qaly are cohort-level totals vs No Screening.
# No model re-runs are required.
#
# Range swept: £1–£10 per person (= £100,000–£1,000,000 per 100,000 migrants),
# covering the reviewer-requested £150k–£800k range with margin on each side.
# Base case points: cough strategies £3/pp (£300k), allsx strategies £5/pp (£500k).
# =============================================================================
cat("\n================================================================================\n")
cat("    IGRA TESTING PROGRAMME COST SENSITIVITY — PARALLEL IGRA STRATEGIES          \n")
cat("================================================================================\n\n")

parallel_keys_all <- names(strategies)[sapply(names(strategies), function(nm)
  !is.null(strategies[[nm]]$prog_cost_type))]

ref_sr       <- strategy_results[[which(sapply(strategy_results, function(x) x$strategy_name == "Passive case finding"))]]
ref_cost_tot <- ref_sr$total_cost
ref_qaly_tot <- ref_sr$total_qaly

prog_cost_pp_seq <- seq(0.5, 60, by = 0.5)  # £0.50–£60 per person (£50k–£6M per 100k)
# Extends beyond reviewer-requested £1.50–£8.00/pp to show that strategies remain
# cost-effective even at 7× the base case programme cost (breakeven ≈ £1,300/pp; far off-chart)

prog_cost_icer_df <- bind_rows(lapply(parallel_keys_all, function(k) {
  strat  <- strategies[[k]]
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strat$name))
  if (length(sr_idx) == 0) return(NULL)
  sr <- strategy_results[[sr_idx]]

  pc_base      <- strat$prog_cost_total         # base programme delivery cost in £ total
  base_ic_cost <- sr$total_cost - ref_cost_tot  # incremental cost vs No Screening at base programme cost
  base_ic_qaly <- sr$total_qaly - ref_qaly_tot  # incremental QALY (fixed — programme cost doesn't affect QALY)
  if (base_ic_qaly <= 0) return(NULL)           # skip dominated at base

  bind_rows(lapply(prog_cost_pp_seq, function(pc_pp) {
    pc          <- pc_pp * n_c               # total programme delivery cost £ for this scenario
    new_ic_cost <- base_ic_cost + (pc - pc_base)
    tibble(
      strategy        = strat$name,
      prog_cost_type  = strat$prog_cost_type,
      prog_cost_pp    = pc_pp,
      prog_cost_total = pc,
      icer            = new_ic_cost / base_ic_qaly
    )
  }))
}))

write.csv(prog_cost_icer_df, "output/csv/igra_programme_cost_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_programme_cost_sensitivity.csv\n")

# Print breakeven programme cost for each strategy (ICER = £25,000/QALY)
cat("\nBreakeven programme delivery cost (ICER = £25,000/QALY) for parallel IGRA strategies:\n")
cat(sprintf("  %-34s  %12s  %12s  %12s\n", "Strategy", "ICER (base)", "Breakeven £/pp", "Breakeven £/100k"))
cat(strrep("-", 80), "\n")
for (k in parallel_keys_all) {
  strat  <- strategies[[k]]
  df_k   <- prog_cost_icer_df %>% filter(strategy == strat$name)
  if (nrow(df_k) == 0) next
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strat$name))
  sr     <- strategy_results[[sr_idx]]
  base_icer_k    <- (sr$total_cost - ref_cost_tot) / (sr$total_qaly - ref_qaly_tot)
  base_ic_qaly_k <- sr$total_qaly - ref_qaly_tot
  # Breakeven: ICER(pc) = 25000 => pc = pc_base + (25000 - base_icer_k) * base_ic_qaly
  breakeven_pc   <- strat$prog_cost_total + (25000 - base_icer_k) * base_ic_qaly_k
  cat(sprintf("  %-34s  %11s  %12.2f  %12.0f\n",
              strat$name,
              paste0("\u00a3", format(round(base_icer_k), big.mark = ",")),
              breakeven_pc / n_c,
              breakeven_pc))
}
cat("\n")

# --- Plot: ICER vs programme delivery cost per person ---
# Focus on Ultra confirmatory strategies (the frontier ones)
plot_strats <- c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                 "Parallel Cough+T-STOP (Ultra)", "Parallel Sx+T-SPOT (Ultra)")

strat_cols_ov <- c(
  "Parallel Cough+QFT (Ultra)"    = "#7a8a99",  # steel gray (ext. dominated)
  "Parallel Sx+QFT (Ultra)"       = LANCET[7],  # dark red   (frontier)
  "Parallel Cough+T-STOP (Ultra)" = LANCET[8],  # gray       (non-frontier)
  "Parallel Sx+T-SPOT (Ultra)"    = LANCET[1]   # navy       (non-frontier)
)

prog_cost_plot_df <- prog_cost_icer_df %>%
  filter(strategy %in% plot_strats) %>%
  mutate(strategy = factor(strategy, levels = plot_strats))

# Base case annotation points (cough = £3/pp, allsx = £5/pp)
base_pts <- prog_cost_plot_df %>%
  mutate(is_base = (prog_cost_type == "cough" & abs(prog_cost_pp - 3) < 0.1) |
                   (prog_cost_type == "allsx" & abs(prog_cost_pp - 5) < 0.1)) %>%
  filter(is_base)

# Reviewer-requested range shading: £1.50–£8.00/pp
reviewer_lo <- 1.5
reviewer_hi <- 8.0

p_prog_cost <- ggplot(prog_cost_plot_df, aes(x = prog_cost_pp, y = icer,
                                              colour = strategy, group = strategy)) +
  # Shade reviewer-requested range
  annotate("rect", xmin = reviewer_lo, xmax = reviewer_hi, ymin = -Inf, ymax = Inf,
           fill = "#0099b4", alpha = 0.06) +
  geom_line(linewidth = 1.1) +
  # Base case vertical lines
  geom_vline(xintercept = 3, linetype = "dotdash", colour = "#0099b4", linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 5, linetype = "dotdash", colour = "#ed0000", linewidth = 0.5, alpha = 0.8) +
  # WTP thresholds
  geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey30", linewidth = 0.8) +
  geom_hline(yintercept = 35000, linetype = "dotted", colour = "grey50", linewidth = 0.7) +
  annotate("text", x = 59, y = 26200, label = "\u00a325,000/QALY (NICE lower)", hjust = 1, size = 3.0, colour = "grey30") +
  annotate("text", x = 59, y = 36200, label = "\u00a335,000/QALY (NICE upper)", hjust = 1, size = 3.0, colour = "grey50") +
  # Base case labels
  annotate("text", x = 3.5, y = 1200, label = "Base (cough)\n£3/pp", size = 2.6,
           colour = "#0099b4", hjust = 0) +
  annotate("text", x = 5.5, y = 1200, label = "Base (all-sx)\n£5/pp", size = 2.6,
           colour = "#ed0000", hjust = 0) +
  # Breakeven off-chart annotation
  annotate("text", x = 55, y = 3500,
           label = "Breakeven \u2248 £1,300/pp\n(off-chart; all strategies\nremain CE far beyond\nplausible programme cost range)",
           size = 2.5, colour = "grey40", hjust = 1) +
  # Base case dots
  geom_point(data = base_pts, aes(colour = strategy), size = 3, shape = 21,
             fill = "white", stroke = 1.5) +
  scale_colour_manual(values = strat_cols_ov, name = "Strategy") +
  scale_x_continuous(
    breaks = c(0.5, seq(5, 60, by = 5)),
    labels = function(x) paste0("\u00a3", x, "\n(", format(round(x * 100), big.mark = ","), "k/100k)")
  ) +
  scale_y_continuous(
    labels  = function(y) paste0("\u00a3", format(round(y), big.mark = ",")),
    limits  = c(0, 40000),
    breaks  = seq(0, 40000, by = 5000)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = NULL,
    x     = "IGRA testing programme delivery cost per person [\u00a3/pp] (equivalent per 100,000 migrants)",
    y     = "Sequential ICER (\u00a3/QALY)"
  ) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 2))


cat("\n================================================================================\n")

# =============================================================================
# IGRA UPTAKE SENSITIVITY (PARALLEL STRATEGIES)
#
# Base case assumes 100% of LTBI-positive migrants accept the IGRA test
# (LatentDiagnosed = n_ltbi × IGRA sensitivity). In practice uptake in UK
# migrant screening programmes is typically 60-80% (UKHSA 2024).
#
# This sensitivity sweeps uptake from 50% to 100% in 5% steps.
# For each level, LatentDiagnosed is scaled proportionally and
# LatentUndiagnosed compensated. Costs are held fixed (programme delivery cost is a
# programme-level cost independent of individual uptake). Pairwise ICERs
# vs No Screening are recomputed for all 8 parallel IGRA strategies.
# =============================================================================
cat("\n================================================================================\n")
cat("         IGRA UPTAKE SENSITIVITY (PARALLEL STRATEGIES)                          \n")
cat("================================================================================\n\n")

# Identify parallel IGRA strategies by presence of prog_cost_type field
parallel_keys <- names(strategies)[sapply(names(strategies), function(nm)
  !is.null(strategies[[nm]]$prog_cost_type))]

igra_uptake_vals <- seq(0.50, 1.00, by = 0.05)

# Run No Screening once as reference
res_no_screen_uptake <- run_strategy(strategies[["Passive_case_finding"]])

cat(sprintf("Sweeping IGRA uptake: %.0f%% to %.0f%% (%d steps) across %d parallel strategies\n\n",
            min(igra_uptake_vals)*100, max(igra_uptake_vals)*100,
            length(igra_uptake_vals), length(parallel_keys)))

uptake_rows <- lapply(igra_uptake_vals, function(u) {
  strat_results <- mclapply(parallel_keys, function(nm) {
    s <- strategies[[nm]]
    orig_ld <- s$init["LatentDiagnosed"]
    orig_lu <- s$init["LatentUndiagnosed"]
    # Scale from 100%-uptake value so SA always sweeps the full 0-100% range,
    # regardless of the base-case uptake level set in config (igra_uptake_base).
    full_ld      <- if (!is.null(s$ltbi_100pct)) s$ltbi_100pct else orig_ld
    n_ltbi_strat <- orig_ld + orig_lu   # conserves total LTBI
    s$init["LatentDiagnosed"]   <- full_ld * u
    s$init["LatentUndiagnosed"] <- n_ltbi_strat - full_ld * u
    # test_cost is held constant across uptake levels. programme_cost (£300k–£500k/100k)
    # is the dominant component and represents fixed programme infrastructure (staffing,
    # logistics, equipment). At lower uptake the IGRA reagent variable cost is overstated,
    # so ICERs are slightly inflated — this provides a pessimistic bound on cost-effectiveness.
    run_strategy(s)
  }, mc.cores = max(1L, detectCores() - 1L))

  lapply(strat_results, function(res) {
    inc_cost <- res$total_cost  - res_no_screen_uptake$total_cost
    inc_qaly <- res$total_qaly  - res_no_screen_uptake$total_qaly
    tibble(
      strategy = res$strategy_name,
      uptake   = u,
      icer_vs_noscreening = if (inc_qaly > 0) inc_cost / inc_qaly else NA_real_,
      cost_pp  = res$total_cost  / n_c,
      qaly_pp  = res$total_qaly  / n_c
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Print summary at base-case uptake levels (70%, 85%, 100%)
for (u_show in c(0.70, 0.85, 1.00)) {
  cat(sprintf("--- Uptake = %.0f%% ---\n", u_show * 100))
  sub <- uptake_rows %>%
    filter(abs(uptake - u_show) < 0.001) %>%
    arrange(icer_vs_noscreening)
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-38s  ICER vs No Screening: %s\n",
                sub$strategy[i],
                if (is.na(sub$icer_vs_noscreening[i])) "N/A"
                else paste0("\u00a3", format(round(sub$icer_vs_noscreening[i]), big.mark = ","))))
  }
  cat("\n")
}

write.csv(uptake_rows, "output/csv/igra_uptake_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/igra_uptake_sensitivity.csv\n\n")

# Plot: sequential ICER vs uptake — Parallel Sx+QFT (Ultra) only.
# Sequential ICER = (cost_parallel - cost_sx_cxr) / (qaly_parallel - qaly_sx_cxr).
# Symptom screen+CXR cost/QALY are invariant to IGRA uptake (no IGRA in pathway),
# so they are fixed at base-case values from icer_table.
sx_cxr_base <- icer_table[icer_table$strategy == "Symptom screen+CXR", ]
sx_cxr_cost <- sx_cxr_base$cost_per_person
sx_cxr_qaly <- sx_cxr_base$qaly_per_person

uptake_plot_df <- uptake_rows %>%
  filter(strategy == "Parallel Sx+QFT (Ultra)") %>%
  mutate(
    uptake_pct    = uptake * 100,
    seq_icer      = (cost_pp - sx_cxr_cost) / (qaly_pp - sx_cxr_qaly)
  )

# Label the three anchor points: 50% (low), 75% (PSA base), 100% (original base)
uptake_labels <- uptake_plot_df %>%
  filter(uptake_pct %in% c(50, 75, 100)) %>%
  mutate(label = paste0("\u00a3", format(round(seq_icer), big.mark = ",")))

p_uptake <- ggplot(uptake_plot_df,
                   aes(x = uptake_pct, y = seq_icer,
                       colour = "Parallel Sx+QFT (Ultra)")) +
  geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 51, y = 25800,
           label = "\u00a325,000/QALY (NICE threshold)", hjust = 0, size = 3.0,
           colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text(data = uptake_labels,
            aes(label = label),
            vjust = -0.9, size = 3.0, fontface = "bold",
            show.legend = FALSE) +
  geom_vline(xintercept = 75, linetype = "dotted", colour = "grey60", linewidth = 0.7) +
  annotate("text", x = 75, y = 1200,
           label = "PSA base case", hjust = 0.5, vjust = 0, size = 2.8, colour = "grey50") +
  scale_colour_manual(values = c("Parallel Sx+QFT (Ultra)" = LANCET[7]), name = "Strategy") +
  scale_x_continuous(breaks = seq(50, 100, by = 10),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = function(y) paste0("\u00a3", format(round(y), big.mark = ",")),
                     limits = c(0, 30000),
                     breaks = seq(0, 30000, by = 5000)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "IGRA test uptake (%)",
    y = "Sequential ICER (\u00a3/QALY)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9)
  )

ggsave("output/supplementary/igra_uptake_supplementary.tiff", p_uptake,
       width = 10, height = 7, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/igra_uptake_supplementary.tiff\n")

cat("\n================================================================================\n")

# =============================================================================
# LTBI TREATMENT COMPLETION SENSITIVITY
#
# Base case: p_LatentTreated_LatentCompleted = 0.253/mo
#   Source: Surey et al. 2021 (3HR arm, n=25, London UK): 76% completion / 3 months
#   With p_LTFU = 0.057/mo (also Surey 2021: 4/25 LFU over 12wk) and
#   p_discontinued = 0.02/mo, conditional completion = 0.253/(0.253+0.02+0.057) = 76.7%
#   → consistent with 76% trial primary outcome. Base case is now self-calibrating.
#
# Real-world lower bound: 0.185/mo
#   Source: UKHSA TB in England 2025 (Prevention: England 2024): 55.5% completion
#   in contacts / 3 months = 0.185/mo. Tests sensitivity to lower real-world adherence.
#   Note: UKHSA also reports 11.7% overall completion but this figure has large
#   missing data and is not suitable as a SA bound (see config.csv row 48 note).
#
# Direction of effect: lower completion → fewer people reach LatentCompleted →
#   slightly less QALY benefit from LTBI treatment → frontier ICERs increase.
#   Expected to remain <<£25k.
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI TREATMENT COMPLETION \u2014 TRIAL vs REAL-WORLD SENSITIVITY               \n")
cat("================================================================================\n\n")

p_ltbi_completion_rw <- as.numeric(config_list["p_ltbi_completion_rw"]) / 3

p_ltfu_base <- as.numeric(config_list["p_LatentTreated_LatentLtfu"])
p_disc_base <- as.numeric(config_list["p_LatentTreated_LatentDiscontinued"])
p_comp_base <- as.numeric(config_list["p_LatentTreated_LatentCompleted"])
cond_completion_base <- 100 * p_comp_base / (p_comp_base + p_disc_base + p_ltfu_base)

cat(sprintf("Base case:  p_LatentTreated_LatentCompleted = %.4f/mo  (Surey 2021: 76%% / 3mo)\n",
            p_comp_base))
cat(sprintf("            p_LatentTreated_LatentLtfu       = %.4f/mo  (Surey 2021: 4/25 LFU 3HR arm)\n",
            p_ltfu_base))
cat(sprintf("            Conditional completion in model  = %.1f%%  (target: 76%% trial rate)\n\n",
            cond_completion_base))
cat(sprintf("Sensitivity: p_LatentTreated_LatentCompleted = %.4f/mo  (UKHSA 2025: 55.5%% / 3mo)\n\n",
            p_ltbi_completion_rw))

params_ltbi_rw <- paramsData
params_ltbi_rw[["p_LatentTreated_LatentCompleted"]] <- p_ltbi_completion_rw

strategy_results_ltbi_rw <- mclapply(strategies, run_strategy,
                                      params   = params_ltbi_rw,
                                      mc.cores = max(1L, detectCores() - 1L))
icer_table_ltbi_rw <- calculate_icer(strategy_results_ltbi_rw)

cat("ICER TABLE \u2014 Sensitivity: real-world LTBI completion (55.5%; base case: 76%):\n")
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_ltbi_rw)) {
  icer_str     <- if (is.na(icer_table_ltbi_rw$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_table_ltbi_rw$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_ltbi_rw$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_table_ltbi_rw$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_ltbi_rw$strategy[i],
              paste0("\u00a3", format(round(icer_table_ltbi_rw$cost_per_person[i]), big.mark = ",")),
              icer_table_ltbi_rw$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_table_ltbi_rw$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_ltbi_rw$dominance[i]))
}
cat("\n")

# Focused frontier comparison (base vs real-world vs calibrated) — 3 scenarios
frontier_names_ltbi <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

cat(sprintf("FRONTIER COMPARISON \u2014 LTBI Completion: Base (~%.0f%% conditional) vs Real-World (55.5%% UKHSA 2025):\n",
            cond_completion_base))
cat(sprintf("  %-32s  %26s  %26s  %12s\n",
    "Strategy",
    sprintf("Base (~%.0f%%/0.253/mo)", cond_completion_base),
    "RW (55.5%/0.185/mo)",
    "Delta"))
cat(strrep("-", 104), "\n")
for (st in frontier_names_ltbi) {
  r_base <- icer_table[icer_table$strategy == st, ]
  r_rw   <- icer_table_ltbi_rw[icer_table_ltbi_rw$strategy == st, ]
  icer_b <- if (is.na(r_base$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_base$sequential_icer), big.mark = ","), "/QALY")
  icer_r <- if (is.na(r_rw$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_rw$sequential_icer), big.mark = ","), "/QALY")
  delta  <- if (is.na(r_base$sequential_icer) || is.na(r_rw$sequential_icer)) "\u2014" else
            paste0(ifelse(r_rw$sequential_icer - r_base$sequential_icer >= 0, "+", ""),
                   "\u00a3", format(round(r_rw$sequential_icer - r_base$sequential_icer),
                                    big.mark = ","))
  cat(sprintf("  %-32s  %26s  %26s  %12s\n", st, icer_b, icer_r, delta))
}
cat("\n================================================================================\n\n")

# Save combined CSV (2 scenarios: base trial rate, real-world UKHSA 2025)
icer_table_ltbi_combined <- bind_rows(
  icer_table %>%
    mutate(scenario = "trial_76pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c),
  icer_table_ltbi_rw %>%
    mutate(scenario = "realworld_55pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c)
)
write.csv(icer_table_ltbi_combined, "output/csv/ltbi_completion_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/ltbi_completion_sensitivity.csv\n\n")

# Frontier comparison plot: base -> real-world arrows
frontier_names <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

frontier_ltbi_base <- icer_table %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_ltbi_rw <- icer_table_ltbi_rw %>%
  filter(strategy %in% frontier_names) %>%
  select(strategy, q_rw = qaly, c_rw = cost) %>%
  mutate(q_rw = q_rw / n_c, c_rw = c_rw / n_c)

frontier_ltbi_arrows <- frontier_ltbi_base %>%
  left_join(frontier_ltbi_rw, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = frontier_names))

cond_label_base <- sprintf("Base case (~%.0f%% treatment completion)", cond_completion_base)
cond_label_rw   <- "Sensitivity: real-world completion (55.5%, UKHSA 2025)"

frontier_ltbi_long <- bind_rows(
  frontier_ltbi_arrows %>% transmute(strategy, q = q_base, c = c_base, cond = cond_label_base),
  frontier_ltbi_arrows %>% transmute(strategy, q = q_rw,   c = c_rw,   cond = cond_label_rw)
)

label_ltbi_df <- frontier_ltbi_arrows %>%
  transmute(strategy, q = q_rw, c = c_rw)

p_ltbi_compare <- ggplot() +
  geom_path(data = frontier_ltbi_arrows %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "grey55", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_ltbi_arrows %>% arrange(q_rw),
            aes(x = q_rw, y = c_rw),
            colour = LANCET[1], linetype = "solid", linewidth = 0.6, alpha = 0.8) +
  geom_segment(
    data = frontier_ltbi_arrows,
    aes(x = q_base, y = c_base, xend = q_rw, yend = c_rw),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = LANCET[1], linewidth = 0.50
  ) +
  geom_point(data = frontier_ltbi_long,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = label_ltbi_df,
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey20",
    nudge_y = 4, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = setNames(c("grey55", LANCET[1]),
                                        c(cond_label_base, cond_label_rw)), name = NULL) +
  scale_shape_manual(values  = setNames(c(16, 17),
                                        c(cond_label_base, cond_label_rw)), name = NULL) +
  theme_minimal(base_size = 12) +
  labs(
    x     = "QALYs per person",
    y     = "Cost per person (\u00a3)",
    title = NULL
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

ggsave("output/supplementary/ltbi_completion_supplementary.tiff", p_ltbi_compare,
       width = 11, height = 7, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/ltbi_completion_compare.tiff\n")

# =============================================================================
# SENSITIVITY ANALYSIS: LTBI TREATMENT EFFICACY
#
# Base: p_LatentCompleted_ActiveUndiagnosed = 0.000128/mo (HR 0.14; Berrocal-Almanza 2022)
#   Effective lifetime failure rate = 0.000128 / (0.08 + 0.000128) ≈ 0.16%
# SA:  NICE NG33 estimates 10-40% residual reactivation risk post-treatment.
#   Competing-risk formula: p_fail = target × p_exit / (1 - target)
#   where p_exit = p_LatentCompleted_Uninfected = 0.08/mo
#   20% failure: p_fail = 0.20 × 0.08 / 0.80 = 0.02/mo
#   40% failure: p_fail = 0.40 × 0.08 / 0.60 = 0.0533/mo
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI TREATMENT EFFICACY \u2014 RESIDUAL REACTIVATION RISK SENSITIVITY           \n")
cat("================================================================================\n\n")

p_lc_exit  <- 0.08          # p_LatentCompleted_Uninfected (from config)
p_lc_base  <- as.numeric(config_list["p_LatentCompleted_ActiveUndiagnosed"])
base_risk  <- p_lc_base / (p_lc_exit + p_lc_base)
cat(sprintf("Base case: p_fail = %.6f/mo  (lifetime failure rate = %.2f%%)\n",
            p_lc_base, base_risk * 100))

sa_ltbi_eff <- list(
  list(label = "20% residual risk (NICE NG33 lower)", target = 0.20),
  list(label = "40% residual risk (NICE NG33 upper)", target = 0.40)
)

ltbi_eff_rows <- list()
for (sa in sa_ltbi_eff) {
  p_fail <- sa$target * p_lc_exit / (1 - sa$target)
  actual_risk <- p_fail / (p_lc_exit + p_fail)
  cat(sprintf("SA: %-38s  p_fail = %.4f/mo  (residual risk = %.1f%%)\n",
              sa$label, p_fail, actual_risk * 100))

  params_eff <- paramsData
  params_eff[["p_LatentCompleted_ActiveUndiagnosed"]] <- p_fail

  res_eff  <- mclapply(strategies, run_strategy, params = params_eff,
                       mc.cores = max(1L, detectCores() - 1L))
  icer_eff <- calculate_icer(res_eff)

  par_row      <- icer_eff %>% filter(strategy == "Parallel Sx+QFT (Ultra)")
  par_dominance <- if (nrow(par_row) > 0) par_row$dominance else "not found"
  par_on_frontier <- par_dominance %in% c("ref", "non-dominated")
  frontier_strats <- icer_eff %>%
    filter(dominance %in% c("ref", "non-dominated")) %>%
    arrange(cost_per_person) %>%
    pull(strategy)
  ltbi_eff_rows[[length(ltbi_eff_rows) + 1]] <- data.frame(
    scenario          = sa$label,
    p_fail            = round(p_fail, 6),
    residual_risk     = round(actual_risk, 4),
    parallel_on_frontier = par_on_frontier,
    parallel_dominance   = par_dominance,
    frontier_strategies  = paste(frontier_strats, collapse = " | ")
  )
}

ltbi_eff_summary <- bind_rows(ltbi_eff_rows)
write.csv(ltbi_eff_summary, "output/csv/ltbi_efficacy_sensitivity.csv", row.names = FALSE)
cat("\nSaved: output/csv/ltbi_efficacy_sensitivity.csv\n")
for (i in seq_len(nrow(ltbi_eff_summary))) {
  if (ltbi_eff_summary$parallel_on_frontier[i]) {
    cat(sprintf("  %s:  Parallel Sx+QFT on frontier\n", ltbi_eff_summary$scenario[i]))
  } else {
    cat(sprintf("  %s:  Parallel Sx+QFT %s (off frontier)\n",
                ltbi_eff_summary$scenario[i], ltbi_eff_summary$parallel_dominance[i]))
    cat(sprintf("    Frontier: %s\n", ltbi_eff_summary$frontier_strategies[i]))
  }
}
cat("\n")

# =============================================================================
# SENSITIVITY ANALYSIS: LTBI PREVALENCE
#
# Base case: 17.8% (Berrocal-Almanza et al. 2022, IGRA-screened UK migrants from
#   high-burden countries; pooled estimate from 11 UK studies)
# Sensitivity: 15.1% (UKHSA TB in England 2025 Prevention report, post-entry
#   migrant screening; more recent national estimate)
#
# Direction of effect: lower LTBI prevalence → fewer LTBI detected →
#   smaller QALY gains from LTBI treatment → frontier ICERs increase
#   (parallel IGRA strategies marginally less attractive, but remain <<£25k)
#
# Init vector adjustment:
#   - Sequential strategies: LatentDiagnosed scaled proportionally (×15100/17800);
#     LatentUndiagnosed = 15100 − LatentDiagnosed; Uninfected = n_c − 1000 − 15100
#   - Parallel IGRA strategies: LatentDiagnosed = 15100 × igra_sensitivity;
#     LatentUndiagnosed = 15100 − LatentDiagnosed; Uninfected = 83900
#   - Passive_case_finding: LatentUndiagnosed = 15100; Uninfected = 83900
#   - ActiveUndiagnosed: unchanged at 1% (prev_active; separate parameter)
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI PREVALENCE SENSITIVITY \u2014 17.8% (base) vs 15.1% (UKHSA 2025)          \n")
cat("================================================================================\n\n")

prev_ltbi_ukhsa    <- as.numeric(config_list["prev_ltbi_ukhsa"])
n_ltbi_ukhsa       <- n_c * prev_ltbi_ukhsa                         # 15,100 per 100k
n_uninfected_ukhsa <- n_c * (1 - prev_ltbi_ukhsa - prev_active)     # 83,900 per 100k
scale_ltbi_prev    <- prev_ltbi_ukhsa / prev_ltbi                    # 15.1 / 17.8 = 0.8483

cat(sprintf("Base:        prev_ltbi = %.1f%% (Berrocal-Almanza 2022) \u2192 %g LTBI per 100k\n",
            prev_ltbi * 100, n_ltbi_global))
cat(sprintf("Sensitivity: prev_ltbi = %.1f%% (UKHSA 2025)             \u2192 %g LTBI per 100k\n\n",
            prev_ltbi_ukhsa * 100, n_ltbi_ukhsa))

strategies_ltbi_prev <- lapply(strategies, function(s) {
  if (is.null(s$excel_name) || is.na(s$excel_name)) {
    # Passive_case_finding: no LTBI detected
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  } else if (!is.null(s$prog_cost_type)) {
    # Parallel IGRA: detection = n_ltbi_ukhsa × igra_sensitivity
    sens <- if (grepl("tspt", s$excel_name, ignore.case = TRUE)) igra_sens_tspt else igra_sens_qft
    ltbi_det_new <- n_ltbi_ukhsa * sens
    s$init["LatentDiagnosed"]   <- ltbi_det_new
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa - ltbi_det_new
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  } else {
    # Sequential strategies: scale LatentDiagnosed proportionally
    ltbi_det_new <- s$init["LatentDiagnosed"] * scale_ltbi_prev
    s$init["LatentDiagnosed"]   <- ltbi_det_new
    s$init["LatentUndiagnosed"] <- n_ltbi_ukhsa - ltbi_det_new
    s$init["Uninfected"]        <- n_uninfected_ukhsa
  }
  s
})

strategy_results_ltbi_prev <- mclapply(strategies_ltbi_prev, run_strategy,
                                        params   = paramsData,
                                        mc.cores = max(1L, detectCores() - 1L))
icer_table_ltbi_prev <- calculate_icer(strategy_results_ltbi_prev)

cat("ICER TABLE \u2014 Sensitivity: LTBI prevalence 15.1% (UKHSA 2025):\n")
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_table_ltbi_prev)) {
  icer_str     <- if (is.na(icer_table_ltbi_prev$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_table_ltbi_prev$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table_ltbi_prev$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_table_ltbi_prev$sequential_icer[i]),
                                               big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_table_ltbi_prev$strategy[i],
              paste0("\u00a3", format(round(icer_table_ltbi_prev$cost_per_person[i]), big.mark = ",")),
              icer_table_ltbi_prev$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_table_ltbi_prev$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_table_ltbi_prev$dominance[i]))
}
cat("\n")

# Frontier comparison: base (17.8%) vs UKHSA 2025 (15.1%)
frontier_names_prev <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

cat("FRONTIER COMPARISON \u2014 LTBI Prevalence: Base (17.8%) vs UKHSA 2025 (15.1%):\n")
cat(sprintf("  %-32s  %18s  %24s  %12s\n",
    "Strategy", "Base ICER (17.8%)", "UKHSA 2025 ICER (15.1%)", "Delta"))
cat(strrep("-", 96), "\n")
for (st in frontier_names_prev) {
  r_base <- icer_table[icer_table$strategy == st, ]
  r_prev <- icer_table_ltbi_prev[icer_table_ltbi_prev$strategy == st, ]
  icer_b <- if (is.na(r_base$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_base$sequential_icer), big.mark = ","), "/QALY")
  icer_p <- if (is.na(r_prev$sequential_icer)) "Ref" else
            paste0("\u00a3", format(round(r_prev$sequential_icer), big.mark = ","), "/QALY")
  delta  <- if (is.na(r_base$sequential_icer) || is.na(r_prev$sequential_icer)) "\u2014" else
            paste0(ifelse(r_prev$sequential_icer - r_base$sequential_icer >= 0, "+", ""),
                   "\u00a3", format(round(r_prev$sequential_icer - r_base$sequential_icer),
                                    big.mark = ","))
  cat(sprintf("  %-32s  %18s  %24s  %12s\n", st, icer_b, icer_p, delta))
}
cat("\n================================================================================\n\n")

# Save CSV
icer_table_prev_combined <- bind_rows(
  icer_table %>%
    mutate(scenario = "base_17.8pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c),
  icer_table_ltbi_prev %>%
    mutate(scenario = "ukhsa2025_15.1pct",
           cost_per_person = cost / n_c,
           qaly_per_person = qaly / n_c)
)
write.csv(icer_table_prev_combined, "output/csv/ltbi_prevalence_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/ltbi_prevalence_sensitivity.csv\n")

# LTBI prevalence sensitivity plot: base (17.8%) vs UKHSA 2025 (15.1%)
frontier_ltbi_prev_base <- icer_table %>%
  filter(strategy %in% frontier_names_prev) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_ltbi_prev_new <- icer_table_ltbi_prev %>%
  filter(strategy %in% frontier_names_prev) %>%
  select(strategy, q_new = qaly, c_new = cost) %>%
  mutate(q_new = q_new / n_c, c_new = c_new / n_c)

frontier_prev_arrows <- frontier_ltbi_prev_base %>%
  left_join(frontier_ltbi_prev_new, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = frontier_names_prev))

lbl_base_prev <- "Base case (17.8% LTBI prevalence)"
lbl_new_prev  <- "Sensitivity: UKHSA 2025 (15.1% LTBI prevalence)"

frontier_prev_long <- bind_rows(
  frontier_prev_arrows %>% transmute(strategy, q = q_base, c = c_base, cond = lbl_base_prev),
  frontier_prev_arrows %>% transmute(strategy, q = q_new,  c = c_new,  cond = lbl_new_prev)
)

label_prev_df <- frontier_prev_arrows %>%
  transmute(strategy, q = q_new, c = c_new)

p_ltbi_prev <- ggplot() +
  geom_path(data = frontier_prev_arrows %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = "grey55", linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_prev_arrows %>% arrange(q_new),
            aes(x = q_new, y = c_new),
            colour = LANCET[1], linetype = "solid", linewidth = 0.6, alpha = 0.8) +
  geom_segment(
    data = frontier_prev_arrows,
    aes(x = q_base, y = c_base, xend = q_new, yend = c_new),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = LANCET[1], linewidth = 0.50
  ) +
  geom_point(data = frontier_prev_long,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = label_prev_df,
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey20",
    nudge_y = 4, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = setNames(c("grey55", LANCET[1]),
                                        c(lbl_base_prev, lbl_new_prev)), name = NULL) +
  scale_shape_manual(values  = setNames(c(16, 17),
                                        c(lbl_base_prev, lbl_new_prev)), name = NULL) +
  scale_y_continuous(labels = function(x) paste0("\u00a3", x)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "QALYs per person",
    y = "Cost per person (\u00a3)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

ggsave("output/supplementary/ltbi_prevalence_supplementary.tiff", p_ltbi_prev,
       width = 11, height = 7, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/ltbi_prevalence_sensitivity.tiff\n")

# Combined 2-panel figure: LTBI completion (A) / LTBI prevalence (B)
# p_ltbi_compare = LTBI treatment completion (SF5); p_ltbi_prev = LTBI prevalence (SF6)
p_ltbi_combined <- p_ltbi_compare / p_ltbi_prev + plot_annotation(tag_levels = "A")
ggsave("output/supplementary/ltbi_sa_combined_supplementary.tiff", p_ltbi_combined,
       width = 10, height = 14, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/ltbi_sa_combined_supplementary.tiff\n")

# =============================================================================
# LTBI TREATMENT UTILITY DECREMENT SENSITIVITY ANALYSIS
#
# Base case: qaly_LatentTreated = 0.92 (no decrement beyond SAEs; standard UK CEA
#   practice; consistent with Pareek 2011, Zenner 2022).
# SA: Apply annual utility decrement of 0.0133 during LTBI treatment months.
#   Source: Bauer et al. 2015 (Qual Life Res 24:1337-1349) — longitudinal cohort;
#   predominantly 9H isoniazid; decrement relative to healthy controls.
#   Used identically by Dale et al. 2022 (Am J Epidemiol 191:255-270) who found
#   this decrement dominated all onshore LTBI strategies in Australian migrants.
#
# Implementation: qaly_LatentTreated = 0.92 - 0.0133 = 0.9067.
#   All other QALY states unchanged.
#   Note: 3HR is the current UK standard regimen (~3 months); Bauer 2015 used 9H
#   isoniazid. The decrement here is therefore conservative (applied over fewer cycles
#   in practice). The SA tests whether the qualitative finding from Dale holds in
#   a UK context with more favourable WTP thresholds and transmission savings.
# =============================================================================
cat("\n================================================================================\n")
cat("    LTBI TREATMENT UTILITY DECREMENT (Bauer 2015 / Dale 2022)                   \n")
cat("================================================================================\n\n")

qaly_lt_base <- as.numeric(config_list["qaly_LatentTreated"])   # 0.92
decrement_bauer <- 0.0133                                        # Bauer 2015 annual decrement
qaly_lt_sa   <- qaly_lt_base - decrement_bauer                  # 0.9067

cat(sprintf("Base:  qaly_LatentTreated = %.4f (no decrement; Kind 1998 population norm)\n",
            qaly_lt_base))
cat(sprintf("SA:    qaly_LatentTreated = %.4f (decrement = -%.4f; Bauer 2015 / Dale 2022)\n\n",
            qaly_lt_sa, decrement_bauer))

params_ltu        <- paramsData
params_ltu[["qaly_LatentTreated"]] <- qaly_lt_sa

res_ltu  <- mclapply(strategies, run_strategy, params = params_ltu,
                     mc.cores = max(1L, detectCores() - 1L))
icer_ltu <- calculate_icer(res_ltu)

frontier_ltu <- icer_ltu %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost) %>%
  pull(strategy)

par_row_ltu <- icer_ltu %>% filter(strategy == "Parallel Sx+QFT (Ultra)")
par_on_frontier_ltu  <- nrow(par_row_ltu) > 0 &&
  par_row_ltu$dominance %in% c("ref", "non-dominated")

cat("Frontier strategies (utility decrement SA):\n")
for (st in frontier_ltu) {
  r <- icer_ltu %>% filter(strategy == st)
  seq_str <- if (is.na(r$sequential_icer)) "Ref" else
    paste0("£", format(round(r$sequential_icer), big.mark = ","), "/QALY")
  cat(sprintf("  %-32s  seq. ICER %s\n", st, seq_str))
}
cat("\n")

if (par_on_frontier_ltu) {
  cat("  Parallel Sx+QFT (Ultra): ON frontier\n")
} else {
  cat(sprintf("  Parallel Sx+QFT (Ultra): OFF frontier  [dominance: %s]\n",
              if (nrow(par_row_ltu) > 0) par_row_ltu$dominance else "not found"))
}

# Frontier comparison: base vs utility decrement SA
cat("\nFRONTIER COMPARISON — Base vs LTBI treatment utility decrement (Bauer 2015):\n")
cat(sprintf("  %-32s  %22s  %26s\n",
    "Strategy", "Base ICER", "Utility decrement ICER"))
cat(strrep("-", 88), "\n")
all_frontier_strats <- union(
  icer_table %>% filter(dominance %in% c("ref", "non-dominated")) %>% pull(strategy),
  frontier_ltu
)
for (st in all_frontier_strats) {
  r_base <- icer_table %>% filter(strategy == st)
  r_ltu  <- icer_ltu   %>% filter(strategy == st)
  icer_b <- if (nrow(r_base) == 0 || is.na(r_base$sequential_icer)) "—"
            else paste0("£", format(round(r_base$sequential_icer), big.mark = ","), "/QALY")
  icer_l <- if (nrow(r_ltu) == 0 || is.na(r_ltu$sequential_icer)) {
    if (nrow(r_ltu) > 0) paste0("[", r_ltu$dominance, "]") else "—"
  } else paste0("£", format(round(r_ltu$sequential_icer), big.mark = ","), "/QALY")
  cat(sprintf("  %-32s  %22s  %26s\n", st, icer_b, icer_l))
}
cat("\n")

# Save CSV
ltu_combined <- bind_rows(
  icer_table %>% mutate(scenario = "base_no_decrement",
                        cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
  icer_ltu   %>% mutate(scenario = "ltbi_tx_utility_decrement_0.0133",
                        cost_per_person = cost / n_c, qaly_per_person = qaly / n_c)
)
write.csv(ltu_combined, "output/csv/ltbi_tx_utility_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/ltbi_tx_utility_sensitivity.csv\n\n")

# =============================================================================
# ACTIVE TB PREVALENCE SENSITIVITY ANALYSIS
# Base: 1% (1,000/100k; Zenner et al. 2025 ERJ decision tree, high-risk migrants)
# SA1:  0.215% (215/100k; Chen et al. 2025 pooled meta-analysis: 40M+ migrants)
# SA2:  0.44%  (440/100k; upper bound for refugee populations; Chen et al. 2025)
# Method: TP and FN scaled proportionally (fixed detection rate assumption).
# Rationale: entry prevalence does not alter test accuracy; it scales the absolute
#   number of cases detected and missed proportionally.
# =============================================================================
cat("\n================================================================================\n")
cat("         ACTIVE TB PREVALENCE \u2014 ENTRY PREVALENCE SENSITIVITY                    \n")
cat("================================================================================\n\n")
cat(sprintf("Base case:   prev_active = %.1f%% (%d/100k)  [Zenner et al. 2025 ERJ]\n",
            prev_active * 100, round(n_c * prev_active)))
cat("Sensitivity: 215/100k  [Chen et al. 2025 pooled meta-analysis: 40M+ migrants]\n")
cat("             440/100k  [Chen et al. 2025: high-risk subgroup estimate]\n\n")

active_tb_sa_scenarios <- list(
  list(label = "base_1pct",      prev = prev_active,                                      desc = "1%   (1,000/100k; base case)"),
  list(label = "pooled_0215pct", prev = as.numeric(config_list["prev_active_low"]),        desc = "0.215% (215/100k; Chen et al. 2025 pooled)"),
  list(label = "highrisk_044pct", prev = as.numeric(config_list["prev_active_high"]),      desc = "0.44%  (440/100k; high-risk subgroup; Chen et al. 2025)")
)

active_tb_sa_results <- list()

for (sc in active_tb_sa_scenarios) {
  scale_f <- sc$prev / prev_active
  strats_atb <- lapply(strategies, function(s) {
    s2 <- s
    if (!is.null(s2$base_tp) && !is.na(s2$base_tp)) {
      new_tp <- s2$base_tp * scale_f
      new_fn <- n_c * sc$prev - new_tp
      s2$init["ActiveDiagnosed"]   <- max(0, new_tp)
      s2$init["ActiveUndiagnosed"] <- max(0, new_fn)
    } else {
      # Passive_case_finding
      s2$init["ActiveDiagnosed"]   <- 0
      s2$init["ActiveUndiagnosed"] <- n_c * sc$prev
    }
    s2
  })
  res_atb <- mclapply(strats_atb, run_strategy,
                      params   = paramsData,
                      mc.cores = max(1L, detectCores() - 1L))
  icer_atb <- calculate_icer(res_atb)
  active_tb_sa_results[[sc$label]] <- list(icer = icer_atb, desc = sc$desc, prev = sc$prev)
}

# Print frontier comparison
cat("FRONTIER COMPARISON \u2014 Active TB Entry Prevalence:\n")
frontier_strategies <- active_tb_sa_results[["base_1pct"]]$icer %>%
  filter(dominance == "non-dominated") %>%
  pull(strategy)
frontier_strategies <- c("Passive case finding", frontier_strategies[frontier_strategies != "Passive case finding"])

header <- sprintf("  %-40s  %12s  %12s  %12s",
                  "Strategy",
                  "Base 1%",
                  "Pooled 0.215%",
                  "Refugee 0.44%")
cat(header, "\n")
cat(strrep("-", nchar(header)), "\n")

for (strat in frontier_strategies) {
  icers <- sapply(active_tb_sa_results, function(sc) {
    r <- sc$icer %>% filter(strategy == strat)
    if (nrow(r) == 0) return("N/A")
    si <- r$sequential_icer
    if (is.na(si) || si == Inf) return("dominated")
    sprintf("\u00a3%s/QALY", format(round(si), big.mark = ","))
  })
  cat(sprintf("  %-40s  %12s  %12s  %12s\n", strat, icers[1], icers[2], icers[3]))
}

# Save CSV
atb_sa_all <- bind_rows(lapply(names(active_tb_sa_results), function(nm) {
  active_tb_sa_results[[nm]]$icer %>%
    mutate(scenario       = nm,
           active_tb_prev = active_tb_sa_results[[nm]]$prev)
}))
write.csv(atb_sa_all, "output/csv/active_tb_prev_sensitivity.csv", row.names = FALSE)
cat("\nSaved: output/csv/active_tb_prev_sensitivity.csv\n")

# Comparison plot: ICER for frontier strategies across prevalence scenarios.
# Includes Symptom screen+CXR at all 3 prevalence values.
# At 0.215% and 0.44%: Symptom screen+CXR is extendedly dominated — shown as
# grey bar with annotation. At 1% (base): shown in full colour on frontier.
atb_sa_strats_plot <- c("Cough+CXR (TB sx)", "Symptom screen+CXR",
                         "Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)")

atb_plot_df <- bind_rows(lapply(names(active_tb_sa_results), function(nm) {
  sc <- active_tb_sa_results[[nm]]
  sc$icer %>%
    filter(strategy %in% atb_sa_strats_plot) %>%
    mutate(
      scenario = nm,
      prev_pct = sc$prev * 100,
      label    = sc$desc
    )
})) %>%
  filter(sequential_icer != Inf | is.na(sequential_icer)) %>%
  mutate(
    is_dominated = is.na(sequential_icer) | dominance == "extendedly dominated",
    icer_plot    = ifelse(is_dominated, NA_real_, sequential_icer)
  )

if (nrow(atb_plot_df) > 0) {
  # x-axis order: ascending prevalence (0.215 < 0.44 < 1)
  atb_plot_df$prev_label <- factor(
    sprintf("%.3g%%", atb_plot_df$prev_pct),
    levels = c("0.215%", "0.44%", "1%")
  )

  # Colour palette: dominated Symptom screen+CXR shown as grey; others by strategy colour
  atb_fill_vals <- c(
    "Cough+CXR (TB sx)"          = LANCET[4],
    "Symptom screen+CXR"         = LANCET[1],
    "Parallel Cough+QFT (Ultra)" = "#7a8a99",
    "Parallel Sx+QFT (Ultra)"    = LANCET[7]
  )

  # Data for "extendedly dominated" label bars (Symptom screen+CXR at low prevalence)
  atb_dom_label_df <- atb_plot_df %>%
    filter(is_dominated, strategy == "Symptom screen+CXR")

  # Value labels for 0.215% scenario bars (non-dominated only)
  atb_val_labels_215 <- atb_plot_df %>%
    filter(prev_label == "0.215%", !is_dominated, !is.na(icer_plot))

  p_atb_sa <- ggplot(atb_plot_df %>% filter(!is_dominated),
                     aes(x = prev_label, y = icer_plot, fill = strategy)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    # Extendedly dominated bars shown as grey
    geom_col(data = atb_dom_label_df,
             aes(x = prev_label, y = 1000, fill = strategy),
             position = "dodge", width = 0.7, alpha = 0.35,
             inherit.aes = FALSE) +
    geom_text(data = atb_dom_label_df,
              aes(x = prev_label, y = 600, label = "Extendedly\ndominated", group = strategy),
              position = position_dodge(width = 0.7),
              size = 2.5, colour = "grey50", fontface = "italic",
              vjust = 0, inherit.aes = FALSE) +
    # Value labels on 0.215% bars
    geom_text(data = atb_val_labels_215,
              aes(x = prev_label, y = icer_plot,
                  label = paste0("\u00a3", format(round(icer_plot), big.mark = ",")),
                  group = strategy),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 2.8, colour = "grey20",
              inherit.aes = FALSE) +
    geom_hline(yintercept = 25000, linetype = "dashed", colour = "grey30", linewidth = 0.8) +
    annotate("text", x = 3.45, y = 26500, label = "\u00a325,000/QALY (NICE threshold)",
             colour = "grey30", size = 3.0, hjust = 1) +
    scale_y_continuous(labels = scales::label_comma(prefix = "\u00a3"),
                       expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = atb_fill_vals) +
    labs(
      title = NULL,
      x     = "Active TB entry prevalence (%)",
      y     = "Sequential ICER (\u00a3/QALY)",
      fill  = "Strategy"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave("output/supplementary/active_tb_prevalence_supplementary.tiff", p_atb_sa,
         width = 10, height = 6, dpi = 300, compression = "lzw")
  cat("Saved: output/supplementary/active_tb_prev_sensitivity.tiff\n")
}

cat("\n================================================================================\n")

# =============================================================================
# COST_ACTIVEUNDIAGNOSED SENSITIVITY ANALYSIS
#
# Base case: Cost_ActiveUndiagnosed = £0/month (conservative; standard UK TB CEA
#   simplification — Jit 2011, Green 2025, Zenner 2022 all use £0 for pre-diagnosis
#   state costs).
#
# Sensitivity: Cost_ActiveUndiagnosed = £150 one-time cost per person
#   (spread as a monthly cost over mean dwell time ~6.7 months = ~£22/month)
#   Evidence: Jones et al. BJGP 2007 — 2-4 GP visits before TB diagnosis in UK;
#   PSSRU Unit Costs 2024: GP face-to-face £49 × 3 visits = £147 ≈ £150.
#   Aldridge et al. 2019 (London): 45% GP-first, 35% A&E-first routes to diagnosis.
#
# Direction: adding cost to ActiveUndiagnosed raises Passive case finding costs
#   more than screening arms (more person-months undiagnosed) → incremental costs
#   for screening strategies fall → ICERs decrease (conservative estimate direction).
#   Effect expected to be small (~5-8% ICER reduction based on analytic calculation).
#
# Source: config parameter Cost_ActiveUndiagnosed_scenario = 150
# =============================================================================
cat("\n================================================================================\n")
cat("    COST_ACTIVEUNDIAGNOSED SENSITIVITY (base £0 vs scenario £150/person)        \n")
cat("================================================================================\n\n")

cost_au_scenario <- as.numeric(config_list["Cost_ActiveUndiagnosed_scenario"])
cost_au_base     <- as.numeric(config_list["Cost_ActiveUndiagnosed"])

cat(sprintf("Base case:   Cost_ActiveUndiagnosed = \u00a3%.0f/month\n", cost_au_base))
cat(sprintf("Sensitivity: Cost_ActiveUndiagnosed = \u00a3%.0f one-time cost per undiagnosed person\n",
            cost_au_scenario))
cat(sprintf("             (spread over mean dwell ~6.7 months = \u00a3%.1f/month)\n\n",
            cost_au_scenario / (1 / (as.numeric(config_list["p_ActiveUndiagnosed_ActiveDiagnosed"]) +
                                    as.numeric(config_list["p_ActiveUndiagnosed_Dead"])))))

# Model the scenario cost as a monthly state cost (total cost / dwell time)
mean_dwell_au <- 1 / (as.numeric(config_list["p_ActiveUndiagnosed_ActiveDiagnosed"]) +
                      as.numeric(config_list["p_ActiveUndiagnosed_Dead"]))
monthly_cost_au_scenario <- cost_au_scenario / mean_dwell_au

params_cost_au <- paramsData
params_cost_au[["Cost_ActiveUndiagnosed"]] <- monthly_cost_au_scenario

strategy_results_cost_au <- mclapply(strategies, run_strategy,
                                      params   = params_cost_au,
                                      mc.cores = max(1L, detectCores() - 1L))

icer_cost_au <- calculate_icer(strategy_results_cost_au)

frontier_cost_au <- icer_cost_au %>%
  filter(dominance %in% c("ref", "non-dominated")) %>%
  arrange(cost_per_person)

cat("Frontier strategies — Cost_ActiveUndiagnosed scenario (£150/person):\n")
cat(sprintf("  %-40s  %10s  %12s\n", "Strategy", "Cost/pp", "Seq ICER"))
for (i in seq_len(nrow(frontier_cost_au))) {
  r <- frontier_cost_au[i, ]
  icer_str <- if (is.na(r$sequential_icer)) "(ref)" else
    paste0("\u00a3", format(round(r$sequential_icer), big.mark = ","))
  cat(sprintf("  %-40s  \u00a3%8.0f  %12s\n",
              r$strategy, r$cost_per_person, icer_str))
}

# Compare ICERs base vs scenario
base_front <- icer_table %>%
  filter(dominance %in% c("ref", "non-dominated"), strategy != "Passive case finding") %>%
  select(strategy, sequential_icer_base = sequential_icer)
sa_front <- frontier_cost_au %>%
  filter(strategy != "Passive case finding") %>%
  select(strategy, sequential_icer_sa = sequential_icer)

compare_df <- left_join(base_front, sa_front, by = "strategy") %>%
  mutate(
    icer_change    = sequential_icer_sa - sequential_icer_base,
    pct_change     = 100 * icer_change / sequential_icer_base
  )

cat("\nICER comparison — base (£0) vs scenario (£150/person undiagnosed cost):\n")
cat(sprintf("  %-40s  %12s  %12s  %8s\n", "Strategy", "Base ICER", "SA ICER", "Change"))
for (i in seq_len(nrow(compare_df))) {
  r <- compare_df[i, ]
  if (is.na(r$sequential_icer_base) || is.na(r$sequential_icer_sa)) next
  cat(sprintf("  %-40s  \u00a3%10s  \u00a3%10s  %+6.1f%%\n",
              r$strategy,
              format(round(r$sequential_icer_base), big.mark = ","),
              format(round(r$sequential_icer_sa),   big.mark = ","),
              r$pct_change))
}

# Save CSV
cost_au_sa_out <- icer_cost_au %>%
  mutate(scenario = sprintf("undiag_cost_£%.0f_per_person", cost_au_scenario))
write.csv(cost_au_sa_out, "output/csv/cost_activeundiag_sensitivity.csv", row.names = FALSE)
cat("\nSaved: output/csv/cost_activeundiag_sensitivity.csv\n")

cat("\n================================================================================\n")

# =============================================================================
# POST-TB SEQUELAE SENSITIVITY ANALYSIS
#
# Base case: qaly_ActiveCompleted = 0.92 (Kind 1998 age-16-35 population norm;
#   assumes full recovery after TB treatment).
# SA: three values across Ranzani et al. 2021 SR range (0.88–0.91): 0.88, 0.89, 0.91.
#   Sequelae: bronchiectasis, COPD, fatigue, reduced exercise capacity.
# Effect direction: lower post-TB utility reduces QALY gains from active TB
#   treatment → makes active-TB-detection strategies (Cough+CXR, symptom screen)
#   appear less cost-effective relative to base case.
# =============================================================================
cat("\n================================================================================\n")
cat("         POST-TB SEQUELAE — qaly_ActiveCompleted SENSITIVITY                   \n")
cat("================================================================================\n\n")

qaly_ac_base <- paramsData[["qaly_ActiveCompleted"]]

# Three values: full range from Ranzani et al. 2021 SR (0.88–0.91; post-TB sequelae)
qaly_ac_sa_values <- c(0.88, 0.89, 0.91)

cat(sprintf("Base: qaly_ActiveCompleted = %.2f (Kind 1998; full recovery assumed)\n", qaly_ac_base))
cat(sprintf("SA:   qaly_ActiveCompleted = %.2f / %.2f / %.2f (Ranzani 2021 SR range; post-TB sequelae)\n\n",
            qaly_ac_sa_values[1], qaly_ac_sa_values[2], qaly_ac_sa_values[3]))

sequelae_sa_results <- lapply(qaly_ac_sa_values, function(qac) {
  pd <- paramsData
  pd[["qaly_ActiveCompleted"]] <- qac
  res <- mclapply(strategies, run_strategy,
                  params   = pd,
                  mc.cores = max(1L, detectCores() - 1L))
  icer <- calculate_icer(res)
  icer$qaly_ac_value <- qac
  icer
})

frontier_strats_ac <- c("Passive case finding", "Cough+CXR (TB sx)",
                         "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)")

cat("Frontier ICERs — post-TB sequelae SA (qaly_ActiveCompleted range):\n")
cat(sprintf("  %-40s  %12s  %12s  %12s  %12s\n",
            "Strategy", "Base (0.92)", "Low (0.88)", "Mid (0.89)", "High (0.91)"))
cat(strrep("-", 90), "\n")
for (st in frontier_strats_ac) {
  icer_base_str <- {
    r <- icer_table[icer_table$strategy == st, ]
    if (nrow(r) == 0 || is.na(r$sequential_icer)) "(ref)" else
      paste0("\u00a3", format(round(r$sequential_icer), big.mark = ","))
  }
  sa_strs <- sapply(sequelae_sa_results, function(d) {
    r <- d[d$strategy == st, ]
    if (nrow(r) == 0 || is.na(r$sequential_icer)) "(ref)" else
      paste0("\u00a3", format(round(r$sequential_icer), big.mark = ","))
  })
  cat(sprintf("  %-40s  %12s  %12s  %12s  %12s\n",
              st, icer_base_str, sa_strs[1], sa_strs[2], sa_strs[3]))
}
cat("\n")

# Save CSV
sequelae_sa_df <- bind_rows(sequelae_sa_results) %>%
  mutate(scenario      = sprintf("qaly_ac_%.2f", qaly_ac_value),
         cost_per_person = cost / n_c,
         qaly_per_person = qaly / n_c)
write.csv(sequelae_sa_df, "output/csv/post_tb_sequelae_sensitivity.csv", row.names = FALSE)
cat("Saved: output/csv/post_tb_sequelae_sensitivity.csv\n")

# --- Supplementary figure: post-TB sequelae SA line plot ---
sequelae_plot_df <- sequelae_sa_df %>%
  filter(strategy %in% c("Cough+CXR (TB sx)", "Parallel Sx+QFT (Ultra)")) %>%
  filter(!is.na(sequential_icer)) %>%
  mutate(
    strategy    = factor(strategy, levels = c("Cough+CXR (TB sx)", "Parallel Sx+QFT (Ultra)")),
    qaly_ac_value = as.numeric(qaly_ac_value)
  )

cat("\n================================================================================\n")

# =============================================================================
# SEQUENTIAL LTBI DETECTION AT 17.8% PREVALENCE SENSITIVITY
#
# Background: sequential strategies use LTBI detection from Zenner 2025 Excel
#   (~100–150/100k), which embeds a cohort-level LTBI prevalence of ~0.1–0.15%.
#   The model background LTBI prevalence is 17.8% (Berrocal-Almanza 2022).
#   This SA asks: if sequential IGRA pathways were applied to the full 17.8%
#   LTBI pool, what would they detect? Formula: n_ltbi × igra_sensitivity
#   (QFT: 0.83; T-SPOT: 0.88). No specificity correction applied (no FP
#   deduction) → this is an upper bound on sequential IGRA LTBI yield.
#   Parallel strategies apply specificity and uptake in addition; this SA is
#   not directly comparable but indicates whether the asymmetry is material.
#   Strategies without IGRA (CXR-only, symptom screen) retain Excel LTBI = ~0.
#   Active TB TP/FN held fixed.
#
# Purpose: assess whether sequential ICERs change materially if sequential
#   LTBI yield is recalculated at the background 17.8% prevalence, i.e., whether
#   the sequential/parallel asymmetry affects the efficient frontier.
# =============================================================================
cat("\n================================================================================\n")
cat("   SEQUENTIAL LTBI DETECTION AT 17.8% PREVALENCE — STRUCTURAL SA              \n")
cat("================================================================================\n\n")

cat(sprintf("Base: sequential LTBI from Zenner 2025 Excel (~100-150/100k; ~0.15%% LTBI)\n"))
cat(sprintf("SA:   sequential LTBI = n_ltbi (%.0f) x IGRA sensitivity per pathway\n",
            n_ltbi_global))
cat(sprintf("      QFT-Plus: %.0f/100k (17,800 x %.2f)  T-SPOT: %.0f/100k (17,800 x %.2f)\n\n",
            n_ltbi_global * igra_sens_qft, igra_sens_qft,
            n_ltbi_global * igra_sens_tspt, igra_sens_tspt))

strategies_seq_ltbi <- lapply(strategies, function(s) {
  if (is.null(s$excel_name) || is.na(s$excel_name)) {
    # Passive case finding — no LTBI detected
    return(s)
  }
  if (!is.null(s$prog_cost_type)) {
    # Parallel IGRA — already uses population-based formula; retain as-is
    return(s)
  }
  # Sequential strategy: determine if IGRA is in the pathway
  nm_lo <- tolower(s$excel_name)
  if (grepl("qft|qfn|quantiferon", nm_lo)) {
    ltbi_sa <- n_ltbi_global * igra_sens_qft
  } else if (grepl("tspt|t-spot|tspot|elispot", nm_lo)) {
    ltbi_sa <- n_ltbi_global * igra_sens_tspt
  } else {
    # No IGRA: CXR, symptom screen, culture — retain Excel LTBI (near zero)
    return(s)
  }
  s$init["LatentDiagnosed"]   <- ltbi_sa
  s$init["LatentUndiagnosed"] <- n_ltbi_global - ltbi_sa
  s
})

strategy_results_seq_ltbi <- mclapply(strategies_seq_ltbi, run_strategy,
                                       params   = paramsData,
                                       mc.cores = max(1L, detectCores() - 1L))
icer_seq_ltbi <- calculate_icer(strategy_results_seq_ltbi)

cat("Frontier comparison — sequential LTBI detection base vs 17.8% population prevalence:\n")
cat(sprintf("  %-40s  %18s  %22s  %10s\n",
            "Strategy", "Base ICER (Excel)", "SA ICER (17.8% LTBI)", "Delta"))
cat(strrep("-", 96), "\n")

all_strats_check <- union(
  icer_table$strategy[icer_table$dominance %in% c("ref", "non-dominated")],
  icer_seq_ltbi$strategy[icer_seq_ltbi$dominance %in% c("ref", "non-dominated")]
)
for (st in all_strats_check) {
  r_base <- icer_table[icer_table$strategy == st, ]
  r_sa   <- icer_seq_ltbi[icer_seq_ltbi$strategy == st, ]
  icer_b <- if (nrow(r_base)==0||is.na(r_base$sequential_icer)) "Ref/NA" else
    paste0("\u00a3", format(round(r_base$sequential_icer), big.mark=","), "/QALY")
  icer_s <- if (nrow(r_sa)==0||is.na(r_sa$sequential_icer)) "Ref/NA" else
    paste0("\u00a3", format(round(r_sa$sequential_icer), big.mark=","), "/QALY")
  delta  <- if (nrow(r_base)==0||nrow(r_sa)==0||
                is.na(r_base$sequential_icer)||is.na(r_sa$sequential_icer)) "\u2014" else
    paste0(ifelse(r_sa$sequential_icer - r_base$sequential_icer >= 0, "+", ""),
           "\u00a3", format(round(r_sa$sequential_icer - r_base$sequential_icer), big.mark=","))
  cat(sprintf("  %-40s  %18s  %22s  %10s\n", st, icer_b, icer_s, delta))
}
cat("\n")

new_front <- icer_seq_ltbi$strategy[icer_seq_ltbi$dominance %in% c("ref", "non-dominated")]
cat(sprintf("SA efficient frontier (%d strategies): %s\n\n",
            length(new_front), paste(new_front, collapse=" \u2192 ")))

# Save CSV
seq_ltbi_sa_df <- bind_rows(
  icer_table %>% mutate(scenario = "base_excel_ltbi",
                        cost_per_person = cost / n_c,
                        qaly_per_person = qaly / n_c),
  icer_seq_ltbi %>% mutate(scenario = "sequential_ltbi_178pct",
                            cost_per_person = cost / n_c,
                            qaly_per_person = qaly / n_c)
)
write.csv(seq_ltbi_sa_df, "output/csv/sequential_ltbi_sa.csv", row.names = FALSE)
cat("Saved: output/csv/sequential_ltbi_sa.csv\n")

cat("\n================================================================================\n")

# =============================================================================
# TRANSMISSION ANALYSIS (SOCIETAL PERSPECTIVE)
# =============================================================================

if (!exists("strategy_results")) {
  source("scripts/model_functions.R")
  strategies       <- create_diagnostic_strategies()
  strategy_results <- mclapply(strategies, run_strategy,
                               mc.cores = max(1L, parallel::detectCores() - 1L))
}

# PRIMARY ANALYSIS: TRANSMISSION PREVENTION (SOCIETAL PERSPECTIVE)
#
# The NHS perspective (above) is the robustness check. This transmission prevention
# analysis is the primary result.
#
# Base case assumes no onward transmission from undiagnosed active TB cases
# (standard convention in UK TB CEA models; Jit 2011, Dale 2022, Green 2025).
# This analysis quantifies the value of averted secondary transmissions.
#
# Method (Brooks-Pollock et al. 2020, PLOS Comput Biol, PMID 32218567):
#   1. Estimate new active TB cases per strategy over the model horizon:
#      new_atb = p_reactivation × cumulative person-months in all reactivating
#      latent states (LatentUndiagnosed, LatentDiagnosed, LatentNotreated,
#      LatentLtfu at 0.000912/mo; LatentDiscontinued at 0.0008/mo).
#   2. Active TB cases prevented vs No Screening = no_screening_atb - strategy_atb.
#   3. Secondary active TB cases avoided = cases_prevented × beta_tx
#      (base: beta = 0.205, conservative base = 5 contacts/case/yr × 4.1%
#       transmission probability per contact; anchored to Brooks-Pollock et al.
#       2020 UK estimate (β=0.41 at 10 contacts);
#       DSA: 5 (base)/10 (mid)/15 (high) contacts/yr; PSA: contacts ~ Uniform(5,15)).
#   3b. Secondary LTBI cases avoided = cases_prevented × n_contacts × p_ltbi_per_contact
#       (Green et al. 2025: 23% LTBI prevalence in close contacts; SA: 18.3% UKHSA 2024).
#   4. Cost saving = active TB savings (Green et al. 2025: £6,055/case) +
#                    LTBI savings (Snoad et al. 2025: £458/LTBI case managed).
#   5. Subtract total savings from strategy costs → recalculate ICERs.
#
# This approach is consistent with Pareek et al. 2011 (Lancet ID) and
# Green et al. 2025 (ERJ Open Res) which apply a fixed secondary-case
# multiplier to active TB cases prevented.
#
# Outputs: icer_societal_perspective.csv
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("         PRIMARY ANALYSIS: TRANSMISSION PREVENTION (SOCIETAL PERSPECTIVE)         \n")
cat("================================================================================\n\n")

p_tx_per_contact      <- as.numeric(config_list["p_tx_per_contact"])
n_contacts_base       <- as.numeric(config_list["n_contacts_base"])
n_contacts_dsa_mid    <- as.numeric(config_list["n_contacts_dsa_mid"])
n_contacts_dsa_high   <- as.numeric(config_list["n_contacts_dsa_high"])
cost_per_secondary_tb   <- as.numeric(config_list["cost_secondary_tb"])
p_ltbi_per_contact      <- as.numeric(config_list["p_ltbi_per_contact"])
cost_per_secondary_ltbi <- as.numeric(config_list["cost_secondary_ltbi"])
beta_tx_base            <- n_contacts_base * p_tx_per_contact
p_react_standard      <- 0.000912  # monthly reactivation rate phase 1 (Berrocal-Almanza 2022, config row 41)
p_react_phase2_tx     <- as.numeric(config_list["p_react_phase2"])  # 0.000333, phase 2 rate
p_react_disc          <- 0.0008    # LatentDiscontinued reactivation (config row 74)
# Scale LatentDiscontinued rate for phase 2 proportionally (same phase switch logic as main model)
p_react_disc_phase2   <- p_react_disc * (p_react_phase2_tx / p_react_standard)

cat(sprintf("Beta (secondary cases per active TB case, UK): %.2f (Brooks-Pollock 2020, 95%% CrI 0.30-0.60)\n", beta_tx_base))
cat(sprintf("Contacts parameterisation: %d contacts/case/yr × %.3f transmission probability = beta %.2f\n",
            n_contacts_base, p_tx_per_contact, beta_tx_base))
cat(sprintf("Cost per secondary active TB case: £%s (Green et al. 2025 ERJ Open Res, PMC12183743; NHS direct costs, 2024 GBP)\n",
            format(cost_per_secondary_tb, big.mark = ",")))
cat(sprintf("Secondary LTBI per contact: %.2f (Green et al. 2025); cost per LTBI case managed: £%s (Snoad et al. 2025 medRxiv)\n",
            p_ltbi_per_contact, format(cost_per_secondary_ltbi, big.mark = ",")))
cat(sprintf("(DSA: contacts = %g / %g / %g per untreated case/yr; PSA: Uniform(%g, %g))\n\n",
            n_contacts_base, n_contacts_dsa_mid, n_contacts_dsa_high,
            n_contacts_base, n_contacts_dsa_high))

# Step 1: Estimate new active TB cases per strategy over the 55-year horizon.
# New active TB cases arise from reactivation out of all reactivating latent states.
# Approximation: new_cases_from_state = p_reactivation × cumulative person-months in state
# (standard Markov cohort incidence approximation; cycle length = 1 month).
# Two-phase rates applied: phase 1 (cycles 1–60) = 0.000912/mo;
#                          phase 2 (cycles 61–660) = 0.000333/mo (Horsburgh 2004 NEJM).
latent_states_standard <- c("LatentUndiagnosed", "LatentDiagnosed",
                             "LatentNotreated",  "LatentLtfu")

new_atb_cases <- sapply(strategy_results, function(x) {
  sm      <- x$state_membership
  n_cy    <- nrow(sm)
  cut     <- min(60L, n_cy)         # phase 1: cycles 1-60
  # Phase 1 person-months
  std_p1  <- sum(sm[seq_len(cut),          latent_states_standard])
  disc_p1 <- sum(sm[seq_len(cut),          "LatentDiscontinued"])
  # Phase 2 person-months (cycles 61+)
  std_p2  <- if (n_cy > cut) sum(sm[(cut + 1L):n_cy, latent_states_standard]) else 0
  disc_p2 <- if (n_cy > cut) sum(sm[(cut + 1L):n_cy, "LatentDiscontinued"])   else 0

  std_p1  * p_react_standard   + disc_p1 * p_react_disc +
  std_p2  * p_react_phase2_tx  + disc_p2 * p_react_disc_phase2
})
names(new_atb_cases) <- sapply(strategy_results, function(x) x$strategy_name)

# Step 2: Active TB cases prevented vs No Screening
atb_no_screening  <- new_atb_cases["Passive case finding"]
atb_cases_prevent <- atb_no_screening - new_atb_cases   # positive = screening prevents

# Step 3: Secondary active TB cases avoided = active TB cases prevented × beta
secondary_cases_avoided <- atb_cases_prevent * beta_tx_base

# Step 3b: Secondary LTBI cases avoided = active TB cases prevented × n_contacts × p_ltbi_per_contact
secondary_ltbi_avoided <- atb_cases_prevent * n_contacts_base * p_ltbi_per_contact

# Step 4: Total cost savings = active TB savings + LTBI contact management savings
transmission_savings <- secondary_cases_avoided * cost_per_secondary_tb +
                        secondary_ltbi_avoided   * cost_per_secondary_ltbi

cat("New active TB cases (from reactivation) and transmission savings by strategy:\n")
cat(sprintf("  %-28s %12s %14s %14s %14s %16s\n",
            "Strategy", "New ATB", "ATB prev.", "Sec.ATB avoid.", "Sec.LTBI avoid.", "Savings (£, 100k)"))
cat(strrep("-", 104), "\n")
for (nm in names(new_atb_cases)) {
  cat(sprintf("  %-28s %12.0f %14.0f %14.1f %15.1f %16s\n",
              nm,
              new_atb_cases[nm],
              atb_cases_prevent[nm],
              secondary_cases_avoided[nm],
              secondary_ltbi_avoided[nm],
              format(round(transmission_savings[nm]), big.mark = ",")))
}
cat("\n")

# Save transmission prevention counts to CSV
tx_counts_df <- data.frame(
  strategy                  = names(new_atb_cases),
  new_atb_cases             = as.numeric(new_atb_cases),
  atb_cases_prevented       = as.numeric(atb_cases_prevent),
  secondary_cases_avoided   = as.numeric(secondary_cases_avoided),
  secondary_ltbi_avoided    = as.numeric(secondary_ltbi_avoided),
  transmission_saving_total = as.numeric(transmission_savings),
  transmission_saving_pp    = as.numeric(transmission_savings) / n_c,
  stringsAsFactors = FALSE
)
write.csv(tx_counts_df, "output/csv/transmission_prevention_counts.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_prevention_counts.csv\n")

# Step 5: Build modified strategy_results with adjusted total_cost, then recalculate ICERs
strategy_results_transmission <- lapply(strategy_results, function(x) {
  x_mod <- x
  savings <- transmission_savings[x$strategy_name]
  if (!is.na(savings)) {
    x_mod$total_cost <- x$total_cost - savings
  }
  x_mod
})

icer_societal_perspective <- calculate_icer(strategy_results_transmission)

cat(sprintf("ICER TABLE — PRIMARY ANALYSIS: Transmission Prevention / Societal Perspective (beta = %.2f, Brooks-Pollock 2020):\n", beta_tx_base))
cat(sprintf("%-22s %11s %10s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 108), "\n")
for (i in 1:nrow(icer_societal_perspective)) {
  icer_str     <- if (is.na(icer_societal_perspective$icer[i])) "Ref"
                  else paste0("\u00a3", format(round(icer_societal_perspective$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_societal_perspective$sequential_icer[i])) "-"
                  else paste0("\u00a3", format(round(icer_societal_perspective$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %14s %15s  %s\n",
              icer_societal_perspective$strategy[i],
              paste0("\u00a3", format(round(icer_societal_perspective$cost_per_person[i]), big.mark = ",")),
              icer_societal_perspective$qaly_per_person[i],
              paste0("\u00a3", format(round(icer_societal_perspective$inc_cost[i] / n_c), big.mark = ",")),
              icer_str,
              seq_icer_str,
              icer_societal_perspective$dominance[i]))
}
cat("\n================================================================================\n\n")

write.csv(icer_societal_perspective %>%
  mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
  "output/csv/icer_societal_perspective.csv", row.names = FALSE)
cat("Saved: output/csv/icer_societal_perspective.csv\n")


# CE scatter delegated to run_plots_and_tables.R (reads icer_societal_perspective.csv).
source("scripts/run_plots_and_tables.R")

# Frontier-only arrow plot: NHS perspective → societal perspective (transmission prevention)
# Shows frontier strategies with horizontal arrows (costs shift left = cheaper).
# QALYs unchanged; only costs change with transmission savings.

# Union of base-case + transmission frontier strategies:
# Frontier composition is unchanged between NHS and societal perspectives
# (same 4 strategies: Passive case finding, Cough+CXR, Symptom screen+CXR, Parallel Sx+QFT).
# Union taken dynamically; no strategies added or dropped.
tx_frontier_names <- union(
  frontier_names,
  icer_societal_perspective %>%
    filter(dominance %in% c("ref", "non-dominated")) %>%
    arrange(cost) %>% pull(strategy)
)

frontier_base_tx <- icer_table %>%
  filter(strategy %in% tx_frontier_names) %>%
  select(strategy, q_base = qaly, c_base = cost) %>%
  mutate(q_base = q_base / n_c, c_base = c_base / n_c)

frontier_scen_tx <- icer_societal_perspective %>%
  filter(strategy %in% tx_frontier_names) %>%
  select(strategy, q_scen = qaly, c_scen = cost) %>%
  mutate(q_scen = q_scen / n_c, c_scen = c_scen / n_c)

frontier_arrows_tx <- left_join(frontier_base_tx, frontier_scen_tx, by = "strategy") %>%
  mutate(strategy = factor(strategy, levels = tx_frontier_names))

lbl_base_tx <- "NHS perspective robustness check (no secondary cases)"
lbl_scen_tx <- sprintf("Primary analysis: societal perspective (\u03b2 = %.2f; \u00a3%s/secondary case)",
                        beta_tx_base,
                        format(cost_per_secondary_tb, big.mark = ","))

frontier_long_tx <- bind_rows(
  frontier_arrows_tx %>% transmute(strategy, q = q_base, c = c_base,
                                   cond = lbl_base_tx),
  frontier_arrows_tx %>% transmute(strategy, q = q_scen, c = c_scen,
                                   cond = lbl_scen_tx)
)

label_df_tx <- frontier_arrows_tx %>%
  transmute(strategy, q = q_base, c = c_base,
            is_igra = grepl("^Parallel", strategy, ignore.case = FALSE))

tx_colours <- c(LANCET[4], LANCET[4])
names(tx_colours) <- c(lbl_base_tx, lbl_scen_tx)
tx_shapes <- c(16L, 17L)
names(tx_shapes) <- c(lbl_base_tx, lbl_scen_tx)

p_transmission_compare <- ggplot() +
  geom_path(data = frontier_arrows_tx %>% arrange(q_base),
            aes(x = q_base, y = c_base),
            colour = LANCET[4], linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_path(data = frontier_arrows_tx %>% arrange(q_scen),
            aes(x = q_scen, y = c_scen),
            colour = LANCET[4], linetype = "solid", linewidth = 0.6, alpha = 0.7) +
  geom_segment(
    data = frontier_arrows_tx,
    aes(x = q_base, y = c_base, xend = q_scen, yend = c_scen),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    colour = "grey40", linewidth = 0.5
  ) +
  geom_point(data = frontier_long_tx,
             aes(x = q, y = c, colour = cond, shape = cond),
             size = 4, alpha = 0.95) +
  # IGRA strategy labels in bold; non-IGRA in grey
  ggrepel::geom_text_repel(
    data = label_df_tx %>% filter(!is_igra),
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = "grey30",
    nudge_y = -8, box.padding = 0.3,
    segment.color = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_df_tx %>% filter(is_igra),
    aes(x = q, y = c, label = strategy),
    size = 3.2, colour = LANCET[7], fontface = "bold",
    nudge_x = 0.004, nudge_y = c(18, 6),
    box.padding = 0.4, point.padding = 0.3,
    segment.color = LANCET[7], max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = tx_colours, name = NULL) +
  scale_shape_manual(values  = tx_shapes,  name = NULL) +
  theme_minimal(base_size = 13) +
  labs(
    x        = "QALYs per person",
    y        = "Cost per person (\u00a3)",
    title    = NULL,
    subtitle = sprintf(
      "Arrows show cost reduction when secondary cases are included | Frontier strategies only\n\u03b2 = %.2f secondary cases per active TB case | Cost per secondary case: \u00a3%s\n\u25b2 Parallel IGRA strategies (bold) benefit most — more LTBI detected \u2192 more secondary cases prevented",
      beta_tx_base, format(cost_per_secondary_tb, big.mark = ","))
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(color = "grey40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10)
  )

# scenario_transmission_compare removed (session 55): content covered by
# nhs_vs_societal_dumbbell_supplementary and scenario_frontier_comparison_supplementary.


# nhs_vs_societal_dumbbell removed: two of three strategies show NHS ≈ societal (flat),
# making the figure uninformative. The £14,344→£12,738 shift for Parallel Sx+QFT is
# fully described in Para 2 of the narrative and in the ICER table.
# File deleted; removed from verify_project.R required list.

cat("\n================================================================================\n")

# =============================================================================
# TRANSMISSION DSA
# =============================================================================
# TRANSMISSION DSA: contacts/case/yr (5/10/15) × cost (£6,055/£12,110)
#
# Contacts DSA range: 5 (conservative) / 10 (base) / 15 (high) per untreated
# active TB case per year; consistent with WHO/CDC contact investigation literature.
# Implied beta: contacts × p_tx_per_contact (0.041), anchored so that
# 10 contacts × 0.041 = beta 0.41 (Brooks-Pollock 2020 UK estimate).
# cost DSA: £6,055 base vs £12,110 upper (Green 2025 ERJ Open Res DSA range).
# Active TB cases per strategy held fixed from base case run above.
# =============================================================================

tx_dsa_contacts <- c(n_contacts_base, n_contacts_dsa_mid, n_contacts_dsa_high)
tx_dsa_results <- lapply(tx_dsa_contacts, function(nc) {
  beta_eff <- nc * p_tx_per_contact
  savings  <- atb_cases_prevent * beta_eff * cost_per_secondary_tb
  sr_mod   <- lapply(strategy_results, function(x) {
    x$total_cost <- x$total_cost - savings[x$strategy_name]
    x
  })
  tbl <- calculate_icer(sr_mod)
  tbl$n_contacts <- nc
  tbl$beta_eff   <- beta_eff
  tbl
})
tx_dsa_df <- bind_rows(tx_dsa_results)

# Save DSA transmission table
write.csv(tx_dsa_df %>% mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c),
          "output/csv/transmission_dsa_results.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_dsa_results.csv\n")

# Print sequential ICERs for frontier strategies across contact values
cat("\nTransmission contacts DSA — sequential ICERs for Parallel Sx+QFT (Ultra):\n")
cat(sprintf("  %-30s %12s %22s %15s\n",
    "Contacts/case/yr", "Implied beta", "Parallel Sx+QFT cost/pp", "Seq. ICER"))
for (nc in tx_dsa_contacts) {
  row <- tx_dsa_df[tx_dsa_df$n_contacts == nc &
                   tx_dsa_df$strategy == "Parallel Sx+QFT (Ultra)", ]
  cat(sprintf("  %-30s %12s %22s %15s\n",
    sprintf("%d contacts%s", nc, if (nc == n_contacts_base) " (base)" else ""),
    sprintf("\u03b2 = %.3f", nc * p_tx_per_contact),
    paste0("\u00a3", format(round(row$cost / n_c), big.mark = ",")),
    if (length(row$sequential_icer) == 0 || is.na(row$sequential_icer))
      "\u2014" else paste0("\u00a3", format(round(row$sequential_icer), big.mark = ","))))
}
cat("\n")

# Grouped bar chart: Cough+CXR + Parallel Cough+QFT + Parallel Sx+QFT × 3 contact values
tx_dsa_strats <- c("Cough+CXR (TB sx)", "Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)")
tx_dsa_cols <- c(
  "Cough+CXR (TB sx)"          = LANCET[4],   # teal
  "Parallel Cough+QFT (Ultra)" = "#7a8a99",   # steel gray (ext. dominated)
  "Parallel Sx+QFT (Ultra)"    = LANCET[7]    # dark red
)
tx_dsa_short <- c(
  "Cough+CXR (TB sx)"          = "Cough+CXR (TB sx)",
  "Parallel Cough+QFT (Ultra)" = "Parallel Cough+QFT (Ultra)",
  "Parallel Sx+QFT (Ultra)"    = "Parallel Sx+QFT (Ultra)"
)
tx_dsa_xlabs <- c(
  "5"  = "5 contacts\n(base)",
  "10" = "10 contacts\n(mid)",
  "15" = "15 contacts\n(high)"
)

# transmission_dsa_bar_supplementary removed: data available in output/csv/transmission_dsa_results.csv
# and ST4g sheet in Sensitivity_Analyses_supplementary.xlsx.

# Cost sensitivity: base (£6,055) vs upper (£12,110) — Green 2025 ERJ Open Res DSA range
cat("\nCost per secondary TB case sensitivity (5 contacts base; Parallel Sx+QFT (Ultra)):\n")
cat(sprintf("  %-18s %30s %15s\n", "Cost/case", "Parallel Sx+QFT cost/pp", "Seq. ICER"))
for (cost_test in c(cost_per_secondary_tb, as.numeric(config_list["cost_secondary_tb_upper"]))) {
  sav_test <- atb_cases_prevent * beta_tx_base * cost_test
  sr_cost <- lapply(strategy_results, function(x) {
    x$total_cost <- x$total_cost - sav_test[x$strategy_name]; x })
  tbl_cost <- calculate_icer(sr_cost)
  row_u    <- tbl_cost[tbl_cost$strategy == "Parallel Sx+QFT (Ultra)", ]
  cat(sprintf("  \u00a3%-17s %30s %15s\n",
    format(cost_test, big.mark = ","),
    paste0("\u00a3", format(round(row_u$cost / n_c), big.mark = ",")),
    if (length(row_u$sequential_icer) == 0 || is.na(row_u$sequential_icer)) "\u2014"
    else paste0("\u00a3", format(round(row_u$sequential_icer), big.mark = ","))))
}
cat("\n")


# =============================================================================
# PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
# =============================================================================

if (!exists("strategies")) {
  source("scripts/model_functions.R")
  strategies <- create_diagnostic_strategies()
}

if (!exists("atb_cases_prevent")) {
  tx <- read.csv("output/csv/transmission_prevention_counts.csv",
                 stringsAsFactors = FALSE)
  atb_cases_prevent <- setNames(tx$atb_cases_prevented, tx$strategy)
}

# Transmission parameters — load from config if not inherited from run_transmission.R
if (!exists("n_contacts_base")) {
  p_tx_per_contact      <- as.numeric(config_list["p_tx_per_contact"])
  n_contacts_base       <- as.numeric(config_list["n_contacts_base"])
  n_contacts_dsa_high   <- as.numeric(config_list["n_contacts_dsa_high"])
  cost_per_secondary_tb   <- as.numeric(config_list["cost_secondary_tb"])
  p_ltbi_per_contact      <- as.numeric(config_list["p_ltbi_per_contact"])
  cost_per_secondary_ltbi <- as.numeric(config_list["cost_secondary_ltbi"])
}

# PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
#
# 1,000 MC simulations; params drawn from config.csv distributions; fixed seed.
# Outputs: CE plane, CEAC, ICER 95% CrIs.
# =============================================================================

set.seed(42)  # fixed seed ensures results are reproducible across runs
psa_df <- run_psa(strategies_base = strategies, n_sim = 1000)

# Save PSA raw results for supplementary materials and audit trail
write.csv(psa_df, "output/csv/psa_results.csv", row.names = FALSE)

# =============================================================================
# TRANSMISSION PSA: PROBABILISTIC SENSITIVITY ANALYSIS WITH CONTACTS
#
# n_contacts ~ Uniform(5, 15) per untreated active TB case per year.
# Implied beta_i = n_contacts_i × p_tx_per_contact (0.041), covering the
# full range from 5 (conservative) to 15 (high) contacts per case.
# This corresponds to beta ~ Uniform(0.21, 0.62), consistent with and
# slightly wider than Brooks-Pollock 2020 95% CrI (0.30-0.60).
# Active TB cases prevented held fixed from base case run; only contacts
# vary per simulation. Savings include both active TB and LTBI secondary
# components (consistent with deterministic societal analysis). Merged
# with main PSA cost draws.
# =============================================================================
cat("\nRunning transmission PSA (1,000 simulations, contacts ~ Uniform(5, 15))...\n")

# No separate seed — transmission PSA continues from RNG state after main PSA.
# Both analyses are fully reproducible from the single seed set before run_psa().
psa_tx_results <- list()
for (i in 1:1000) {
  if (i %% 200 == 0) cat(sprintf("  Transmission PSA simulation %d/1000\n", i))
  n_contacts_i <- runif(1, min = n_contacts_base, max = n_contacts_dsa_high)
  beta_i       <- n_contacts_i * p_tx_per_contact
  savings_i    <- atb_cases_prevent * beta_i * cost_per_secondary_tb +
                  atb_cases_prevent * n_contacts_i * p_ltbi_per_contact * cost_per_secondary_ltbi
  psa_tx_results[[i]] <- tibble(
    sim        = i,
    strategy   = names(atb_cases_prevent),
    n_contacts = n_contacts_i,
    beta_tx    = beta_i,
    savings    = as.numeric(savings_i)
  )
}
psa_tx_savings <- bind_rows(psa_tx_results)

psa_tx_df <- psa_df %>%
  left_join(psa_tx_savings, by = c("sim", "strategy")) %>%
  mutate(cost_tx = cost - savings)

ref_tx <- psa_tx_df %>% filter(strategy == "Passive case finding") %>%
  select(sim, ref_cost_tx = cost_tx, ref_qaly = qaly)

psa_tx_icer <- psa_tx_df %>%
  filter(strategy != "Passive case finding") %>%
  left_join(ref_tx, by = "sim") %>%
  mutate(
    inc_cost_tx = cost_tx - ref_cost_tx,
    inc_qaly    = qaly - ref_qaly,
    icer_tx     = ifelse(inc_qaly > 0, inc_cost_tx / inc_qaly, NA_real_),
    ce_25k      = icer_tx < 25000,
    ce_35k      = icer_tx < 35000
  )

frontier_tx_summary <- psa_tx_icer %>%
  filter(strategy %in% c("Cough+CXR (TB sx)",
                          "Parallel Cough+QFT (Ultra)",
                          "Parallel Sx+QFT (Ultra)")) %>%
  group_by(strategy) %>%
  summarise(
    median_icer = median(icer_tx, na.rm = TRUE),
    lo_icer     = quantile(icer_tx, 0.025, na.rm = TRUE),
    hi_icer     = quantile(icer_tx, 0.975, na.rm = TRUE),
    prob_ce_25k = mean(ce_25k, na.rm = TRUE),
    prob_ce_35k = mean(ce_35k, na.rm = TRUE),
    .groups     = "drop"
  )

cat("\nTransmission PSA — frontier strategies (contacts ~ Uniform(5, 15); implied beta ~ Uniform(0.21, 0.62)):\n")
cat(sprintf("  %-22s %12s %24s %12s %12s\n",
    "Strategy", "Median ICER", "95% CrI", "P(CE \u00a325k)", "P(CE \u00a335k)"))
cat(strrep("-", 88), "\n")
for (i in seq_len(nrow(frontier_tx_summary))) {
  r <- frontier_tx_summary[i, ]
  cat(sprintf("  %-22s %12s %24s %11.0f%% %11.0f%%\n",
    r$strategy,
    paste0("\u00a3", format(round(r$median_icer), big.mark = ",")),
    paste0("\u00a3", format(round(r$lo_icer), big.mark = ","), " to \u00a3",
           format(round(r$hi_icer), big.mark = ",")),
    r$prob_ce_25k * 100,
    r$prob_ce_35k * 100))
}
write.csv(frontier_tx_summary, "output/csv/transmission_psa_summary.csv", row.names = FALSE)
cat("Saved: output/csv/transmission_psa_summary.csv\n")
cat("\n================================================================================\n")

# Summarise mean and standard deviation of cost and QALYs across simulations
psa_summary <- psa_df %>%
  group_by(strategy) %>%
  summarise(
    mean_cost = mean(cost),
    sd_cost = sd(cost),
    mean_qaly = mean(qaly),
    sd_qaly = sd(qaly),
    .groups = "drop"
  )

cat("\n")
cat("================================================================================\n")
cat("                    PSA RESULTS (1,000 simulations)                             \n")
cat("================================================================================\n\n")
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "Strategy", "Mean Cost", "SD Cost", "Mean QALYs", "SD QALYs"))
cat(sprintf("%-20s %14s %14s %14s %14s\n",
            "--------------------", "--------------", "--------------", "--------------", "--------------"))
for (i in 1:nrow(psa_summary)) {
  cat(sprintf("%-20s %13s %13s %14.1f %14.1f\n",
              psa_summary$strategy[i],
              paste0("£", format(round(psa_summary$mean_cost[i]), big.mark = ",")),
              paste0("£", format(round(psa_summary$sd_cost[i]), big.mark = ",")),
              psa_summary$mean_qaly[i],
              psa_summary$sd_qaly[i]))
}
cat("\n================================================================================\n")

# -------------------- ICER credible intervals from PSA ------------------------
# ICER 95% CrI from 2.5th/97.5th percentiles; P(CE) = proportion of sims with ICER < £25k.
ref_psa <- psa_df %>% filter(strategy == "Passive case finding") %>%
  select(sim, ref_cost = cost, ref_qaly = qaly)

icer_psa <- psa_df %>%
  filter(strategy != "Passive case finding") %>%
  left_join(ref_psa, by = "sim") %>%
  mutate(
    inc_cost = cost - ref_cost,
    inc_qaly = qaly - ref_qaly,
    icer = ifelse(inc_qaly != 0, inc_cost / inc_qaly, NA)
  )

icer_ci <- icer_psa %>%
  group_by(strategy) %>%
  summarise(
    mean_icer = mean(icer, na.rm = TRUE),
    median_icer = median(icer, na.rm = TRUE),
    icer_lo = quantile(icer, 0.025, na.rm = TRUE),
    icer_hi = quantile(icer, 0.975, na.rm = TRUE),
    mean_inc_cost = mean(inc_cost),
    mean_inc_qaly = mean(inc_qaly),
    prob_ce_25k = mean(inc_qaly > 0 & (inc_cost / inc_qaly) < 25000, na.rm = TRUE),
    prob_ce_35k = mean(inc_qaly > 0 & (inc_cost / inc_qaly) < 35000, na.rm = TRUE),
    prob_dominated = mean(inc_cost > 0 & inc_qaly <= 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_inc_cost)

cat("\n")
cat("================================================================================\n")
cat("                    ICER WITH 95% CREDIBLE INTERVALS (PSA)                     \n")
cat("================================================================================\n\n")
cat(sprintf("%-20s %12s %12s %24s %10s %10s %8s\n",
            "Strategy", "Mean ICER", "Median ICER", "95% CrI", "P(CE £25k)", "P(CE £35k)", "P(Dom)"))
cat(sprintf("%-20s %12s %12s %24s %10s %10s %8s\n",
            "--------------------", "------------", "------------",
            "------------------------", "----------", "----------", "--------"))
for (i in 1:nrow(icer_ci)) {
  mean_str <- paste0("£", format(round(icer_ci$mean_icer[i]), big.mark = ","))
  med_str  <- paste0("£", format(round(icer_ci$median_icer[i]), big.mark = ","))
  lo_str   <- format(round(icer_ci$icer_lo[i]), big.mark = ",")
  hi_str   <- format(round(icer_ci$icer_hi[i]), big.mark = ",")
  cri_str  <- paste0("(£", lo_str, " to £", hi_str, ")")
  cat(sprintf("%-20s %12s %12s %24s %7.1f%% %7.1f%% %7.1f%%\n",
              icer_ci$strategy[i], mean_str, med_str, cri_str,
              100 * icer_ci$prob_ce_25k[i],
              100 * icer_ci$prob_ce_35k[i],
              100 * icer_ci$prob_dominated[i]))
}
cat("\nP(CE) = Probability cost-effective at £25,000/QALY WTP threshold")
cat("\nP(Dom) = Probability dominated (more costly and fewer QALYs vs No Screening)\n")
cat("\n================================================================================\n")

