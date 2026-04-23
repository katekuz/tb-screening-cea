# run_plots_and_tables.R — regenerate all figures and supplementary tables from output CSVs
# Requires: output/csv/ populated by MasterTBModel.R + run_dsa.R + run_transmission.R + run_psa.R
# Usage:    Rscript scripts/run_plots_and_tables.R
#           Or sourced by MasterTBModel.R / run_dsa.R / run_transmission.R

if (!file.exists("input/config.csv"))
  stop("Run from project root: Rscript scripts/run_plots_and_tables.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(ggsci)
  library(patchwork)
  library(cowplot)
  library(tidyr)
  library(forcats)
})

# =============================================================================
# SHARED PALETTE — single source of truth for all figures
# pal_lancet("lanonc", alpha = 0.85): index 1=navy, 2=red, 3=green, 4=teal, 5=purple, 8=grey
# =============================================================================
LANCET <- pal_lancet("lanonc", alpha = 0.85)(8)
NAVY   <- LANCET[1]

STRAT_COLS <- c(
  "Passive case finding"       = LANCET[8],   # gray
  "Cough+CXR (TB sx)"          = LANCET[4],   # teal
  "Symptom screen+CXR"         = LANCET[1],   # navy
  "Parallel Sx+QFT (Ultra)"    = LANCET[7],   # dark red
  "Parallel Cough+QFT (Ultra)" = "#7a8a99"    # steel gray (SA only, ext. dominated)
)
FRONTIER <- names(STRAT_COLS)[1:4]
DOM_GREY <- "#c0cad2"

# =============================================================================
# FIGURE 1 — CE scatter (body_ce_scatter.tiff)
# All 43 strategies; societal perspective; efficient frontier highlighted
# =============================================================================
cat("Generating CE scatter...\n")

icer <- read.csv("output/csv/icer_societal_perspective.csv", check.names = FALSE)

icer <- icer %>%
  mutate(
    layer      = case_when(
      strategy %in% FRONTIER ~ "frontier",
      TRUE                   ~ "other"
    ),
    label_text = case_when(
      strategy == "Passive case finding" ~ "Passive case finding\n(comparator)",
      strategy == "Cough+CXR (TB sx)"   ~ "Cough+CXR\u2020",
      TRUE                               ~ strategy
    ),
    fill_group = case_when(layer == "frontier" ~ strategy, TRUE ~ "Non-frontier"),
    xval = qaly_per_person,
    yval = cost_per_person
  )

ref   <- icer %>% filter(strategy == "Passive case finding")
x_ref <- ref$xval;  y_ref <- ref$yval
frontier_path <- icer %>% filter(layer == "frontier") %>% arrange(xval)

xr   <- diff(range(icer$xval));  yr <- diff(range(icer$yval))
f_mid_x <- mean(frontier_path$xval[3:4]);  f_mid_y <- mean(frontier_path$yval[3:4])
psq <- icer %>% filter(strategy == "Parallel Sx+QFT (Ultra)")
psq_x <- psq$xval;  psq_y <- psq$yval

p_scatter <- ggplot() +
  geom_path(data = frontier_path,
            aes(x = xval, y = yval, linetype = "Efficient frontier"),
            colour = "grey32", linewidth = 0.9) +
  geom_point(data = icer %>% filter(layer == "other"),
             aes(x = xval, y = yval, fill = fill_group),
             shape = 21, colour = "white", size = 5, stroke = 0.5) +
  geom_point(data = icer %>% filter(layer == "frontier"),
             aes(x = xval, y = yval, fill = fill_group),
             shape = 21, colour = "white", size = 10, stroke = 0.7) +
  ggrepel::geom_label_repel(
    data = icer %>% filter(strategy == "Cough+CXR (TB sx)"),
    aes(x = xval, y = yval, label = label_text, colour = strategy),
    size = 5.5, fontface = "bold", show.legend = FALSE,
    nudge_y = -yr * 0.18, direction = "x",
    label.size = 0.25, label.r = unit(0.15, "lines"),
    label.padding = unit(0.3, "lines"),
    fill = alpha("white", 0.85),
    segment.colour = "grey35", segment.size = 0.5, segment.alpha = 0.9,
    min.segment.length = 0.1, seed = 101) +
  ggrepel::geom_label_repel(
    data = icer %>% filter(strategy == "Symptom screen+CXR"),
    aes(x = xval, y = yval, label = label_text, colour = strategy),
    size = 5.5, fontface = "bold", show.legend = FALSE,
    nudge_y = yr * 0.18, direction = "x",
    label.size = 0.25, label.r = unit(0.15, "lines"),
    label.padding = unit(0.3, "lines"),
    fill = alpha("white", 0.85),
    segment.colour = "grey35", segment.size = 0.5, segment.alpha = 0.9,
    min.segment.length = 0.1, seed = 102) +
  annotate("segment",
           x = psq_x, xend = psq_x - xr * 0.10, y = psq_y, yend = psq_y,
           colour = "grey35", linewidth = 0.5, alpha = 0.9) +
  annotate("label", x = psq_x - xr * 0.10, y = psq_y,
           label = "Parallel Sx+QFT (Ultra)",
           colour = STRAT_COLS["Parallel Sx+QFT (Ultra)"],
           fill = alpha("white", 0.85), size = 5.5, fontface = "bold",
           label.r = unit(0.15, "lines"), label.padding = unit(0.3, "lines"),
           hjust = 1, vjust = 0.5) +
  annotate("segment",
           x = x_ref, xend = x_ref, y = y_ref, yend = y_ref + yr * 0.13,
           colour = "grey35", linewidth = 0.5, alpha = 0.9) +
  annotate("label", x = x_ref, y = y_ref + yr * 0.14,
           label = "Passive case finding\n(comparator)",
           colour = STRAT_COLS["Passive case finding"],
           fill = alpha("white", 0.85), size = 5.5, fontface = "bold",
           label.r = unit(0.15, "lines"), label.padding = unit(0.3, "lines"),
           hjust = 0.5, vjust = 0) +
  scale_fill_manual(
    values = c(STRAT_COLS, "Non-frontier" = DOM_GREY),
    limits = c(FRONTIER, "Non-frontier"),
    breaks = c(FRONTIER, "Non-frontier"),
    labels = c("Passive case finding (comparator)",
               "Cough+CXR (TB sx)",
               "Symptom screen+CXR",
               "Parallel Sx+QFT (Ultra)",
               "All other strategies"),
    name = NULL) +
  scale_colour_manual(values = c(STRAT_COLS, "Non-frontier" = "white"), guide = "none") +
  scale_linetype_manual(values = c("Efficient frontier" = "solid"), name = NULL) +
  scale_x_continuous(labels = function(x) formatC(x, format = "f", digits = 3)) +
  scale_y_continuous(
    labels = function(y)
      paste0("\u00a3", formatC(as.integer(round(y)), format = "d", big.mark = ","))) +
  coord_cartesian(
    xlim = c(min(icer$xval) - xr * 0.12, max(icer$xval) + xr * 0.06),
    ylim = c(min(icer$yval) - yr * 0.08, max(icer$yval) + yr * 0.12)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey93", linewidth = 0.35),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.text       = element_text(size = 10),
    legend.key.size   = unit(0.55, "cm"),
    legend.spacing.x  = unit(0.3, "cm"),
    legend.margin     = margin(t = 6, r = 4, b = 4, l = 4),
    axis.title        = element_text(size = 11, colour = "grey20"),
    axis.text         = element_text(size = 10, colour = "grey30")) +
  labs(x = "Lifetime QALYs per person", y = "Lifetime cost per person (\u00a3)") +
  guides(
    fill = guide_legend(order = 1, nrow = 2,
      override.aes = list(
        shape  = rep(21L, 5),
        fill   = c(unname(STRAT_COLS[1:4]), DOM_GREY),
        colour = rep("white", 5),
        size   = c(5, 5, 5, 5, 4), stroke = rep(0.7, 5))),
    linetype = guide_legend(order = 2,
      override.aes = list(colour = "grey32", linewidth = 1.0)))

ggsave("output/body/ce_scatter.tiff", p_scatter, width = 16, height = 12, dpi = 300,
       compression = "lzw")
cat("Saved: output/body/ce_scatter.tiff\n")

# =============================================================================
# FIGURE 2 — CEAC panel (body_ceac_panel.tiff)
# Societal (solid) and NHS (dashed) overlay; frontier strategies only
# =============================================================================
cat("Generating CEAC panel...\n")

set.seed(20260318)

cfg <- read.csv("input/config.csv", stringsAsFactors = FALSE)
config_list <- setNames(cfg$value, cfg$parameter)

p_tx_per_contact        <- as.numeric(config_list["p_tx_per_contact"])
cost_per_secondary_tb   <- as.numeric(config_list["cost_secondary_tb"])
p_ltbi_per_contact      <- as.numeric(config_list["p_ltbi_per_contact"])
cost_per_secondary_ltbi <- as.numeric(config_list["cost_secondary_ltbi"])
n_contacts_psa_min      <- as.numeric(config_list["n_contacts_base"])
n_contacts_psa_max      <- as.numeric(config_list["n_contacts_dsa_high"])

n_sim     <- 1000
wtp_range <- seq(0, 50000, by = 500)

strat_cols_ceac <- STRAT_COLS[names(STRAT_COLS) %in% FRONTIER]

compute_ceac <- function(psa_filt, wtp_range) {
  bind_rows(lapply(wtp_range, function(w) {
    psa_filt %>%
      mutate(nmb = qaly * w - cost) %>%
      group_by(sim) %>%
      mutate(is_best = nmb == max(nmb)) %>%
      ungroup() %>%
      group_by(strategy) %>%
      summarise(prob_ce = mean(is_best), .groups = "drop") %>%
      mutate(wtp = w)
  })) %>%
    mutate(strategy = factor(strategy, levels = FRONTIER))
}

psa_df    <- read.csv("output/csv/psa_results.csv", stringsAsFactors = FALSE)
tx_counts <- read.csv("output/csv/transmission_prevention_counts.csv", stringsAsFactors = FALSE)
psa_df    <- psa_df %>% mutate(cost = cost / 1e5, qaly = qaly / 1e5)

atb_prevent  <- setNames(tx_counts$atb_cases_prevented, tx_counts$strategy)
contacts_sim <- runif(n_sim, min = n_contacts_psa_min, max = n_contacts_psa_max)

tx_savings <- expand.grid(sim = 1:n_sim, strategy = names(atb_prevent),
                          stringsAsFactors = FALSE)
tx_savings$savings <- mapply(function(s, i) {
  n_c_i  <- contacts_sim[i]
  beta_i <- n_c_i * p_tx_per_contact
  atb_i  <- max(0, atb_prevent[s])
  (atb_i * beta_i * cost_per_secondary_tb +
   atb_i * n_c_i * p_ltbi_per_contact * cost_per_secondary_ltbi) / 1e5
}, tx_savings$strategy, tx_savings$sim)

psa_tx  <- psa_df %>%
  left_join(tx_savings, by = c("sim", "strategy")) %>%
  mutate(cost = cost - coalesce(savings, 0)) %>%
  filter(strategy %in% FRONTIER)
psa_nhs <- psa_df %>% filter(strategy %in% FRONTIER)

ceac_both <- bind_rows(
  compute_ceac(psa_tx,  wtp_range) %>% mutate(perspective = "Societal (primary)"),
  compute_ceac(psa_nhs, wtp_range) %>% mutate(perspective = "NHS (robustness)")
) %>%
  mutate(
    strategy    = factor(strategy,    levels = FRONTIER),
    perspective = factor(perspective, levels = c("Societal (primary)", "NHS (robustness)"))
  )

p_ceac <- ggplot(ceac_both,
                 aes(x = wtp, y = prob_ce, colour = strategy, linetype = perspective)) +
  geom_vline(xintercept = 25000, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 35000, colour = "grey60", linewidth = 0.4) +
  annotate("text", x = 25000, y = 0.04, label = "\u00a325k",
           colour = "grey40", size = 3.5, hjust = -0.15) +
  annotate("text", x = 35000, y = 0.04, label = "\u00a335k",
           colour = "grey60", size = 3.5, hjust = -0.15) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    labels = function(x) ifelse(x == 0, "\u00a30", paste0("\u00a3", x / 1000, "k")),
    breaks = seq(0, 50000, by = 10000), limits = c(0, 50000)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = strat_cols_ceac, name = NULL) +
  scale_linetype_manual(
    values = c("Societal (primary)" = "solid", "NHS (robustness)" = "dashed"),
    labels = c("Societal (primary)" = "Societal perspective (includes transmission savings)",
               "NHS (robustness)"   = "NHS perspective"),
    name = NULL) +
  guides(
    colour   = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
    linetype = guide_legend(order = 2, keywidth = unit(2.2, "cm"),
                            override.aes = list(colour = "grey20", linewidth = 1.0))) +
  theme_minimal(base_size = 12) +
  labs(x = "Willingness-to-pay (\u00a3/QALY)", y = "Cost-effectiveness probability") +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.box = "vertical", legend.text = element_text(size = 10.5),
        legend.spacing.y = unit(2, "pt"))

ggsave("output/body/ceac_panel.tiff", p_ceac, width = 9, height = 6, dpi = 300,
       compression = "lzw")
cat("Saved: output/body/ceac_panel.tiff\n")

# =============================================================================
# SUPPLEMENTARY — PSA CE plane (psa_ce_plane_supplementary.tiff)
# Incremental cost vs incremental QALY scatter; 1,000 sims; frontier vs passive;
# NHS perspective (consistent with CEAC body figure)
# =============================================================================
cat("Generating PSA CE plane...\n")

# psa_df already loaded above (NHS perspective, per-person units)
# Compute incremental vs Passive case finding per simulation
passive_psa <- psa_nhs %>%
  filter(strategy == "Passive case finding") %>%
  select(sim, cost_ref = cost, qaly_ref = qaly)

psa_ce_plane <- psa_nhs %>%
  filter(strategy %in% FRONTIER, strategy != "Passive case finding") %>%
  left_join(passive_psa, by = "sim") %>%
  mutate(
    inc_cost = cost - cost_ref,
    inc_qaly = qaly - qaly_ref,
    strategy = factor(strategy, levels = FRONTIER)
  )

cols_ce_plane <- strat_cols_ceac[names(strat_cols_ceac) != "Passive case finding"]

# Per-facet P(CE) annotation at £25k/QALY (fraction of sims below £25k WTP line)
.pce_annot <- psa_ce_plane %>%
  group_by(strategy) %>%
  summarise(
    pce_25k = mean(inc_cost < 25000 * inc_qaly, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(label = paste0(round(pce_25k * 100, 1), "% CE at \u00a325k/QALY"))

# WTP line label positions: per-strategy 85th percentile of inc_qaly
.wtp_lbl <- psa_ce_plane %>%
  group_by(strategy) %>%
  summarise(x_ref = quantile(inc_qaly, 0.85, na.rm = TRUE), .groups = "drop") %>%
  tidyr::crossing(
    data.frame(wtp = c(25000, 35000),
               lbl = c("\u00a325k/QALY", "\u00a335k/QALY"),
               col_wtp = c("grey40", "grey60"),
               stringsAsFactors = FALSE)
  ) %>%
  mutate(y_ref = x_ref * wtp)

p_ce_plane <- ggplot(psa_ce_plane,
    aes(x = inc_qaly, y = inc_cost, colour = strategy)) +
  geom_abline(slope = 25000, intercept = 0, colour = "grey40",
              linewidth = 0.4, linetype = "solid") +
  geom_abline(slope = 35000, intercept = 0, colour = "grey60",
              linewidth = 0.4, linetype = "solid") +
  geom_text(data = .wtp_lbl,
            aes(x = x_ref, y = y_ref, label = lbl, colour = NULL),
            colour = .wtp_lbl$col_wtp,
            size = 3.0, hjust = 0, vjust = -0.3, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_text(data = .pce_annot,
            aes(label = label), x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.4, size = 3.0,
            colour = "grey30", inherit.aes = FALSE) +
  scale_colour_manual(values = cols_ce_plane, name = NULL) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0",
    paste0(ifelse(x < 0, "\u2212", ""), round(abs(x), 3)))) +
  scale_y_continuous(labels = function(x) {
    sign_str <- ifelse(x < 0, "\u2212\u00a3", "\u00a3")
    paste0(sign_str, format(round(abs(x)), big.mark = ","))
  }) +
  facet_wrap(~strategy, ncol = 3, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    x        = "Incremental QALYs vs passive case finding (per person)",
    y        = "Incremental cost vs passive case finding (\u00a3 per person)",
    title    = NULL,
    subtitle = NULL
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    strip.text       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(colour = "grey40", size = 10.5)
  )

ggsave("output/supplementary/psa_ce_plane_supplementary.tiff", p_ce_plane,
       width = 13, height = 5, dpi = 300, compression = "lzw")
cat("Saved: output/supplementary/psa_ce_plane_supplementary.tiff\n")

# =============================================================================
# FIGURE 3 — Tornado panel (body_tornado_panel.tiff) + supplementary SF10-SF12
# DSA one-way NMB deviations; top 6 parameters per strategy for panel
# =============================================================================
cat("Generating tornado diagrams...\n")

col_coughcxr <- LANCET[4]
col_symscr   <- LANCET[1]
col_parallel <- LANCET[7]

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
  "Prob: Active not treated to diagnosed"            = "Active TB: untreated \u2192 detected",
  "Prob: Active not treated to Dead"                 = "Active TB mortality (untreated)",
  # Active TB treatment
  "Prob: Active treatment to completed"              = "Active TB treatment completion",
  "Prob: Active treatment to discontinued"           = "Active TB treatment discontinuation",
  "Prob: Active treatment to lost to follow-up"      = "Active TB treatment: lost to follow-up",
  "Prob: Active treated to Dead"                     = "Active TB mortality (on treatment)",
  "Prob: Active discontinued to undiagnosed"         = "Active TB: discontinued \u2192 undetected",
  "Prob: Active discontinued to Dead"                = "Active TB mortality (discontinued)",
  # Cost
  "Cost - Active TB under treatment"                 = "Cost: active TB treatment (monthly)"
)

plot_tornado <- function(csv_file, subtitle_str, bar_col, cap_col, out_file,
                         top_n = 15, base_size = 13, save = TRUE,
                         label_width = 40, show_caption = FALSE, show_title = TRUE) {
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  base_nmb <- df$nmb_base[1]
  top_n <- min(top_n, nrow(df))
  t_df <- df %>%
    slice_head(n = top_n) %>%
    mutate(
      param_short = dplyr::coalesce(dsa_short_labels[param_label], param_label),
      param_short = ifelse(nchar(param_short) > label_width,
                           paste0(substr(param_short, 1, label_width - 3), "..."),
                           param_short),
      param_short = make.unique(param_short, sep = " "),
      param_short = factor(param_short, levels = rev(param_short)),
      nmb_lo_dev  = nmb_lo - base_nmb,
      nmb_hi_dev  = nmb_hi - base_nmb,
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
      }, n.breaks = 5) +
    theme_minimal(base_size = base_size) +
    labs(
      x        = "Change in NMB from base case (2023/24 \u00a3)",
      y        = "",
      title    = if (show_title) "One-way deterministic sensitivity analysis" else NULL,
      subtitle = subtitle_str,
      caption  = if (show_caption) paste0(
        "NMB = net monetary benefit; WTP = willingness-to-pay threshold (NICE \u00a325,000\u2013\u00a335,000/QALY); ",
        "DSA = deterministic sensitivity analysis.\n",
        "All costs in 2023/24 GBP. Parameters varied individually across 95% plausible range; ",
        "all others held at base-case values.") else NULL) +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 2),
      plot.subtitle      = element_text(color = "grey40", size = base_size - 1),
      plot.caption       = element_text(color = "grey40",
                                        size = max(base_size - 3, 7), hjust = 0),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = base_size - 1))
  if (save) {
    ggsave(out_file, p, width = 14, height = 8, dpi = 300)
    cat(sprintf("Saved: %s\n", out_file))
  }
  invisible(p)
}

# Supplementary individual tornado figures omitted: the 3-panel body figure covers all three strategies.

# 3-panel body figure (top 6 per strategy)
p_tA <- plot_tornado("output/csv/dsa_results_coughcxr.csv",
  subtitle_str = "A. Cough+CXR (TB sx)",
  bar_col = col_coughcxr, cap_col = NAVY,
  top_n = 6, base_size = 12, label_width = 36, save = FALSE,
  show_caption = FALSE, show_title = FALSE) + labs(x = "")

p_tB <- plot_tornado("output/csv/dsa_results_symscrCXR.csv",
  subtitle_str = "B. Symptom screen+CXR",
  bar_col = col_symscr, cap_col = NAVY,
  top_n = 6, base_size = 12, label_width = 36, save = FALSE,
  show_caption = FALSE, show_title = FALSE)

p_tC <- plot_tornado("output/csv/dsa_results.csv",
  subtitle_str = "C. Parallel Sx+QFT (Ultra)",
  bar_col = col_parallel, cap_col = NAVY,
  top_n = 6, base_size = 12, label_width = 36, save = FALSE,
  show_caption = FALSE, show_title = FALSE) +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/body/tornado_panel.tiff",
       (p_tA | p_tB | p_tC), width = 13, height = 5.5, dpi = 300,
       compression = "lzw")
cat("Saved: output/body/tornado_panel.tiff\n")

# =============================================================================
# FIGURE 4 — Transmission prevention counts (body_transmission_prevention_counts.tiff)
# Grouped bar: reactivations prevented + secondary active TB + secondary LTBI avoided
# =============================================================================
cat("Generating transmission prevention counts...\n")

COMP_COLS <- c(
  "Active TB reactivations prevented (direct)" = LANCET[1],
  "Secondary active TB cases avoided"          = LANCET[4],
  "Secondary LTBI cases avoided"               = LANCET[8]
)

tx_counts2 <- read.csv("output/csv/transmission_prevention_counts.csv",
                       stringsAsFactors = FALSE)
strat_order <- c("Cough+CXR (TB sx)", "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)")

dat <- tx_counts2 %>%
  filter(strategy %in% strat_order) %>%
  transmute(
    strategy       = factor(strategy, levels = rev(strat_order)),
    atb_prevented  = pmax(0, round(atb_cases_prevented)),
    secondary_atb  = pmax(0, round(secondary_cases_avoided)),
    secondary_ltbi = pmax(0, round(secondary_ltbi_avoided)),
    total          = atb_prevented + secondary_atb + secondary_ltbi
  )

long <- dat %>%
  pivot_longer(c(atb_prevented, secondary_atb, secondary_ltbi),
               names_to = "component", values_to = "cases") %>%
  mutate(component = factor(component,
    levels = c("atb_prevented", "secondary_atb", "secondary_ltbi"),
    labels = names(COMP_COLS)))

p_tx <- ggplot(long, aes(x = cases, y = strategy, fill = component)) +
  geom_col(position = "stack", width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(data = dat, aes(x = total, y = strategy,
    label = formatC(total, format = "d", big.mark = ","),
    colour = strategy), inherit.aes = FALSE,
    hjust = -0.15, size = 4.0, fontface = "bold") +
  geom_text(data = long %>% filter(cases >= 60),
    aes(x = cases, y = strategy,
        label = formatC(cases, format = "d", big.mark = ",")),
    position = position_stack(vjust = 0.5), inherit.aes = FALSE,
    colour = "white", size = 3.2) +
  scale_fill_manual(values = COMP_COLS, name = NULL) +
  scale_colour_manual(
    values = c("Cough+CXR (TB sx)"       = LANCET[4],
               "Symptom screen+CXR"      = LANCET[1],
               "Parallel Sx+QFT (Ultra)" = LANCET[7]),
    guide = "none") +
  scale_x_continuous(limits = c(0, 3200), breaks = seq(0, 3000, 500),
    expand = c(0, 0),
    labels = function(x) formatC(x, format = "d", big.mark = ",")) +
  labs(x = "Cases prevented per 100,000 migrants screened", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9.5, colour = "grey40"),
    legend.position = "bottom", legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"), legend.spacing.x = unit(0.3, "cm"),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    axis.text.y = element_text(size = 11, colour = "grey20"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, colour = "grey20"),
    plot.margin = margin(8, 40, 8, 8)) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))

ggsave("output/body/transmission_prevention_counts.tiff",
       p_tx, width = 10, height = 4.5, dpi = 300, compression = "lzw")
cat("Saved: output/body/transmission_prevention_counts.tiff\n")

# Perspective comparison figure omitted: cost shift is meaningful only for Parallel Sx+QFT
# (Ultra) (−£15.73/person); non-IGRA strategies shift <£0.01. Content is covered by the
# transmission prevention bar chart (Figure 4).

# =============================================================================
# EXCEL TABLE BUILDERS (openxlsx — replaces Python build_*.py scripts)
# =============================================================================
suppressPackageStartupMessages(library(openxlsx))

.today   <- format(Sys.Date(), "%Y-%m-%d")
.N_COH   <- 100000L
.WTP_PTS <- c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000)

# ── Shared style constructors ─────────────────────────────────────────────────
.NAVY <- "#00468B"; .YELLOW <- "#FFF2CC"; .ALT <- "#EBF3FB"
.s_hdr  <- function() createStyle(fgFill=.NAVY, fontColour="white", textDecoration="bold",
                                    halign="center", valign="center", wrapText=TRUE, border="TopBottomLeftRight", borderStyle="thin")
.s_hi   <- function() createStyle(fgFill=.YELLOW, textDecoration="bold", border="TopBottomLeftRight", borderStyle="thin", wrapText=TRUE)
.s_hin  <- function() createStyle(fgFill=.YELLOW, textDecoration="bold", border="TopBottomLeftRight", borderStyle="thin", halign="center")
.s_row  <- function() createStyle(border="TopBottomLeftRight", borderStyle="thin", wrapText=TRUE)
.s_num  <- function() createStyle(border="TopBottomLeftRight", borderStyle="thin", halign="center")
.s_alt  <- function() createStyle(fgFill=.ALT, border="TopBottomLeftRight", borderStyle="thin", wrapText=TRUE)
.s_altn <- function() createStyle(fgFill=.ALT, border="TopBottomLeftRight", borderStyle="thin", halign="center")
.s_note <- function() createStyle(fontSize=8, textDecoration="italic", fontColour="#555555", wrapText=TRUE)
.s_ttl  <- function() createStyle(fgFill=.NAVY, fontColour="white", textDecoration="bold",
                                    fontSize=11, halign="left", valign="center", wrapText=TRUE)
.s_sub  <- function() createStyle(fontSize=9, fontColour="#333333", textDecoration="italic",
                                    wrapText=TRUE, valign="top")

.DOM_LBL <- c("ref"="Reference","non-dominated"="Non-dominated","dominated"="Dominated",
               "ext_dominated"="Ext. dominated","simply dominated"="Simply dominated",
               "extendedly dominated"="Ext. dominated")
.gbp <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v)|!is.finite(v),"\u2014",paste0("\u00a3",formatC(round(v),format="d",big.mark=",")))
}
.qf  <- function(v) { v <- suppressWarnings(as.numeric(v)); ifelse(is.na(v)|!is.finite(v),"\u2014",sprintf("%.4f",v)) }
.pf  <- function(v) { v <- suppressWarnings(as.numeric(v)); ifelse(is.na(v)|!is.finite(v),"\u2014",sprintf("%.1f%%",v*100)) }
.icf <- function(v, dom) {
  if (!is.na(dom) && dom == "ref") return("Reference")
  .gbp(v)
}
.igra_lbl <- function(s) ifelse(grepl("^Parallel", s), "Yes", "No")

.style_icer_rows <- function(wb, sh, df, hdr_row, text_cols, num_cols) {
  for (ri in seq_len(nrow(df))) {
    is_f <- df$strategy[ri] %in% FRONTIER
    rr   <- hdr_row + ri
    addStyle(wb, sh, if(is_f).s_hi()  else if(ri%%2==0).s_alt()  else .s_row(),
             rows=rr, cols=text_cols, gridExpand=TRUE, stack=FALSE)
    addStyle(wb, sh, if(is_f).s_hin() else if(ri%%2==0).s_altn() else .s_num(),
             rows=rr, cols=num_cols,  gridExpand=TRUE, stack=FALSE)
  }
}

# Write one ICER results sheet (common structure for ST3a/b + icer_table)
.write_icer_sheet <- function(wb, sh, csv_path, title, subtitle="", has_detection=TRUE) {
  df <- tryCatch(read.csv(csv_path,stringsAsFactors=FALSE,check.names=FALSE),error=function(e)NULL)
  if (is.null(df)) { cat(sprintf("  [SKIP] %s not found\n", csv_path)); return() }

  addWorksheet(wb, sh)

  # Sort: ref -> non-dominated -> rest
  .ok <- c("ref"=0,"non-dominated"=1); df$sk <- .ok[df$dominance]; df$sk[is.na(df$sk)] <- 2
  df <- df[order(df$sk, suppressWarnings(as.numeric(df$cost_per_person))), ]; df$sk <- NULL

  ref_d <- suppressWarnings(as.numeric(df$deaths[df$dominance=="ref"][1]))
  df$d_avert <- ifelse(is.na(df$deaths)|is.na(ref_d), "\u2014",
                        formatC(round(ref_d - suppressWarnings(as.numeric(df$deaths))),format="d",big.mark=","))

  hdrs <- c("Strategy","Cost/person (\u00a3)","QALYs/person","Inc. cost (\u00a3)","Inc. QALYs",
             "Pairwise ICER (\u00a3/QALY)","Sequential ICER (\u00a3/QALY)","Dominance",
             "Deaths averted/100k")
  cw   <- c(34,14,12,14,12,18,18,16,15)
  if (has_detection && "ltbi_detected" %in% names(df)) {
    hdrs <- c(hdrs,"LTBI detected/100k","Active TB detected/100k"); cw <- c(cw,16,18)
  }
  nc <- length(hdrs)

  # Title
  writeData(wb,sh,title,startRow=1,startCol=1,colNames=FALSE)
  addStyle(wb,sh,.s_ttl(),rows=1,cols=1); mergeCells(wb,sh,cols=1:nc,rows=1)
  setRowHeights(wb,sh,rows=1,heights=20)

  hdr_row <- 2
  if (nchar(subtitle)>0) {
    writeData(wb,sh,subtitle,startRow=2,startCol=1,colNames=FALSE)
    addStyle(wb,sh,.s_sub(),rows=2,cols=1); mergeCells(wb,sh,cols=1:nc,rows=2)
    setRowHeights(wb,sh,rows=2,heights=24); hdr_row <- 3
  }
  writeData(wb,sh,as.data.frame(t(hdrs)),startRow=hdr_row,startCol=1,colNames=FALSE)
  addStyle(wb,sh,.s_hdr(),rows=hdr_row,cols=1:nc,gridExpand=TRUE)
  setRowHeights(wb,sh,rows=hdr_row,heights=30)

  display <- data.frame(
    s1 = df$strategy,
    s2 = sapply(df$cost_per_person,.gbp),
    s3 = sapply(df$qaly_per_person,.qf),
    s4 = mapply(function(v) .gbp(suppressWarnings(as.numeric(v))/.N_COH), df$inc_cost),
    s5 = mapply(function(v) { x<-suppressWarnings(as.numeric(v))/.N_COH; ifelse(is.na(x)|!is.finite(x),"\u2014",sprintf("%.6f",x)) }, df$inc_qaly),
    s6 = mapply(.icf, df$icer, df$dominance),
    s7 = mapply(.icf, df$sequential_icer, df$dominance),
    s8 = ifelse(is.na(.DOM_LBL[df$dominance]), df$dominance, .DOM_LBL[df$dominance]),
    s9 = df$d_avert,
    stringsAsFactors=FALSE, check.names=FALSE
  )
  if (has_detection && "ltbi_detected" %in% names(df)) {
    display$s10 <- sapply(df$ltbi_detected,.gbp)
    display$s11 <- sapply(df$active_detected,.gbp)
  }
  writeData(wb,sh,display,startRow=hdr_row+1,startCol=1,colNames=FALSE)
  setRowHeights(wb,sh,rows=(hdr_row+1):(hdr_row+nrow(display)),heights=17)
  .style_icer_rows(wb,sh,df,hdr_row,c(1,8),setdiff(1:nc,c(1,8)))
  setColWidths(wb,sh,cols=1:nc,widths=cw)

  nr <- hdr_row+nrow(display)+2
  writeData(wb,sh,"Highlighted rows = efficient frontier. Sequential ICER = incremental cost/QALY vs next cheapest non-dominated strategy. All costs 2023/24 GBP. 3.5% discount. NICE WTP \u00a325,000\u2013\u00a335,000/QALY.",startRow=nr,startCol=1,colNames=FALSE)
  addStyle(wb,sh,.s_note(),rows=nr,cols=1); mergeCells(wb,sh,cols=1:nc,rows=nr)
  setRowHeights(wb,sh,rows=nr,heights=36)
}

# Simple SA sheet (tabular dump with description header)
.write_sa_sheet <- function(wb, sh, csv_path, description) {
  df <- tryCatch(read.csv(csv_path,stringsAsFactors=FALSE,check.names=FALSE),error=function(e)NULL)
  if (is.null(df)) { cat(sprintf("  [SKIP] %s\n",csv_path)); return() }
  addWorksheet(wb, sh)
  writeData(wb,sh,description,startRow=1,startCol=1,colNames=FALSE)
  addStyle(wb,sh,.s_sub(),rows=1,cols=1); mergeCells(wb,sh,cols=1:ncol(df),rows=1)
  setRowHeights(wb,sh,rows=1,heights=28)
  writeDataTable(wb,sh,df,startRow=2,startCol=1,tableStyle="TableStyleLight9",withFilter=TRUE)
  setColWidths(wb,sh,cols=seq_len(ncol(df)),widths="auto")
  if ("strategy" %in% names(df)) {
    for (ri in seq_len(nrow(df))) {
      if (df$strategy[ri] %in% FRONTIER)
        addStyle(wb,sh,.s_hi(),rows=ri+2,cols=1:ncol(df),gridExpand=TRUE,stack=FALSE)
    }
  }
}

# =============================================================================
# output/body/icer_table.xlsx — review table
# =============================================================================
cat("Building output/body/icer_table.xlsx...\n")
.wb_it <- createWorkbook(); modifyBaseFont(.wb_it,fontSize=10,fontName="Calibri")
.write_icer_sheet(.wb_it,"NHS perspective","output/csv/icer_nhs_perspective.csv",
  "ICER table \u2014 NHS perspective (robustness check)",
  paste0("NHS direct costs only \u00b7 no transmission savings \u00b7 55-year horizon \u00b7 n=100,000 \u00b7 2023/24 GBP \u00b7 ",.today))
.write_icer_sheet(.wb_it,"Societal perspective","output/csv/icer_societal_perspective.csv",
  "ICER table \u2014 Societal perspective (primary analysis)",
  paste0("Includes secondary TB prevention savings (\u00a320.96/person) \u00b7 \u03b2=0.205 \u00b7 \u00a36,055/secondary case \u00b7 ",.today))
saveWorkbook(.wb_it,"output/body/icer_table.xlsx",overwrite=TRUE)
cat("Saved: output/body/icer_table.xlsx\n")

# =============================================================================
# output/supplementary/CE_Results_supplementary.xlsx — NHS + societal + PSA summary
# =============================================================================
cat("Building output/supplementary/CE_Results_supplementary.xlsx...\n")
.wb3 <- createWorkbook(); modifyBaseFont(.wb3,fontSize=10,fontName="Calibri")
.write_icer_sheet(.wb3,"NHS perspective (robustness)","output/csv/icer_nhs_perspective.csv",
  "ST3a. Cost-effectiveness results \u2014 all 43 strategies (NHS perspective, robustness check)",
  paste0("NHS direct costs \u00b7 no transmission savings \u00b7 55-year horizon \u00b7 n=100,000 \u00b7 2023/24 GBP \u00b7 ",.today))
.write_icer_sheet(.wb3,"Primary \u2014 Transmission","output/csv/icer_societal_perspective.csv",
  "ST3b. Cost-effectiveness results \u2014 primary analysis: transmission prevention (societal perspective)",
  paste0("Transmission savings included (\u03b2=0.205 \u00b7 \u00a36,055/secondary case) \u00b7 2023/24 GBP \u00b7 ",.today))

.psa_ci <- tryCatch(read.csv("output/csv/icer_confidence_intervals.csv",stringsAsFactors=FALSE),error=function(e)NULL)
if (!is.null(.psa_ci)) {
  addWorksheet(.wb3,"PSA summary")
  .ph <- c("Strategy","Median ICER (\u00a3/QALY)","95% CrI low","95% CrI high",
            "Mean inc. cost/person (\u00a3)","Mean inc. QALYs/person",
            "P(CE) at \u00a325,000","P(CE) at \u00a335,000","P(dominated)")
  writeData(.wb3,"PSA summary","ST3c. Probabilistic sensitivity analysis \u2014 ICER credible intervals (1,000 simulations)",
            startRow=1,startCol=1,colNames=FALSE)
  addStyle(.wb3,"PSA summary",.s_ttl(),rows=1,cols=1); mergeCells(.wb3,"PSA summary",cols=1:9,rows=1); setRowHeights(.wb3,"PSA summary",rows=1,heights=20)
  writeData(.wb3,"PSA summary",paste0("1,000 MC simulations \u00b7 Beta/Gamma distributions \u00b7 2023/24 GBP \u00b7 ",.today),startRow=2,startCol=1,colNames=FALSE)
  addStyle(.wb3,"PSA summary",.s_sub(),rows=2,cols=1); mergeCells(.wb3,"PSA summary",cols=1:9,rows=2); setRowHeights(.wb3,"PSA summary",rows=2,heights=18)
  writeData(.wb3,"PSA summary",as.data.frame(t(.ph)),startRow=3,startCol=1,colNames=FALSE)
  addStyle(.wb3,"PSA summary",.s_hdr(),rows=3,cols=1:9,gridExpand=TRUE); setRowHeights(.wb3,"PSA summary",rows=3,heights=30)
  setColWidths(.wb3,"PSA summary",cols=1:9,widths=c(34,18,14,14,22,22,16,16,12))
  .ps <- .psa_ci[order(suppressWarnings(as.numeric(.psa_ci$median_icer))), ]
  .pd <- data.frame(s=.ps$strategy, a=sapply(.ps$median_icer,.gbp), b=sapply(.ps$icer_lo,.gbp),
                     c=sapply(.ps$icer_hi,.gbp), d=sapply(.ps$mean_inc_cost_pp,.gbp),
                     e=sapply(.ps$mean_inc_qaly_pp,.qf), f=sapply(.ps$prob_ce_25k,.pf),
                     g=sapply(.ps$prob_ce_35k,.pf), h=sapply(.ps$prob_dominated,.pf),
                     stringsAsFactors=FALSE, check.names=FALSE)
  writeData(.wb3,"PSA summary",.pd,startRow=4,startCol=1,colNames=FALSE)
  setRowHeights(.wb3,"PSA summary",rows=4:(3+nrow(.ps)),heights=17)
  for (ri in seq_len(nrow(.ps))) {
    rr <- 3+ri; is_f <- .ps$strategy[ri] %in% FRONTIER
    addStyle(.wb3,"PSA summary",if(is_f).s_hi() else if(ri%%2==0).s_alt() else .s_row(),rows=rr,cols=1,stack=FALSE)
    addStyle(.wb3,"PSA summary",if(is_f).s_hin() else if(ri%%2==0).s_altn() else .s_num(),rows=rr,cols=2:9,gridExpand=TRUE,stack=FALSE)
  }
}
saveWorkbook(.wb3,"output/supplementary/CE_Results_supplementary.xlsx",overwrite=TRUE)
cat("Saved: output/supplementary/CE_Results_supplementary.xlsx\n")

# =============================================================================
# output/supplementary/Sensitivity_Analyses_supplementary.xlsx — all SA sheets
# =============================================================================
cat("Building output/supplementary/Sensitivity_Analyses_supplementary.xlsx...\n")
.wb4 <- createWorkbook(); modifyBaseFont(.wb4,fontSize=10,fontName="Calibri")
.write_sa_sheet(.wb4,"PSA \u2014 ICER credible intervals","output/csv/icer_confidence_intervals.csv",
  paste0("ST4a \u00b7 PSA base case \u2014 ICER 95% CrI and P(CE) for all 43 strategies \u00b7 1,000 sims \u00b7 Beta/Gamma \u00b7 ",.today))
# ST4b: raw PSA simulations (43,001 rows) replaced by a repository note; data in output/csv/psa_results.csv.
addWorksheet(.wb4,"PSA \u2014 Raw simulations (repo)")
writeData(.wb4,"PSA \u2014 Raw simulations (repo)",
  data.frame(Note = paste0(
    "ST4b \u00b7 PSA raw results (1,000 simulations \u00d7 43 strategies) are available in the ",
    "project repository (output/csv/psa_results.csv). The full dataset (43,000 rows) is not ",
    "reproduced here. Key summary statistics are reported in ST4a and body Table 3.")),
  startRow=1,startCol=1,colNames=FALSE)
.write_sa_sheet(.wb4,"PSA \u2014 Transmission","output/csv/transmission_psa_summary.csv",
  paste0("ST4c \u00b7 Societal perspective PSA \u2014 frontier strategies \u00b7 contacts~Uniform(5,15)/case/year \u00b7 \u00a36,055/secondary case \u00b7 ",.today))
.write_sa_sheet(.wb4,"DSA \u2014 Cough+CXR (TB sx)","output/csv/dsa_results_coughcxr.csv",
  paste0("ST4d \u00b7 One-way DSA \u2014 Cough+CXR (TB sx) vs No Screening \u00b7 NMB at \u00a325,000/QALY \u00b7 ",.today))
.write_sa_sheet(.wb4,"DSA \u2014 Symptom screen+CXR","output/csv/dsa_results_symscrCXR.csv",
  paste0("ST4e \u00b7 One-way DSA \u2014 Symptom screen+CXR vs No Screening \u00b7 NMB at \u00a325,000/QALY \u00b7 ",.today))
.write_sa_sheet(.wb4,"DSA \u2014 Parallel Sx+QFT (Ultra)","output/csv/dsa_results.csv",
  paste0("ST4f \u00b7 One-way DSA \u2014 Parallel Sx+QFT (Ultra) vs No Screening \u00b7 NMB at \u00a325,000/QALY \u00b7 ",.today))
.write_sa_sheet(.wb4,"DSA \u2014 Transmission contacts","output/csv/transmission_dsa_results.csv",
  paste0("ST4g \u00b7 Primary analysis DSA \u2014 close contacts per active TB case \u00b7 \u00a36,055/secondary case \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 LTBI completion","output/csv/ltbi_completion_sensitivity.csv",
  paste0("ST4h \u00b7 Structural SA \u2014 LTBI treatment completion \u00b7 Base 74.8% vs real-world 55.5% \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 LTBI prevalence","output/csv/ltbi_prevalence_sensitivity.csv",
  paste0("ST4i \u00b7 Structural SA \u2014 LTBI prevalence at entry \u00b7 Base 17.8% vs 15.1% \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 Active TB prevalence","output/csv/active_tb_prev_sensitivity.csv",
  paste0("ST4j \u00b7 Structural SA \u2014 Active TB prevalence at entry \u00b7 Base 1.0%, 0.44%, 0.215% \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 IGRA uptake","output/csv/igra_uptake_sensitivity.csv",
  paste0("ST4k \u00b7 Structural SA \u2014 IGRA programme uptake \u00b7 70%\u2013100% \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 IGRA specificity","output/csv/igra_specificity_sensitivity.csv",
  paste0("ST4l \u00b7 Structural SA \u2014 IGRA specificity \u00b7 QFT base 0.96; T-SPOT base 0.93; range 0.90\u20131.00 \u00b7 ",.today))
.write_sa_sheet(.wb4,"SA \u2014 LTBI treatment efficacy","output/csv/ltbi_efficacy_sensitivity.csv",
  paste0("ST4n \u00b7 Structural SA \u2014 LTBI treatment efficacy \u00b7 Base HR 0.14 (Berrocal-Almanza 2022, 0.16% lifetime failure) vs NICE NG33 residual risk 20\u201340% \u00b7 At 20\u201340% residual risk Parallel Sx+QFT (Ultra) is simply dominated; frontier reverts to non-IGRA strategies \u00b7 ",.today))
# ST4m: 7 representative values extracted from full sweep (£0.50–£60 in £0.50 steps).
# Full dataset retained in output/csv/igra_programme_cost_sensitivity.csv.
.prog_full <- tryCatch(read.csv("output/csv/igra_programme_cost_sensitivity.csv",stringsAsFactors=FALSE),error=function(e)NULL)
if(!is.null(.prog_full)){
  .prog_summary <- .prog_full[
    .prog_full$strategy == "Parallel Sx+QFT (Ultra)" &
    .prog_full$prog_cost_type == "allsx" &
    .prog_full$prog_cost_pp %in% c(0, 10, 20, 30, 40, 50, 60), ]
  .prog_summary <- .prog_summary[order(.prog_summary$prog_cost_pp), ]
  addWorksheet(.wb4,"SA \u2014 IGRA programme cost")
  .ttl_prog <- paste0("ST4m \u00b7 Structural SA \u2014 IGRA programme delivery cost per person \u00b7 Parallel Sx+QFT (Ultra), symptom-screen pathway \u00b7 7 representative values shown (full sweep \u00a30.50\u2013\u00a360 in \u00a30.50 steps in project repository) \u00b7 ",.today)
  writeData(.wb4,"SA \u2014 IGRA programme cost",.ttl_prog,startRow=1,startCol=1,colNames=FALSE)
  addStyle(.wb4,"SA \u2014 IGRA programme cost",.s_ttl(),rows=1,cols=1)
  mergeCells(.wb4,"SA \u2014 IGRA programme cost",cols=1:4,rows=1)
  setRowHeights(.wb4,"SA \u2014 IGRA programme cost",rows=1,heights=28)
  .prog_hdr <- c("Strategy","Pathway type","Programme cost per person (\u00a3)","Sequential ICER vs passive (\u00a3/QALY)")
  writeData(.wb4,"SA \u2014 IGRA programme cost",as.data.frame(t(.prog_hdr)),startRow=2,startCol=1,colNames=FALSE)
  addStyle(.wb4,"SA \u2014 IGRA programme cost",createStyle(textDecoration="bold",halign="center",border="TopBottomLeftRight",borderStyle="thin"),rows=2,cols=1:4,gridExpand=TRUE)
  .prog_out <- data.frame(
    strategy      = .prog_summary$strategy,
    type          = .prog_summary$prog_cost_type,
    cost_pp       = paste0("\u00a3", round(.prog_summary$prog_cost_pp)),
    icer          = sapply(.prog_summary$icer, .gbp),
    stringsAsFactors=FALSE)
  writeData(.wb4,"SA \u2014 IGRA programme cost",.prog_out,startRow=3,startCol=1,colNames=FALSE)
  setColWidths(.wb4,"SA \u2014 IGRA programme cost",cols=1:4,widths=c(28,14,24,28))
  setRowHeights(.wb4,"SA \u2014 IGRA programme cost",rows=3:(2+nrow(.prog_out)),heights=17)
  fn_prog <- nrow(.prog_out)+4
  writeData(.wb4,"SA \u2014 IGRA programme cost","Base case ICER (£0 programme cost) = £14,344/QALY. ICER increases by ~£722/QALY per £10 programme delivery cost increment. Frontier composition unchanged across full cost range (£0.50\u2013£60/person).",startRow=fn_prog,startCol=1,colNames=FALSE)
  addStyle(.wb4,"SA \u2014 IGRA programme cost",.s_note(),rows=fn_prog,cols=1)
  mergeCells(.wb4,"SA \u2014 IGRA programme cost",cols=1:4,rows=fn_prog)
  setRowHeights(.wb4,"SA \u2014 IGRA programme cost",rows=fn_prog,heights=36)
}
saveWorkbook(.wb4,"output/supplementary/Sensitivity_Analyses_supplementary.xlsx",overwrite=TRUE)
cat("Saved: output/supplementary/Sensitivity_Analyses_supplementary.xlsx\n")

# =============================================================================
# output/supplementary/PSA_Results_supplementary.xlsx — CEAF P(CE) at WTP thresholds
# =============================================================================
cat("Building output/supplementary/PSA_Results_supplementary.xlsx...\n")
.psa_raw <- tryCatch(read.csv("output/csv/psa_results.csv",stringsAsFactors=FALSE),error=function(e)NULL)
.tx_prev <- tryCatch(read.csv("output/csv/transmission_prevention_counts.csv",stringsAsFactors=FALSE),error=function(e)NULL)

if (!is.null(.psa_raw) && !is.null(.tx_prev)) {
  .cfg_raw  <- read.csv("input/config.csv",stringsAsFactors=FALSE)
  .cfg      <- setNames(suppressWarnings(as.numeric(.cfg_raw$value)), .cfg_raw$parameter)
  .atb_prev <- setNames(.tx_prev$atb_cases_prevented, .tx_prev$strategy)
  .p_tx <- .cfg["p_tx_per_contact"]; .c_tb <- .cfg["cost_secondary_tb"]
  .p_lt <- .cfg["p_ltbi_per_contact"]; .c_lt <- .cfg["cost_secondary_ltbi"]
  .n_lo <- .cfg["n_contacts_base"];  .n_hi <- .cfg["n_contacts_dsa_high"]

  set.seed(20260318)
  .N_SIM   <- max(.psa_raw$sim)
  .nc_draw <- runif(.N_SIM, .n_lo, .n_hi)

  .psa_f   <- .psa_raw[.psa_raw$strategy %in% FRONTIER, ]
  .psa_f$cost_pp <- .psa_f$cost / .N_COH
  .psa_f$qaly_pp <- .psa_f$qaly / .N_COH

  .ceaf <- function(perspective) {
    matrix(sapply(.WTP_PTS, function(wtp) {
      counts <- integer(length(FRONTIER)); names(counts) <- FRONTIER
      for (s in seq_len(.N_SIM)) {
        rows  <- .psa_f[.psa_f$sim == s, ]
        costs <- rows$cost_pp
        if (perspective == "societal") {
          nc_i   <- .nc_draw[s]
          saving <- (.atb_prev[rows$strategy] * nc_i * .p_tx * .c_tb +
                     .atb_prev[rows$strategy] * nc_i * .p_lt * .c_lt) / .N_COH
          costs  <- costs - saving
        }
        nmb  <- rows$qaly_pp * wtp - costs
        best <- rows$strategy[which.max(nmb)]
        if (!is.na(best) && best %in% FRONTIER) counts[best] <- counts[best] + 1L
      }
      counts / .N_SIM
    }), nrow=length(FRONTIER), ncol=length(.WTP_PTS), dimnames=list(FRONTIER, .WTP_PTS))
  }

  cat("  Computing societal CEAF P(CE) across 11 WTP points...\n")
  .soc_m <- .ceaf("societal")
  cat("  Computing NHS CEAF P(CE)...\n")
  .nhs_m <- .ceaf("nhs")

  .wb5 <- createWorkbook(); modifyBaseFont(.wb5,fontSize=10,fontName="Calibri")

  .make_pce_sh <- function(wb, sh, mat, title) {
    addWorksheet(wb, sh)
    setColWidths(wb, sh, cols=1:5, widths=c(12,24,24,24,24))
    writeData(wb,sh,title,startRow=1,startCol=1,colNames=FALSE)
    addStyle(wb,sh,.s_ttl(),rows=1,cols=1); mergeCells(wb,sh,cols=1:5,rows=1); setRowHeights(wb,sh,rows=1,heights=28)
    writeData(wb,sh,as.data.frame(t(c("WTP (\u00a3/QALY)",FRONTIER))),startRow=2,startCol=1,colNames=FALSE)
    addStyle(wb,sh,.s_hdr(),rows=2,cols=1:5,gridExpand=TRUE); setRowHeights(wb,sh,rows=2,heights=28)
    for (wi in seq_along(.WTP_PTS)) {
      rr <- wi+2; wtp <- .WTP_PTS[wi]; is_n <- wtp %in% c(25000,35000)
      bg <- if(is_n)"#FFF3CC" else "white"
      wlbl <- if(wtp==0)"\u00a30" else paste0("\u00a3",formatC(wtp,format="d",big.mark=","))
      writeData(wb,sh,wlbl,startRow=rr,startCol=1,colNames=FALSE)
      addStyle(wb,sh,createStyle(fgFill=bg,halign="center",border="TopBottomLeftRight", borderStyle="thin",textDecoration=if(is_n)"bold" else NULL),rows=rr,cols=1,stack=FALSE)
      for (si in seq_along(FRONTIER)) {
        val <- mat[si,wi]
        writeData(wb,sh,round(val,4),startRow=rr,startCol=si+1,colNames=FALSE)
        addStyle(wb,sh,createStyle(fgFill=bg,numFmt="0.0%",halign="center",border="TopBottomLeftRight", borderStyle="thin",
                                    textDecoration=if(is_n&&val>=0.5)"bold" else NULL,
                                    fontColour=if(is_n&&val>=0.5)"#AD002A" else "black"),
                 rows=rr,cols=si+1,stack=FALSE)
      }
      setRowHeights(wb,sh,rows=rr,heights=18)
    }
    fn <- length(.WTP_PTS)+4
    writeData(wb,sh,"P(CE) = proportion of 1,000 PSA simulations where strategy has highest NMB (CEAF approach). Frontier strategies only. Yellow rows = NICE WTP thresholds (\u00a325,000 and \u00a335,000/QALY). Bold red = P(CE) \u2265 50%.",startRow=fn,startCol=1,colNames=FALSE)
    addStyle(wb,sh,.s_note(),rows=fn,cols=1); mergeCells(wb,sh,cols=1:5,rows=fn); setRowHeights(wb,sh,rows=fn,heights=48)
  }

  .make_pce_sh(.wb5,"Societal",.soc_m,
    paste0("ST5a. P(cost-effective) by WTP \u2014 societal perspective (primary; transmission savings included) \u00b7 ",.today))
  .make_pce_sh(.wb5,"NHS",.nhs_m,
    paste0("ST5b. P(cost-effective) by WTP \u2014 NHS perspective (robustness; no transmission savings) \u00b7 ",.today))

  # Delta sheet
  addWorksheet(.wb5,"Comparison (\u0394 societal\u2212NHS)")
  setColWidths(.wb5,"Comparison (\u0394 societal\u2212NHS)",cols=1:5,widths=c(12,24,24,24,24))
  writeData(.wb5,"Comparison (\u0394 societal\u2212NHS)",
    paste0("ST5c. \u0394P(CE): societal minus NHS \u2014 effect of transmission savings assumption \u00b7 ",.today),
    startRow=1,startCol=1,colNames=FALSE)
  addStyle(.wb5,"Comparison (\u0394 societal\u2212NHS)",.s_ttl(),rows=1,cols=1)
  mergeCells(.wb5,"Comparison (\u0394 societal\u2212NHS)",cols=1:5,rows=1); setRowHeights(.wb5,"Comparison (\u0394 societal\u2212NHS)",rows=1,heights=28)
  writeData(.wb5,"Comparison (\u0394 societal\u2212NHS)",as.data.frame(t(c("WTP (\u00a3/QALY)",FRONTIER))),startRow=2,startCol=1,colNames=FALSE)
  addStyle(.wb5,"Comparison (\u0394 societal\u2212NHS)",.s_hdr(),rows=2,cols=1:5,gridExpand=TRUE); setRowHeights(.wb5,"Comparison (\u0394 societal\u2212NHS)",rows=2,heights=28)
  for (wi in seq_along(.WTP_PTS)) {
    rr <- wi+2; wtp <- .WTP_PTS[wi]; is_n <- wtp %in% c(25000,35000)
    bg <- if(is_n)"#FFF3CC" else "white"
    wlbl <- if(wtp==0)"\u00a30" else paste0("\u00a3",formatC(wtp,format="d",big.mark=","))
    writeData(.wb5,"Comparison (\u0394 societal\u2212NHS)",wlbl,startRow=rr,startCol=1,colNames=FALSE)
    addStyle(.wb5,"Comparison (\u0394 societal\u2212NHS)",createStyle(fgFill=bg,halign="center",border="TopBottomLeftRight", borderStyle="thin",textDecoration=if(is_n)"bold" else NULL),rows=rr,cols=1,stack=FALSE)
    for (si in seq_along(FRONTIER)) {
      d <- .soc_m[si,wi] - .nhs_m[si,wi]
      writeData(.wb5,"Comparison (\u0394 societal\u2212NHS)",round(d,4),startRow=rr,startCol=si+1,colNames=FALSE)
      addStyle(.wb5,"Comparison (\u0394 societal\u2212NHS)",createStyle(fgFill=bg,numFmt="+0.0%;-0.0%;0.0%",halign="center",border="TopBottomLeftRight", borderStyle="thin",textDecoration=if(is_n)"bold" else NULL,fontColour=if(d>0.001)"#005000" else if(d< -0.001)"#AD002A" else "black"),rows=rr,cols=si+1,stack=FALSE)
    }
    setRowHeights(.wb5,"Comparison (\u0394 societal\u2212NHS)",rows=rr,heights=18)
  }
  saveWorkbook(.wb5,"output/supplementary/PSA_Results_supplementary.xlsx",overwrite=TRUE)
  cat("Saved: output/supplementary/PSA_Results_supplementary.xlsx\n")
} else {
  cat("  [SKIP] psa_results.csv or transmission_prevention_counts.csv not found\n")
}

# =============================================================================
# output/body/ResultsTable_paper_YYYY-MM-DD.xlsx — body tables 2 and 3
# =============================================================================
cat("Building output/body/ResultsTable_paper_",.today,".xlsx...\n",sep="")
.icer_soc <- tryCatch(read.csv("output/csv/icer_societal_perspective.csv",stringsAsFactors=FALSE,check.names=FALSE),error=function(e)NULL)

if (!is.null(.icer_soc)) {
  .wbr <- createWorkbook(); modifyBaseFont(.wbr,fontSize=11,fontName="Calibri")
  .om <- c("ref"=0,"non-dominated"=1,"extendedly dominated"=2,"simply dominated"=3)
  .icer_soc$sk <- .om[.icer_soc$dominance]; .icer_soc$sk[is.na(.icer_soc$sk)] <- 4
  .icer_soc <- .icer_soc[order(.icer_soc$sk,suppressWarnings(as.numeric(.icer_soc$cost_per_person))), ]
  .icer_soc$sk <- NULL

  # ── Body Table 2: Frontier strategies ──────────────────────────────────────
  addWorksheet(.wbr,"Body Table 2 \u2014 Frontier")
  setColWidths(.wbr,"Body Table 2 \u2014 Frontier",cols=1:7,widths=c(28,12,13,14,16,17,16))
  writeData(.wbr,"Body Table 2 \u2014 Frontier","Table 2: Base-case cost-effectiveness results \u2014 efficient frontier",startRow=1,startCol=1,colNames=FALSE)
  addStyle(.wbr,"Body Table 2 \u2014 Frontier",createStyle(textDecoration="bold",fontSize=13),rows=1,cols=1); mergeCells(.wbr,"Body Table 2 \u2014 Frontier",cols=1:7,rows=1); setRowHeights(.wbr,"Body Table 2 \u2014 Frontier",rows=1,heights=20)
  writeData(.wbr,"Body Table 2 \u2014 Frontier",paste0("55-year lifetime horizon \u00b7 n=100,000 \u00b7 3.5% discount \u00b7 Societal perspective (primary) \u00b7 2023/24 GBP \u00b7 ",.today),startRow=2,startCol=1,colNames=FALSE)
  addStyle(.wbr,"Body Table 2 \u2014 Frontier",createStyle(textDecoration="italic",fontSize=10),rows=2,cols=1); mergeCells(.wbr,"Body Table 2 \u2014 Frontier",cols=1:7,rows=2); setRowHeights(.wbr,"Body Table 2 \u2014 Frontier",rows=2,heights=14)
  .t2h <- c("Strategy","+IGRA (LTBI)","Cost per person (\u00a3)","QALYs per person","Incremental cost vs passive (\u00a3)","Incremental QALYs vs passive","Sequential ICER (\u00a3/QALY)")
  writeData(.wbr,"Body Table 2 \u2014 Frontier",as.data.frame(t(.t2h)),startRow=3,startCol=1,colNames=FALSE)
  addStyle(.wbr,"Body Table 2 \u2014 Frontier",createStyle(textDecoration="bold",halign="center",border="TopBottomLeftRight", borderStyle="thin",wrapText=TRUE),rows=3,cols=1:7,gridExpand=TRUE); setRowHeights(.wbr,"Body Table 2 \u2014 Frontier",rows=3,heights=36)
  .frows <- .icer_soc[.icer_soc$dominance %in% c("ref","non-dominated"), ]
  for (ri in seq_len(nrow(.frows))) {
    r  <- .frows[ri, ]; rr <- ri+3
    ic <- suppressWarnings(as.numeric(r$inc_cost))/.N_COH
    iq <- suppressWarnings(as.numeric(r$inc_qaly))/.N_COH
    sq <- if(!is.na(r$dominance)&&r$dominance=="ref") "Reference" else .gbp(as.numeric(r$sequential_icer))
    rv <- c(r$strategy,.igra_lbl(r$strategy),.gbp(as.numeric(r$cost_per_person)),.qf(as.numeric(r$qaly_per_person)),
            ifelse(is.na(ic)|!is.finite(ic),"\u2014",.gbp(ic)),
            ifelse(is.na(iq)|!is.finite(iq),"\u2014",sprintf("%.6f",iq)),sq)
    writeData(.wbr,"Body Table 2 \u2014 Frontier",as.data.frame(t(rv)),startRow=rr,startCol=1,colNames=FALSE)
    addStyle(.wbr,"Body Table 2 \u2014 Frontier",createStyle(textDecoration="bold",border="TopBottomLeftRight", borderStyle="thin",wrapText=TRUE),rows=rr,cols=1,stack=FALSE)
    addStyle(.wbr,"Body Table 2 \u2014 Frontier",createStyle(textDecoration="bold",border="TopBottomLeftRight", borderStyle="thin",halign="center"),rows=rr,cols=2:7,gridExpand=TRUE,stack=FALSE)
    setRowHeights(.wbr,"Body Table 2 \u2014 Frontier",rows=rr,heights=18)
  }

  # ── Body Table 3: PSA frontier ─────────────────────────────────────────────
  if (!is.null(.psa_ci)) {
    addWorksheet(.wbr,"Body Table 3 \u2014 PSA Frontier")
    setColWidths(.wbr,"Body Table 3 \u2014 PSA Frontier",cols=1:7,widths=c(28,12,14,14,14,16,16))
    writeData(.wbr,"Body Table 3 \u2014 PSA Frontier","Table 3: PSA results \u2014 efficient frontier strategies",startRow=1,startCol=1,colNames=FALSE)
    addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",createStyle(textDecoration="bold",fontSize=13),rows=1,cols=1); mergeCells(.wbr,"Body Table 3 \u2014 PSA Frontier",cols=1:7,rows=1); setRowHeights(.wbr,"Body Table 3 \u2014 PSA Frontier",rows=1,heights=20)
    writeData(.wbr,"Body Table 3 \u2014 PSA Frontier",paste0("1,000 MC simulations \u00b7 Societal perspective \u00b7 ICER vs passive case finding \u00b7 ",.today),startRow=2,startCol=1,colNames=FALSE)
    addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",createStyle(textDecoration="italic",fontSize=10),rows=2,cols=1); mergeCells(.wbr,"Body Table 3 \u2014 PSA Frontier",cols=1:7,rows=2); setRowHeights(.wbr,"Body Table 3 \u2014 PSA Frontier",rows=2,heights=14)
    .t3h <- c("Strategy","+IGRA (LTBI)","Median ICER (\u00a3/QALY)","95% CrI lower","95% CrI upper","P(CE) at \u00a325,000/QALY","P(CE) at \u00a335,000/QALY")
    writeData(.wbr,"Body Table 3 \u2014 PSA Frontier",as.data.frame(t(.t3h)),startRow=3,startCol=1,colNames=FALSE)
    addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",createStyle(textDecoration="bold",halign="center",border="TopBottomLeftRight", borderStyle="thin",wrapText=TRUE),rows=3,cols=1:7,gridExpand=TRUE); setRowHeights(.wbr,"Body Table 3 \u2014 PSA Frontier",rows=3,heights=36)
    for (ri in seq_along(FRONTIER)) {
      st <- FRONTIER[ri]; rr <- ri+3
      ps <- .psa_ci[.psa_ci$strategy==st, ]
      rv <- if(st=="Passive case finding") c(st,.igra_lbl(st),rep("Reference",5))
            else if(nrow(ps)>0) c(st,.igra_lbl(st),.gbp(ps$median_icer[1]),.gbp(ps$icer_lo[1]),.gbp(ps$icer_hi[1]),.pf(ps$prob_ce_25k[1]),.pf(ps$prob_ce_35k[1]))
            else c(st,.igra_lbl(st),rep("\u2014",5))
      writeData(.wbr,"Body Table 3 \u2014 PSA Frontier",as.data.frame(t(rv)),startRow=rr,startCol=1,colNames=FALSE)
      addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",createStyle(textDecoration="bold",border="TopBottomLeftRight", borderStyle="thin",wrapText=TRUE),rows=rr,cols=1,stack=FALSE)
      addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",createStyle(textDecoration="bold",border="TopBottomLeftRight", borderStyle="thin",halign="center"),rows=rr,cols=2:7,gridExpand=TRUE,stack=FALSE)
      setRowHeights(.wbr,"Body Table 3 \u2014 PSA Frontier",rows=rr,heights=18)
    }
    fn3 <- length(FRONTIER) + 5
    writeData(.wbr,"Body Table 3 \u2014 PSA Frontier",
      "P(CE) = probability that the strategy has positive net monetary benefit vs passive case finding across 1,000 PSA simulations (CEAF approach). This is not the probability of being the optimal strategy: Symptom screen+CXR is cost-effective vs passive in 99.4% of simulations but is never the NMB-maximising strategy at NICE WTP thresholds (see CEAC, Figure 2).",
      startRow=fn3,startCol=1,colNames=FALSE)
    addStyle(.wbr,"Body Table 3 \u2014 PSA Frontier",.s_note(),rows=fn3,cols=1)
    mergeCells(.wbr,"Body Table 3 \u2014 PSA Frontier",cols=1:7,rows=fn3)
    setRowHeights(.wbr,"Body Table 3 \u2014 PSA Frontier",rows=fn3,heights=60)
  }

  saveWorkbook(.wbr,paste0("output/body/ResultsTable_paper_",.today,".xlsx"),overwrite=TRUE)
  cat("Saved: output/body/ResultsTable_paper_",.today,".xlsx\n",sep="")
}

# =============================================================================
# output/supplementary/Model_Parameters_supplementary.xlsx
# Full parameter table from input/config.csv — all sections, PSA ranges
# =============================================================================
cat("Building output/supplementary/Model_Parameters_supplementary.xlsx...\n")

.cfg <- read.csv("input/config.csv", stringsAsFactors = FALSE)
.cfg <- setNames(split(.cfg, seq_len(nrow(.cfg))), .cfg$parameter)
.get <- function(name, field) { r <- .cfg[[name]]; if (is.null(r)) "" else { v <- r[[field]]; if (is.na(v) || is.null(v)) "" else v } }

.fmt_val <- function(v) {
  if (is.na(v) || v == "") return("\u2014")
  n <- suppressWarnings(as.numeric(v))
  if (is.na(n)) return(as.character(v))
  if (n == 0) return("0")
  if (abs(n) < 0.0001) return(sprintf("%.6f", n))
  if (abs(n) < 0.01)   return(sprintf("%.5f", n))
  if (abs(n) < 1)      return(sprintf("%.4f", n))
  if (n == as.integer(n)) return(formatC(n, format = "d", big.mark = ","))
  return(sprintf("%.4f", n))
}

.psa_range <- function(name) {
  dist <- .get(name, "distribution")
  p1   <- suppressWarnings(as.numeric(.get(name, "dist_param1")))
  p2   <- suppressWarnings(as.numeric(.get(name, "dist_param2")))
  if (dist == "fixed" || dist == "" || is.na(p1) || is.na(p2)) return("Fixed")
  tryCatch({
    if (dist == "beta") {
      lo <- qbeta(0.025, p1, p2); hi <- qbeta(0.975, p1, p2)
    } else if (dist == "gamma") {
      lo <- qgamma(0.025, shape = p1, scale = p2); hi <- qgamma(0.975, shape = p1, scale = p2)
    } else return("\u2014")
    paste0(.fmt_val(lo), " \u2013 ", .fmt_val(hi))
  }, error = function(e) "\u2014")
}

.SECTIONS <- list(
  list("Structural",                                  c("n_t","n_c","discount")),
  list("Initial conditions",                          c("init_Uninfected","init_LatentUndiagnosed","init_ActiveUndiagnosed")),
  list("Monthly health state costs (\u00a3)",         c("Cost_Uninfected","Cost_LatentUndiagnosed","Cost_LatentDiagnosed","Cost_LatentTreated","Cost_LatentNotreated","Cost_LatentCompleted","Cost_LatentDiscontinued","Cost_LatentLtfu","Cost_ActiveUndiagnosed","Cost_ActiveDiagnosed","Cost_ActiveTreated","Cost_ActiveNotreated","Cost_ActiveCompleted","Cost_ActiveDiscontinued","Cost_ActiveLtfu","Cost_Dead")),
  list("Scenario costs \u2014 not used in base case", c("Cost_ActiveUndiagnosed_scenario","Cost_ActiveNotreated_scenario")),
  list("Health utilities (QALYs)",                    c("qaly_Uninfected","qaly_LatentUndiagnosed","qaly_LatentDiagnosed","qaly_LatentTreated","qaly_LatentNotreated","qaly_LatentCompleted","qaly_LatentDiscontinued","qaly_LatentLtfu","qaly_ActiveUndiagnosed","qaly_ActiveDiagnosed","qaly_ActiveTreated","qaly_ActiveNotreated","qaly_ActiveCompleted","qaly_ActiveDiscontinued","qaly_ActiveLtfu","qaly_Dead")),
  list("Screening test parameters",                   c("igra_spec_qft","igra_spec_tspt","igra_ltbi_sens_qft","igra_ltbi_sens_tspt","prev_ltbi_ukhsa","prev_active_low","prev_active_high")),
  list("Programme delivery costs (\u00a3 per 100,000 cohort)", c("cost_igra_programme_cough","cost_igra_programme_allsx")),
  list("Transmission scenario parameters",            c("p_tx_per_contact","n_contacts_base","n_contacts_dsa_mid","n_contacts_dsa_high","cost_secondary_tb")),
  list("Background mortality (ONS life tables, per 100,000 per year)", c("ons_mort_25_29","ons_mort_30_34","ons_mort_35_39","ons_mort_40_44","ons_mort_45_49","ons_mort_50_54","ons_mort_55_59","ons_mort_60_64","ons_mort_65_69","ons_mort_70_74","ons_mort_75_79")),
  list("Transition probabilities \u2014 disease progression", c("p_Uninfected_LatentUndiagnosed","p_LatentUndiagnosed_ActiveUndiagnosed","p_react_phase2","p_LatentDiagnosed_ActiveUndiagnosed","p_LatentNotreated_ActiveUndiagnosed","p_LatentDiscontinued_ActiveUndiagnosed","p_LatentLtfu_ActiveUndiagnosed")),
  list("Transition probabilities \u2014 LTBI care cascade", c("p_LatentUndiagnosed_LatentDiagnosed","p_LatentDiagnosed_LatentTreated","p_LatentDiagnosed_LatentNotreated","p_LatentTreated_LatentCompleted","p_LatentTreated_LatentDiscontinued","p_LatentTreated_LatentLtfu","p_LatentCompleted_Uninfected","p_LatentDiscontinued_LatentDiagnosed","p_LatentLtfu_LatentDiagnosed")),
  list("Transition probabilities \u2014 active TB care cascade", c("p_ActiveUndiagnosed_ActiveDiagnosed","p_ActiveDiagnosed_ActiveTreated","p_ActiveDiagnosed_ActiveNotreated","p_ActiveTreated_ActiveCompleted","p_ActiveTreated_ActiveDiscontinued","p_ActiveTreated_ActiveLtfu","p_ActiveCompleted_Uninfected","p_ActiveDiscontinued_ActiveUndiagnosed","p_ActiveLtfu_ActiveUndiagnosed","p_ActiveNotreated_ActiveDiagnosed")),
  list("Transition probabilities \u2014 mortality", c("p_Uninfected_Dead","p_LatentUndiagnosed_Dead","p_LatentDiagnosed_Dead","p_LatentTreated_Dead","p_LatentNotreated_Dead","p_LatentCompleted_Dead","p_LatentDiscontinued_Dead","p_LatentLtfu_Dead","p_ActiveUndiagnosed_Dead","p_ActiveDiagnosed_Dead","p_ActiveTreated_Dead","p_ActiveNotreated_Dead","p_ActiveCompleted_Dead","p_ActiveDiscontinued_Dead","p_ActiveLtfu_Dead"))
)

.NAVY2  <- "#1F4E79"; .SECBG <- "#D6E4F0"; .ALTBG <- "#EBF3FB"
.wbk    <- createWorkbook(); modifyBaseFont(.wbk, fontSize = 9, fontName = "Calibri")
addWorksheet(.wbk, "Parameters")
setColWidths(.wbk, "Parameters", cols = 1:7, widths = c(4, 42, 9, 9, 12, 22, 42))

# Title rows
writeData(.wbk, "Parameters", "TB Screening CEA \u2014 Full Model Parameter Table", startRow = 1, startCol = 1, colNames = FALSE)
addStyle(.wbk, "Parameters", createStyle(fontSize = 13, fontColour = .NAVY2, textDecoration = "bold"), rows = 1, cols = 1)
mergeCells(.wbk, "Parameters", cols = 1:7, rows = 1); setRowHeights(.wbk, "Parameters", rows = 1, heights = 18)

writeData(.wbk, "Parameters", paste0("55-year lifetime horizon \u00b7 n=100,000 \u00b7 3.5% discount \u00b7 2023/24 GBP \u00b7 ", .today), startRow = 2, startCol = 1, colNames = FALSE)
addStyle(.wbk, "Parameters", createStyle(fontSize = 9, fontColour = "#404040", textDecoration = "italic"), rows = 2, cols = 1)
mergeCells(.wbk, "Parameters", cols = 1:7, rows = 2); setRowHeights(.wbk, "Parameters", rows = 2, heights = 14)

# Header row
.khdrs <- c("#", "Parameter", "Value", "SD", "Distribution", "PSA range (2.5\u201397.5%)", "Source")
writeData(.wbk, "Parameters", as.data.frame(t(.khdrs)), startRow = 4, startCol = 1, colNames = FALSE)
addStyle(.wbk, "Parameters", createStyle(fgFill = .NAVY2, fontColour = "white", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderStyle = "thin"),
  rows = 4, cols = 1:7, gridExpand = TRUE)
setRowHeights(.wbk, "Parameters", rows = 4, heights = 28)

.krow <- 5L; .kn <- 0L
for (.sec in .SECTIONS) {
  .stitle <- .sec[[1]]; .sparams <- .sec[[2]]
  writeData(.wbk, "Parameters", paste0("  ", .stitle), startRow = .krow, startCol = 1, colNames = FALSE)
  addStyle(.wbk, "Parameters", createStyle(fgFill = .SECBG, fontColour = .NAVY2, textDecoration = "bold",
    fontSize = 9, border = "TopBottomLeftRight", borderStyle = "thin"),
    rows = .krow, cols = 1)
  mergeCells(.wbk, "Parameters", cols = 1:7, rows = .krow)
  setRowHeights(.wbk, "Parameters", rows = .krow, heights = 16)
  .krow <- .krow + 1L

  for (.pname in .sparams) {
    if (!(.pname %in% names(.cfg))) next
    .kn <- .kn + 1L
    .shade <- (.kn %% 2L == 0L)
    .bg <- if (.shade) .ALTBG else "white"
    .sf  <- createStyle(fgFill = .bg, border = "TopBottomLeftRight", borderStyle = "thin", wrapText = TRUE, valign = "center")
    .nf  <- createStyle(fgFill = .bg, border = "TopBottomLeftRight", borderStyle = "thin", halign = "center", valign = "center")

    .desc <- .get(.pname, "parameter_name"); if (.desc == "") .desc <- .pname
    .src  <- .get(.pname, "source_name"); if (nchar(.src) > 150) .src <- paste0(substr(.src, 1, 147), "...")
    .rv   <- c(.kn, .desc, .fmt_val(.get(.pname,"value")), .fmt_val(.get(.pname,"sd")),
               paste0(toupper(substr(.get(.pname,"distribution"),1,1)), substr(.get(.pname,"distribution"),2,99)),
               .psa_range(.pname), .src)
    writeData(.wbk, "Parameters", as.data.frame(t(.rv)), startRow = .krow, startCol = 1, colNames = FALSE)
    addStyle(.wbk, "Parameters", .nf, rows = .krow, cols = c(1,3,4,5,6), gridExpand = TRUE)
    addStyle(.wbk, "Parameters", .sf, rows = .krow, cols = c(2,7), gridExpand = TRUE)
    setRowHeights(.wbk, "Parameters", rows = .krow, heights = 30)
    .krow <- .krow + 1L
  }
}

# Footer
.krow <- .krow + 1L
writeData(.wbk, "Parameters",
  "Notes: PSA range = 2.5th\u201397.5th percentile of the assigned distribution. Fixed = not varied in PSA. SD = standard deviation used to fit PSA distribution. Disease progression uses two-phase reactivation: phase 1 (months 1\u201360) and phase 2 (months 61\u2013660) applied to Latent Undiagnosed only. Background mortality: age-varying ONS all-cause mortality (5-year bands, 25\u201379 years), updated every 60 cycles. Scenario costs not used in base case \u2014 sensitivity analysis only.",
  startRow = .krow, startCol = 1, colNames = FALSE)
addStyle(.wbk, "Parameters", createStyle(fontSize = 9, fontColour = "#404040", textDecoration = "italic", wrapText = TRUE, valign = "top"),
  rows = .krow, cols = 1)
mergeCells(.wbk, "Parameters", cols = 1:7, rows = .krow)
setRowHeights(.wbk, "Parameters", rows = .krow, heights = 48)

saveWorkbook(.wbk, "output/supplementary/Model_Parameters_supplementary.xlsx", overwrite = TRUE)
cat(sprintf("Saved: output/supplementary/Model_Parameters_supplementary.xlsx  (%d parameters)\n", .kn))

# =============================================================================
# output/supplementary/CHEERS_supplementary.xlsx — CHEERS 2022 compliance checklist
# =============================================================================
cat("Building output/supplementary/CHEERS_supplementary.xlsx...\n")

.cheers <- data.frame(stringsAsFactors = FALSE,
  section = c("Title","Abstract","Introduction","Introduction",
    rep("Methods", 16), rep("Results", 4), "Discussion",
    rep("Other", 3)),
  item = 1:28,
  label = c("Title","Abstract","Background and objectives","Health economic analysis plan",
    "Study population","Setting and location","Comparators","Perspective","Time horizon",
    "Discount rate","Selection of outcomes","Measurement of outcomes","Valuation of outcomes",
    "Measurement and valuation of resources and costs","Currency, price date, and conversion",
    "Rationale and description of model","Analytics and assumptions",
    "Characterising heterogeneity","Characterising distributional effects","Characterising uncertainty",
    "Study parameters","Summary of main results","Effect of uncertainty",
    "Effect of engagement with patients and others affected",
    "Study findings, limitations, generalisability, and current knowledge",
    "Source of funding","Conflicts of interest","Availability of data, code, and other materials"),
  where = c("Title","Abstract","Introduction","N/A",
    "Methods \u2014 Model structure","Methods \u2014 Setting","Methods \u2014 Strategies; Table 1",
    "Methods \u2014 Perspective; Table 1","Methods \u2014 Model structure","Methods \u2014 Model structure",
    "Methods \u2014 Outcomes","Methods \u2014 QALYs; ST2","Methods \u2014 QALYs; ST2",
    "Methods \u2014 Costs; ST2","Methods \u2014 Costs","Methods \u2014 Model structure; SF1; Supplement",
    "Methods \u2014 Model structure; Supplement","Sensitivity analyses","Discussion \u2014 Limitations",
    "Methods \u2014 Uncertainty; Results \u2014 SA","ST2 (Model Parameters)","Results \u2014 Table 2; Table 3; Figure 1",
    "Results \u2014 Sensitivity analyses; Figures 3\u20134; supplementary figures",
    "N/A \u2014 modelling study","Discussion","Acknowledgements / Funding","Declarations","Data sharing statement"),
  status = c("Y","Y","Y","N","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","P","P","Y","Y","Y","Y","N","Y","P","P","Y"),
  notes = c(
    "Title states 'cost-effectiveness analysis' and lists the 43 screening strategies vs passive case finding.",
    "5-part structured abstract (Background/Methods/Findings/Interpretation/Funding) \u2264250 words.",
    "TB burden in UK migrants, current NICE guidance, policy gap, study objective stated.",
    "No pre-registered analysis plan. Pre-registration not standard for modelling CEA at time of study.",
    "UK migrants from high TB-burden countries; entry age 30; 100,000 cohort; LTBI prevalence 17.8%.",
    "England, NHS/UK migrant screening context; 2023/24 GBP; NICE reference case.",
    "43 strategies (passive case finding + 34 sequential CXR-based + 8 parallel IGRA) based on Zenner 2025 decision tree.",
    "Primary: societal (NHS costs \u2212 secondary transmission savings). Robustness: NHS perspective (standard NICE reference case).",
    "55-year lifetime horizon (660 monthly cycles); NICE reference case for chronic disease.",
    "3.5% per annum for costs and outcomes; NICE reference case (NICE 2022 Methods Guide).",
    "QALYs (primary benefit); active TB cases prevented, secondary cases avoided (secondary outcomes).",
    "EQ-5D utility weights from UK general population studies (Kind 1998; Dolan 1997; Jit 2011). State-specific weights applied per Markov cycle.",
    "EQ-5D TTO values; UK general population tariff. See ST2 for all utility weights and sources.",
    "NHS reference costs (NHSBC 2023/24); NICE CG117 for active TB treatment; UKHSA 2024 for LTBI; unit costs inflated to 2023/24 GBP using NHS CPI.",
    "2023/24 GBP. Older cost data inflated using NHS CPI index. No currency conversion required (UK study).",
    "Markov cohort model, 16 health states, 660 monthly cycles. Chosen for tractability over individual-level simulation given 43 strategies and 1,000 PSA runs. Model structure diagram: SF1.",
    "Matrix multiplication per cycle; cohort of 100,000; state membership tracked over 55 years. PSA: beta distributions for probabilities/utilities, gamma for costs (n=1,000).",
    "Active TB prevalence SA covers heterogeneous populations (pooled migrants vs refugees). Age-varying background mortality included. Sex disaggregation not modelled (see Limitations).",
    "Not formally modelled. Equity implications discussed qualitatively in Discussion. SAGER sex-disaggregation limitation noted.",
    "PSA (n=1,000; beta/gamma distributions); DSA (one-way, top parameters by NMB impact); scenario SA (6 parameters); transmission DSA (5\u201315 contacts). CEAC at NICE WTP range.",
    "All 88 parameters with point estimates, PSA distributions, and references in ST2 (Model_Parameters_supplementary.xlsx).",
    "CE scatter (Fig 1), ICER table (Table 2/3), P(CE) (Table 3/Fig 2), frontier strategies with costs, QALYs, ICERs reported.",
    "PSA CEAC (Fig 2); tornado DSA (Fig 3); 6 one-way SAs; transmission DSA; results robust across all SAs.",
    "Decision-analytic modelling study. No patient engagement in parameter selection (standard for this study type). Limitation noted.",
    "Key finding: Parallel Sx+QFT cost-effective at NICE thresholds (societal perspective). Limitations: no sex disaggregation, behavioural assumptions, passive case finding calibration.",
    "TBC \u2014 confirm funding source with Kasim before submission.",
    "TBC \u2014 COI statements for all authors to be completed before submission.",
    "R code (MasterTBModel.R) and config.csv available as supplementary files. Parameter values in Model_Parameters_supplementary.xlsx."
  )
)

.STATUS_BG   <- c("Y" = "#E8F5E9", "P" = "#FFF8E1", "N" = "#F5F5F5")
.STATUS_FC   <- c("Y" = "#005000", "P" = "#7B4800", "N" = "#666666")
.STATUS_LABEL <- c("Y" = "Reported", "P" = "Partial", "N" = "N/A")

.wbc <- createWorkbook(); modifyBaseFont(.wbc, fontSize = 10, fontName = "Calibri")
addWorksheet(.wbc, "CHEERS 2022")
setColWidths(.wbc, "CHEERS 2022", cols = 1:7, widths = c(12, 5, 30, 60, 25, 10, 60))

# Title
writeData(.wbc, "CHEERS 2022", "ST1. CHEERS 2022 Checklist \u2014 Cost-Effectiveness of TB Screening Strategies for UK Migrants",
          startRow = 1, startCol = 1, colNames = FALSE)
addStyle(.wbc, "CHEERS 2022", createStyle(fgFill = "#00468B", fontColour = "white", textDecoration = "bold",
  fontSize = 13, halign = "center", valign = "center"), rows = 1, cols = 1)
mergeCells(.wbc, "CHEERS 2022", cols = 1:7, rows = 1); setRowHeights(.wbc, "CHEERS 2022", rows = 1, heights = 28)

writeData(.wbc, "CHEERS 2022",
  "Reporting standard: Husereau et al. (2022) CHEERS 2022 \u2014 Value Health 25(1):1-9. Status: Y=fully reported | P=partial/in progress | N=not applicable.",
  startRow = 2, startCol = 1, colNames = FALSE)
addStyle(.wbc, "CHEERS 2022", createStyle(fontSize = 9, fontColour = "#444444", textDecoration = "italic", wrapText = TRUE, valign = "top"),
  rows = 2, cols = 1)
mergeCells(.wbc, "CHEERS 2022", cols = 1:7, rows = 2); setRowHeights(.wbc, "CHEERS 2022", rows = 2, heights = 36)

# Column headers
writeData(.wbc, "CHEERS 2022", as.data.frame(t(c("Section","Item","Label","Description","Where addressed","Status","Notes"))),
          startRow = 3, startCol = 1, colNames = FALSE)
addStyle(.wbc, "CHEERS 2022", createStyle(fgFill = "#00468B", fontColour = "white", textDecoration = "bold",
  fontSize = 10, halign = "center", valign = "center", wrapText = TRUE,
  border = "TopBottomLeftRight", borderStyle = "thin"), rows = 3, cols = 1:7, gridExpand = TRUE)
setRowHeights(.wbc, "CHEERS 2022", rows = 3, heights = 22)

# Data rows — section merging tracked
.prev_sec <- ""; .sec_start <- 4L; .sec_ranges <- list()
for (.i in seq_len(nrow(.cheers))) {
  .rr <- .i + 3L; .st <- .cheers$status[.i]
  .bg <- .STATUS_BG[.st]
  .sf_c <- createStyle(fgFill = .bg, border = "TopBottomLeftRight", borderStyle = "thin", wrapText = TRUE, valign = "top")
  .nf_c <- createStyle(fgFill = .bg, border = "TopBottomLeftRight", borderStyle = "thin", halign = "center", valign = "top")
  .st_c <- createStyle(fgFill = .bg, border = "TopBottomLeftRight", borderStyle = "thin", halign = "center", valign = "top",
    textDecoration = "bold", fontColour = .STATUS_FC[.st])

  if (.cheers$section[.i] != .prev_sec) {
    if (.prev_sec != "") .sec_ranges[[length(.sec_ranges) + 1]] <- c(.prev_sec, .sec_start, .rr - 1L)
    .sec_start <- .rr; .prev_sec <- .cheers$section[.i]
  }
  .rv <- list(.cheers$section[.i], .cheers$item[.i], .cheers$label[.i],
              "", .cheers$where[.i], .STATUS_LABEL[.st], .cheers$notes[.i])
  for (.ci in seq_along(.rv)) {
    writeData(.wbc, "CHEERS 2022", .rv[[.ci]], startRow = .rr, startCol = .ci, colNames = FALSE)
    addStyle(.wbc, "CHEERS 2022",
      if (.ci %in% c(1,2,6)) .nf_c else if (.ci == 6) .st_c else .sf_c,
      rows = .rr, cols = .ci)
  }
  # Re-apply status style to col 6
  addStyle(.wbc, "CHEERS 2022", .st_c, rows = .rr, cols = 6)
  # Description col 4: write separately (too long for list)
  writeData(.wbc, "CHEERS 2022",
    c("Identify the study as an economic evaluation and specify the interventions being compared.",
      "Provide a structured summary that highlights context, key methods, results, and alternative analyses.",
      "Give the context for the study, the study question and its relevance for health policy or practice decisions.",
      "Indicate whether a health economic analysis plan was developed and where it can be accessed.",
      "Describe characteristics of the study population with reference to the decision problem.",
      "Provide relevant context to the decision problem.",
      "Describe the interventions or strategies being compared and why chosen.",
      "State the perspective(s) adopted by the study and why chosen.",
      "State the time horizon for the study and why appropriate.",
      "Report the choice of discount rate(s) used and why appropriate.",
      "Describe what outcomes were used as the measure(s) of benefit in the evaluation and their relevance.",
      "Describe the population and methods used to measure and value outcomes.",
      "Describe the population and methods used to value outcomes.",
      "Describe how costs were valued.",
      "Report the dates of the cost data used, the currency, and methods for converting costs.",
      "Describe and give reasons for the type of decision-analytic model used, key input parameters, assumptions, and limitations.",
      "Describe any methods for transforming inputs into outputs.",
      "If applicable, describe any methods for estimating how the results vary for subgroups.",
      "If applicable, describe how impacts are distributed across different individuals.",
      "Describe the methods used to characterise uncertainty in the analysis.",
      "Report the values, ranges, references, and probability distributions for all parameters.",
      "Report the mean values for the main categories of costs and outcomes of interest.",
      "Describe how uncertainty about analytic judgements, inputs, or projections affect findings.",
      "Report on any engagement with patients or members of the public and their effects on the analysis.",
      "Report key findings, limitations, ethical or equity considerations, and how these affect practice/policy.",
      "Describe how the study was funded and any role the funder played in the study.",
      "Report authors' conflicts of interest according to journal or professional requirements.",
      "Report how data, analytic code, and other materials used can be accessed."
    )[.i],
    startRow = .rr, startCol = 4, colNames = FALSE)
  addStyle(.wbc, "CHEERS 2022", .sf_c, rows = .rr, cols = 4)
  setRowHeights(.wbc, "CHEERS 2022", rows = .rr, heights = 40)
}
.sec_ranges[[length(.sec_ranges) + 1]] <- c(.prev_sec, .sec_start, nrow(.cheers) + 3L)

# Merge and style section column
for (.sr in .sec_ranges) {
  .r1 <- as.integer(.sr[2]); .r2 <- as.integer(.sr[3])
  if (.r1 < .r2) mergeCells(.wbc, "CHEERS 2022", cols = 1, rows = .r1:.r2)
  addStyle(.wbc, "CHEERS 2022", createStyle(fgFill = "#00468B", fontColour = "white", textDecoration = "bold",
    fontSize = 10, halign = "center", valign = "center", wrapText = TRUE,
    border = "TopBottomLeftRight", borderStyle = "thin"), rows = .r1, cols = 1)
}

# Summary row
.n_y <- sum(.cheers$status == "Y"); .n_p <- sum(.cheers$status == "P"); .n_n <- sum(.cheers$status == "N")
.sum_row <- nrow(.cheers) + 5L
writeData(.wbc, "CHEERS 2022",
  sprintf("Summary: %d/28 items fully reported (Y) | %d/28 partial (P) | %d/28 not applicable (N). Items 26\u201327 (funding/COI) to be completed before submission.",
          .n_y, .n_p, .n_n),
  startRow = .sum_row, startCol = 1, colNames = FALSE)
addStyle(.wbc, "CHEERS 2022", createStyle(textDecoration = "bold", fontSize = 10, wrapText = TRUE, valign = "top"),
  rows = .sum_row, cols = 1)
mergeCells(.wbc, "CHEERS 2022", cols = 1:7, rows = .sum_row)
setRowHeights(.wbc, "CHEERS 2022", rows = .sum_row, heights = 28)

saveWorkbook(.wbc, "output/supplementary/CHEERS_supplementary.xlsx", overwrite = TRUE)
cat("Saved: output/supplementary/CHEERS_supplementary.xlsx\n")

cat("\nAll figures and tables complete.\n")
