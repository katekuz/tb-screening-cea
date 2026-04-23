# TB Screening CEA — Markov cohort model (16 states, 43 strategies, lifetime horizon)
# MSc Modelling for Global Health, University of Oxford

library(ggplot2)    # plotting
library(patchwork)  # multi-panel figure assembly (must load before ggplot2 + overrides)
library(tidyr)      # data reshaping
library(dplyr)      # data manipulation
library(tibble)     # modern data frames
library(readxl)     # reading Excel initial conditions
library(ggrepel)    # non-overlapping text labels on plots
library(parallel)   # mclapply — PSA parallelisation
library(ggsci)      # Lancet colour palette

# Suppress scientific notation globally (never use e.g. 1e+05 in plots or output)
options(scipen = 999)

# =============================================================================
# SHARED FIGURE PALETTE — single source of truth for all plots
# pal_lancet("lanonc", alpha = 0.85): slightly muted Lancet colours.
# Index: 1=navy, 2=red, 3=green, 4=teal, 5=purple, 6=peach, 7=dark-red, 8=grey
# =============================================================================
LANCET <- pal_lancet("lanonc", alpha = 0.85)(8)
STRAT_COLS <- c(
  "Passive case finding"       = LANCET[8],   # grey   — reference / no screening
  "Cough+CXR (TB sx)"          = LANCET[4],   # teal   — frontier
  "Symptom screen+CXR"         = LANCET[3],   # green  — frontier
  "Parallel Sx+QFT (Ultra)"    = LANCET[2],   # red    — frontier
  "Parallel Cough+QFT (Ultra)" = LANCET[5]    # purple — extendedly dominated
)
FRONTIER_STRATS <- names(STRAT_COLS)[1:4]


source("scripts/model_functions.R")
# =============================================================================
# BASE CASE ANALYSIS
# Runs the Markov model using central (mean) parameter values from config.csv
# for the no-screening scenario. This produces the reference trajectory
# against which all screening strategies are compared.
# =============================================================================

results <- model(paramsData)

# =============================================================================
# BASE CASE RESULTS AND VISUALISATION
# Conservative estimate — NHS perspective only.
# The primary analysis (transmission prevention, societal perspective) follows.
# =============================================================================

library(scales)

# -------------------- Summary Output -------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                    TB SCREENING MODEL - RESULTS SUMMARY                        \n")
cat("================================================================================\n")
cat("\n")
cat("INITIAL STATE DISTRIBUTION:\n")
cat(sprintf("  Uninfected:          %s (%.1f%%)\n",
    format(results$state_membership[1, "Uninfected"], big.mark=","),
    100 * results$state_membership[1, "Uninfected"] / n_c))
cat(sprintf("  Latent Undiagnosed:  %s (%.1f%%)\n",
    format(results$state_membership[1, "LatentUndiagnosed"], big.mark=","),
    100 * results$state_membership[1, "LatentUndiagnosed"] / n_c))
cat(sprintf("  Active Undiagnosed:  %s (%.1f%%)\n",
    format(results$state_membership[1, "ActiveUndiagnosed"], big.mark=","),
    100 * results$state_membership[1, "ActiveUndiagnosed"] / n_c))
cat("\n")
cat(sprintf("OUTCOMES OVER %d YEARS (with half-cycle correction):\n", n_t %/% 12))
cat(sprintf("  Total discounted costs:  £%s\n", format(round(results$total_cost), big.mark=",")))
cat(sprintf("  Total discounted QALYs:  %s\n", format(round(results$total_qaly, 1), big.mark=",")))
cat(sprintf("  Deaths:                  %s\n", format(round(results$state_membership[n_t, "Dead"]), big.mark=",")))
cat("\n")
cat("================================================================================\n")
cat("\n")

# -------------------- Custom Color Palette -------------------------------------
state_colors <- c(
  # Uninfected — light lavender
  "Uninfected"         = "#e8eef8",
  # Latent TB — purple family (light to dark)
  "LatentUndiagnosed"  = "#bccfe8",
  "LatentDiagnosed"    = "#80b4d8",
  "LatentTreated"      = "#0099b4",
  "LatentNotreated"    = "#4a90c4",
  "LatentCompleted"    = "#d0e8f4",
  "LatentDiscontinued" = "#2070a8",
  "LatentLtfu"         = "#6aaac8",
  # Active TB — rose/pink family (light to dark)
  "ActiveUndiagnosed"  = "#fdd0b8",
  "ActiveDiagnosed"    = "#ed6060",
  "ActiveTreated"      = "#ad002a",
  "ActiveNotreated"    = "#c83040",
  "ActiveCompleted"    = "#ffe8d8",
  "ActiveDiscontinued" = "#9b1020",
  "ActiveLtfu"         = "#d85040",
  # Dead
  "Dead"               = "#2d2d2d"
)

# -------------------- State Trajectories Plot ----------------------------------
state_mat <- as.matrix(results$state_membership)
state_df <- as.data.frame(state_mat)
state_df$month <- 1:nrow(state_df)
long <- pivot_longer(state_df, cols = -month, names_to = "state", values_to = "n")

state_level_order <- c(
  "Uninfected",
  "LatentUndiagnosed", "LatentDiagnosed", "LatentTreated", "LatentCompleted",
  "LatentNotreated", "LatentDiscontinued", "LatentLtfu",
  "ActiveUndiagnosed", "ActiveDiagnosed", "ActiveTreated", "ActiveCompleted",
  "ActiveNotreated", "ActiveDiscontinued", "ActiveLtfu",
  "Dead"
)
long$state <- factor(long$state, levels = state_level_order)

state_labels <- c(
  "Uninfected"         = "Uninfected",
  "LatentUndiagnosed"  = "LTBI: undiagnosed",
  "LatentDiagnosed"    = "LTBI: diagnosed",
  "LatentTreated"      = "LTBI: on treatment",
  "LatentCompleted"    = "LTBI: treatment completed",
  "LatentNotreated"    = "LTBI: not treated",
  "LatentDiscontinued" = "LTBI: treatment discontinued",
  "LatentLtfu"         = "LTBI: lost to follow-up",
  "ActiveUndiagnosed"  = "Active TB: undiagnosed",
  "ActiveDiagnosed"    = "Active TB: diagnosed",
  "ActiveTreated"      = "Active TB: on treatment",
  "ActiveCompleted"    = "Active TB: treatment completed",
  "ActiveNotreated"    = "Active TB: not treated",
  "ActiveDiscontinued" = "Active TB: treatment discontinued",
  "ActiveLtfu"         = "Active TB: lost to follow-up",
  "Dead"               = "Dead"
)

p_states <- ggplot(long %>% filter(n > 0), aes(x = month, y = n, color = state)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = state_colors,
                     breaks = names(state_labels),
                     labels = state_labels) +
  scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_x_continuous(breaks = seq(0, n_t, by = 60), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Time (years)",
    y = "Number of individuals (log\u2081\u2080 scale)",
    color = "Health state"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10)
  )

# -------------------- Display & Save State Trajectories -----------------------
cat("Plots saved to output/ folder\n")

# ================================================================================
# EPIDEMIOLOGY & IN-DEPTH ANALYSIS PLOTS
# ================================================================================

# -------------------- 1. TB Incidence Plot -------------------------------------
# Calculates monthly active TB incidence by counting new transitions from
# latent states into ActiveUndiagnosed. This reflects genuine disease
# progression from latent to active TB rather than the snapshot state counts.
# Re-entry from ActiveDiscontinued/ActiveLtfu is excluded as these individuals
# already had active TB and are not new incident cases.
# Sources of new active cases per cycle:
#   LatentUndiagnosed -> ActiveUndiagnosed (observed incidence; Berrocal-Almanza 2022)
#   LatentDiagnosed   -> ActiveUndiagnosed (progression despite diagnosis)
#   LatentNotreated   -> ActiveUndiagnosed (untreated LTBI progression)
#   LatentDiscontinued -> ActiveUndiagnosed (partial treatment failure)
#   LatentLtfu        -> ActiveUndiagnosed (lost to follow-up; no treatment protection)
active_states <- c("ActiveUndiagnosed", "ActiveDiagnosed", "ActiveTreated",
                   "ActiveNotreated", "ActiveCompleted", "ActiveDiscontinued", "ActiveLtfu")

m_p_base <- results$m_p
new_active_cases <- numeric(n_t)
new_active_cases[1] <- 0
for (t in 2:n_t) {
  # People transitioning INTO ActiveUndiagnosed from latent/other states
  new_active_cases[t] <-
    state_mat[t-1, "LatentUndiagnosed"] * m_p_base["LatentUndiagnosed", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentDiagnosed"] * m_p_base["LatentDiagnosed", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentNotreated"] * m_p_base["LatentNotreated", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentDiscontinued"] * m_p_base["LatentDiscontinued", "ActiveUndiagnosed"] +
    state_mat[t-1, "LatentLtfu"] * m_p_base["LatentLtfu", "ActiveUndiagnosed"]
}

person_months <- rowSums(state_mat[, v_state_names[v_state_names != "Dead"]])
incidence_per_1000_py <- (new_active_cases / person_months) * 1000 * 12

df_incidence <- tibble(month = 1:n_t, incidence = incidence_per_1000_py)

p_incidence <- ggplot(df_incidence, aes(x = month, y = incidence)) +
  geom_line(color = "#ed0000", linewidth = 1.2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#ed0000", fill = "#fdaf91") +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Incidence per 1,000 person-years",
       title = NULL,
       subtitle = "No Screening scenario | Raw per-cycle rate (line) + LOESS smooth (band)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank())

# -------------------- 2. TB Prevalence Plot ------------------------------------
latent_states <- c("LatentUndiagnosed", "LatentDiagnosed", "LatentTreated",
                   "LatentNotreated", "LatentCompleted", "LatentDiscontinued", "LatentLtfu")

df_prevalence <- tibble(
  month = rep(1:n_t, 3),
  category = factor(rep(c("Uninfected", "Latent TB", "Active TB"), each = n_t),
                    levels = c("Active TB", "Latent TB", "Uninfected")),
  count = c(state_mat[, "Uninfected"], rowSums(state_mat[, latent_states]), rowSums(state_mat[, active_states]))
)

p_prevalence <- ggplot(df_prevalence, aes(x = month, y = count, fill = category)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Active TB" = "#ed0000", "Latent TB" = "#80b4d8", "Uninfected" = "#e8eef8")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Number of individuals",
       title = NULL,
       subtitle = "Distribution of uninfected, latent, and active TB", fill = "Status") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 3. Treatment Cascade (Cumulative Flows) ------------------
# The treatment cascade quantifies how many individuals pass through each
# stage of the care pathway over 20 years. Cumulative flows are estimated
# by summing the transition probabilities applied to state membership at
# each cycle — this gives true individual counts rather than person-months.
# Initial state values (month 1) are added to capture those detected at entry.
m_p_base <- results$m_p

# For treatment cascade and cumulative outcomes: use Parallel Sx+QFT (Ultra)
# — the frontier parallel IGRA strategy (seq. ICER £11,324/QALY, NHS). This shows
# the full LTBI treatment cascade (11,526 LTBI detected via specificity-corrected QFT-GIT).
# Inline init: TP=874.6, FN=125.4 per 100k (decision tree; Parallel Sx+QFT Ultra);
# LTBI = n_ltbi × igra_sens_qft − FP (specificity-corrected; Pai 2008)
init_parallel_sx_qft <- setNames(rep(0, length(v_state_names)), v_state_names)
ltbi_det_pcqft <- max(0, n_c * prev_ltbi * igra_sens_qft - (1 - igra_spec_qft) * n_c * prev_uninfected)  # specificity-corrected: ~11,526
init_parallel_sx_qft["ActiveDiagnosed"]   <- 874.6 * (n_c / excel_denom)
init_parallel_sx_qft["ActiveUndiagnosed"] <- 125.4 * (n_c / excel_denom)
init_parallel_sx_qft["LatentDiagnosed"]   <- ltbi_det_pcqft
init_parallel_sx_qft["LatentUndiagnosed"] <- n_c * prev_ltbi - ltbi_det_pcqft
init_parallel_sx_qft["Uninfected"]        <- n_c * prev_uninfected
results_psqft   <- model(paramsData, init_dist = init_parallel_sx_qft)
state_mat_psqft <- as.matrix(results_psqft$state_membership)
state_mat       <- state_mat_psqft   # override: cascade + cumulative use Parallel Sx+QFT (Ultra)

# Person-month sums below are superseded by the transition-flow approach (cum_flow_*).
# Retained for validation against independent derivation.
cum_ltbi_diagnosed <- sum(state_mat[, "LatentDiagnosed"]) +
  sum(state_mat[, "LatentTreated"]) + sum(state_mat[, "LatentNotreated"]) +
  sum(state_mat[, "LatentCompleted"]) + sum(state_mat[, "LatentDiscontinued"]) +
  sum(state_mat[, "LatentLtfu"])
cum_ltbi_treated <- sum(state_mat[, "LatentTreated"]) +
  sum(state_mat[, "LatentCompleted"]) + sum(state_mat[, "LatentDiscontinued"]) +
  sum(state_mat[, "LatentLtfu"])
cum_ltbi_completed    <- sum(state_mat[, "LatentCompleted"])
cum_active_diagnosed  <- sum(state_mat[, "ActiveDiagnosed"]) +
  sum(state_mat[, "ActiveTreated"]) + sum(state_mat[, "ActiveNotreated"]) +
  sum(state_mat[, "ActiveCompleted"]) + sum(state_mat[, "ActiveDiscontinued"]) +
  sum(state_mat[, "ActiveLtfu"])
cum_active_treated    <- sum(state_mat[, "ActiveTreated"]) +
  sum(state_mat[, "ActiveCompleted"]) + sum(state_mat[, "ActiveDiscontinued"]) +
  sum(state_mat[, "ActiveLtfu"])
cum_active_completed  <- sum(state_mat[, "ActiveCompleted"])

# Transition-flow method: count individuals moving between cascade stages
# each month and accumulate over the full time horizon.
cum_flow_ltbi_diag <- 0
cum_flow_ltbi_treat <- 0
cum_flow_ltbi_comp <- 0
cum_flow_active_diag <- 0
cum_flow_active_treat <- 0
cum_flow_active_comp <- 0
for (t in 2:n_t) {
  cum_flow_ltbi_diag  <- cum_flow_ltbi_diag +
    state_mat[t-1, "LatentUndiagnosed"] * m_p_base["LatentUndiagnosed", "LatentDiagnosed"]
  cum_flow_ltbi_treat <- cum_flow_ltbi_treat +
    state_mat[t-1, "LatentDiagnosed"] * m_p_base["LatentDiagnosed", "LatentTreated"]
  cum_flow_ltbi_comp  <- cum_flow_ltbi_comp +
    state_mat[t-1, "LatentTreated"] * m_p_base["LatentTreated", "LatentCompleted"]
  cum_flow_active_diag  <- cum_flow_active_diag +
    state_mat[t-1, "ActiveUndiagnosed"] * m_p_base["ActiveUndiagnosed", "ActiveDiagnosed"]
  cum_flow_active_treat <- cum_flow_active_treat +
    state_mat[t-1, "ActiveDiagnosed"] * m_p_base["ActiveDiagnosed", "ActiveTreated"]
  cum_flow_active_comp  <- cum_flow_active_comp +
    state_mat[t-1, "ActiveTreated"] * m_p_base["ActiveTreated", "ActiveCompleted"]
}
# Add initial state counts (people who start already diagnosed/treated)
cum_flow_ltbi_diag  <- cum_flow_ltbi_diag + state_mat[1, "LatentDiagnosed"]
cum_flow_active_diag <- cum_flow_active_diag + state_mat[1, "ActiveDiagnosed"]

cascade_data <- tibble(
  stage = factor(c("LTBI Diagnosed", "LTBI Treated", "LTBI Completed",
                   "Active Diagnosed", "Active Treated", "Active Completed"),
                 levels = c("LTBI Diagnosed", "LTBI Treated", "LTBI Completed",
                           "Active Diagnosed", "Active Treated", "Active Completed")),
  count = c(cum_flow_ltbi_diag, cum_flow_ltbi_treat, cum_flow_ltbi_comp,
            cum_flow_active_diag, cum_flow_active_treat, cum_flow_active_comp),
  type = factor(c(rep("Latent TB", 3), rep("Active TB", 3)))
)

p_cascade <- ggplot(cascade_data, aes(x = stage, y = count, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(count)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Latent TB" = "#80b4d8", "Active TB" = "#ed0000")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Number of individuals",
       title = NULL,
       subtitle = sprintf("Parallel Sx+QFT (Ultra) strategy | Cumulative individuals passing through each care stage over %d years", n_t %/% 12), fill = "") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# -------------------- 4. Cumulative Outcomes -----------------------------------
df_cumulative <- tibble(
  month = rep(1:n_t, 4),
  outcome = factor(rep(c("Deaths", "Treatment Completed", "Lost to Follow-up", "Treatment Discontinued"), each = n_t),
                   levels = c("Deaths", "Treatment Completed", "Lost to Follow-up", "Treatment Discontinued")),
  count = c(state_mat[, "Dead"],
            state_mat[, "LatentCompleted"] + state_mat[, "ActiveCompleted"],
            state_mat[, "LatentLtfu"] + state_mat[, "ActiveLtfu"],
            state_mat[, "LatentDiscontinued"] + state_mat[, "ActiveDiscontinued"])
)

p_cumulative <- ggplot(df_cumulative, aes(x = month, y = count, color = outcome)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Deaths" = "#adb6b6", "Treatment Completed" = "#ed0000",
                                 "Lost to Follow-up" = "#80b4d8", "Treatment Discontinued" = "#fdaf91")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Cumulative count",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Deaths, treatment completions, LTFU, and discontinuations", color = "Outcome") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 5. Stacked State Distribution ----------------------------
state_categories <- tibble(
  state = v_state_names,
  category = case_when(
    state == "Uninfected" ~ "Uninfected",
    state == "Dead" ~ "Dead",
    grepl("^Latent", state) ~ "Latent TB",
    grepl("^Active", state) ~ "Active TB"
  )
)


df_stacked <- as.data.frame(state_mat) %>%
  mutate(month = 1:n_t) %>%
  pivot_longer(cols = -month, names_to = "state", values_to = "n") %>%
  left_join(state_categories, by = "state") %>%
  group_by(month, category) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Dead", "Active TB", "Latent TB", "Uninfected")))

p_stacked <- ggplot(df_stacked, aes(x = month, y = n, fill = category)) +
  geom_area(alpha = 0.9) +
  scale_fill_manual(values = c("Dead" = "#2d2d2d", "Active TB" = "#ed0000",
                                "Latent TB" = "#80b4d8", "Uninfected" = "#e8eef8")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, n_t, by = 24), labels = function(x) paste0(x/12, "y")) +
  theme_minimal(base_size = 14) +
  labs(x = "Time", y = "Number of individuals",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Stacked area chart showing disease progression", fill = "Category") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "right")

# -------------------- 6. Final State Pie Chart ---------------------------------
final_dist <- df_stacked %>%
  filter(month == n_t) %>%
  mutate(pct = n / sum(n) * 100, label = paste0(category, "\n", round(pct, 1), "%"))

final_dist <- final_dist %>%
  mutate(n_fmt = format(round(n), big.mark = ","))

# Pie chart: ggplot2 coord_polar stacks in REVERSE factor level order, so
# midpoints must be computed over data sorted in reverse factor level order.
pie_cols <- c("Dead" = "#ad002a", "Active TB" = "#ed0000",
              "Latent TB" = "#fdaf91", "Uninfected" = "#e8eef8")
final_dist_pie <- final_dist %>%
  mutate(category = factor(category,
           levels = c("Dead", "Active TB", "Latent TB", "Uninfected"))) %>%
  arrange(desc(category)) %>%           # reverse order for correct cumsum
  mutate(midpoint   = cumsum(n) - n / 2,
         label_text = paste0(category, "\n", n_fmt, "\n(", round(pct, 1), "%)"))
p_pie <- ggplot(final_dist_pie, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = pie_cols) +
  geom_label(aes(y = midpoint, label = label_text),
             x = 1.6, size = 2.8, show.legend = FALSE,
             label.padding = unit(0.15, "lines"), linewidth = 0.3) +
  theme_void(base_size = 14) +
  labs(title = paste0("Final Population Distribution (Year ", n_t %/% 12, ")"),
       subtitle = paste0("Parallel Sx+QFT (Ultra) strategy | Total population: ", prettyNum(n_c, big.mark = ",")),
       fill = "Category") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5), legend.position = "right")

# -------------------- 7. Cost Breakdown Bar Chart ------------------------------
cost_by_state <- tibble(
  state = v_state_names,
  total_cost = colSums(state_mat * rep(as.numeric(results$m_payoffsCost), each = n_t))
) %>%
  left_join(state_categories, by = "state") %>%
  group_by(category) %>%
  summarise(total_cost = sum(total_cost), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Active TB", "Latent TB", "Uninfected", "Dead")),
         pct = total_cost / sum(total_cost) * 100)

p_cost_breakdown <- ggplot(cost_by_state, aes(x = reorder(category, -total_cost), y = total_cost, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("£", format(round(total_cost), big.mark = ","), "\n(", round(pct, 1), "%)")),
            vjust = -0.2, size = 4) +
  scale_fill_manual(values = c("Dead" = "#2d2d2d", "Active TB" = "#ed0000",
                                "Latent TB" = "#80b4d8", "Uninfected" = "#e8eef8")) +
  scale_y_continuous(labels = label_comma(prefix = "£"), expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Total undiscounted costs",
       title = NULL,
       subtitle = "Parallel Sx+QFT (Ultra) strategy | Undiscounted | Which disease states drive costs?") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "none")

# -------------------- 8. Summary Statistics ------------------------------------
cat("\n")
cat("================================================================================\n")
cat("                    EPIDEMIOLOGY SUMMARY STATISTICS                             \n")
cat("================================================================================\n")
cat("\n")
cat("DISEASE BURDEN:\n")
cat(sprintf("  Peak LTBI prevalence:     %s (%.1f%% of cohort)\n",
    format(round(max(rowSums(state_mat[, latent_states]))), big.mark = ","),
    100 * max(rowSums(state_mat[, latent_states])) / n_c))
cat(sprintf("  Peak Active TB:           %s (%.2f%% of cohort)\n",
    format(round(max(rowSums(state_mat[, active_states]))), big.mark = ","),
    100 * max(rowSums(state_mat[, active_states])) / n_c))
cat(sprintf("  Total deaths:             %s (%.2f%% mortality)\n",
    format(round(state_mat[n_t, "Dead"]), big.mark = ","),
    100 * state_mat[n_t, "Dead"] / n_c))
cat("\nTREATMENT OUTCOMES:\n")
cat(sprintf("  LTBI treatment completed: %s\n", format(round(state_mat[n_t, "LatentCompleted"]), big.mark = ",")))
cat(sprintf("  Active TB cured:          %s\n", format(round(state_mat[n_t, "ActiveCompleted"]), big.mark = ",")))
cat(sprintf("  Lost to follow-up:        %s\n",
    format(round(state_mat[n_t, "LatentLtfu"] + state_mat[n_t, "ActiveLtfu"]), big.mark = ",")))
cat("\nECONOMIC ANALYSIS (with half-cycle correction):\n")
total_costs <- results$total_cost
total_qalys <- results$total_qaly
cat(sprintf("  Total discounted costs:   £%s\n", format(round(total_costs), big.mark = ",")))
cat(sprintf("  Total discounted QALYs:   %s\n", format(round(total_qalys, 1), big.mark = ",")))
cat(sprintf("  Cost per person:          £%s\n", format(round(total_costs / n_c), big.mark = ",")))
cat(sprintf("  QALYs per person:         %.2f\n", total_qalys / n_c))
cat("\n================================================================================\n\n")

# -------------------- Save Epidemiology Plots ----------------------------------

cat("Epidemiology plots saved to output/ folder\n")

# =============================================================================
# DIAGNOSTIC STRATEGY COMPARISON — BASE CASE ICER TABLE
# Runs all 43 strategies using central parameter values.
# Results are compared against no screening and ranked by ICER.
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("                    DIAGNOSTIC STRATEGY COMPARISON                              \n")
cat("================================================================================\n\n")

strategies <- create_diagnostic_strategies()

# Run all 43 strategies with base case (mean) parameter values
strategy_results <- mclapply(strategies, run_strategy,
                             mc.cores = max(1L, detectCores() - 1L))

# Calculate ICER table
icer_table <- calculate_icer(strategy_results)

cat("ICER TABLE (sequential incremental analysis with extended dominance):\n")
cat(sprintf("%-22s %11s %10s %11s %11s %14s %15s  %s\n",
            "Strategy", "Cost/person", "QALYs/pp", "Inc.Cost", "Inc.QALY",
            "ICER(vs ref)", "Sequential ICER", "Dominance"))
cat(strrep("-", 110), "\n")
for (i in 1:nrow(icer_table)) {
  icer_str     <- if (is.na(icer_table$icer[i])) "Ref"
                  else paste0("£", format(round(icer_table$icer[i]), big.mark = ","))
  seq_icer_str <- if (is.na(icer_table$sequential_icer[i])) "-"
                  else paste0("£", format(round(icer_table$sequential_icer[i]), big.mark = ","))
  cat(sprintf("%-22s %10s %10.4f %10s %11.4f %14s %15s  %s\n",
              icer_table$strategy[i],
              paste0("£", format(round(icer_table$cost_per_person[i]), big.mark = ",")),
              icer_table$qaly_per_person[i],
              paste0("£", format(round(icer_table$inc_cost[i] / n_c), big.mark = ",")),
              icer_table$inc_qaly[i] / n_c,
              icer_str,
              seq_icer_str,
              icer_table$dominance[i]))
}
cat("\n================================================================================\n")

# -------------------- Efficient frontier summary --------------------------------
frontier_rows <- icer_table[icer_table$dominance %in% c("ref", "non-dominated"), ]
frontier_rows <- frontier_rows[order(frontier_rows$cost), ]

dominated_strategies <- icer_table$strategy[icer_table$dominance == "simply dominated"]
ext_dom_strategies   <- icer_table$strategy[icer_table$dominance == "extendedly dominated"]

cat("\n")
cat("================================================================================\n")
cat("                         EFFICIENT FRONTIER SUMMARY                             \n")
cat("================================================================================\n\n")
cat("COST-EFFECTIVE FRONTIER (strategies on the non-dominated frontier):\n")
cat(sprintf("  %-24s  Cost/pp  Seq. ICER vs previous\n", "Strategy"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(frontier_rows))) {
  seq_str <- if (is.na(frontier_rows$sequential_icer[i])) "Reference"
             else paste0("£", format(round(frontier_rows$sequential_icer[i]), big.mark = ","), "/QALY")
  cat(sprintf("  %-24s  £%-7s  %s\n",
              frontier_rows$strategy[i],
              format(round(frontier_rows$cost_per_person[i]), big.mark = ","),
              seq_str))
}
cat("\n")

if (length(dominated_strategies) > 0) {
  cat(sprintf("SIMPLY DOMINATED (%d strategies — higher cost AND lower QALYs than a cheaper alternative):\n",
              length(dominated_strategies)))
  for (s in dominated_strategies) cat(sprintf("  - %s\n", s))
  cat("\n")
}
if (length(ext_dom_strategies) > 0) {
  cat(sprintf("EXTENDEDLY DOMINATED (%d strategies — worse value than linear interpolation on frontier):\n",
              length(ext_dom_strategies)))
  for (s in ext_dom_strategies) cat(sprintf("  - %s\n", s))
  cat("\n")
}

cat("INTERPRETATION (NICE WTP threshold: £25,000/QALY; updated Dec 2025, effective Apr 2026):\n")
for (i in seq_len(nrow(frontier_rows))) {
  if (is.na(frontier_rows$sequential_icer[i])) next  # skip reference
  wtp <- 25000
  ce_status <- if (frontier_rows$sequential_icer[i] <= wtp) "COST-EFFECTIVE" else "NOT cost-effective at £25k"
  cat(sprintf("  %-24s  seq. ICER £%s → %s\n",
              frontier_rows$strategy[i],
              format(round(frontier_rows$sequential_icer[i]), big.mark = ","),
              ce_status))
}
cat("\n================================================================================\n\n")

# =============================================================================
# EVENT COUNTS — FRONTIER STRATEGIES VS NO SCREENING
#
# Reports the key mechanism statistics for each frontier strategy:
#   - LTBI detected at entry (from init vector — the structural gap between
#     sequential strategies (~100/100k) and parallel IGRA (~11,526 (QFT) / 9,980 (T-SPOT) per 100k))
#   - Active TB detected at entry (TP from Zenner 2025 decision tree)
#   - Deaths over 55-year lifetime horizon (undiscounted; from Markov Dead state)
#   - Deaths averted vs No Screening
#
# =============================================================================
no_screen_idx    <- which(sapply(strategy_results, function(x) x$strategy_name == "Passive case finding"))
no_screen_res    <- strategy_results[[no_screen_idx]]
no_screen_deaths <- no_screen_res$state_membership[nrow(no_screen_res$state_membership), "Dead"]
no_screen_ltbi   <- unname(strategies[["Passive_case_finding"]]$init["LatentDiagnosed"])     # 0
no_screen_atb    <- unname(strategies[["Passive_case_finding"]]$init["ActiveDiagnosed"])     # 0

cat("================================================================================\n")
cat("               EVENT COUNTS — FRONTIER STRATEGIES VS NO SCREENING               \n")
cat("================================================================================\n\n")
cat(sprintf("  %-34s  %10s  %10s  %10s  %12s\n",
            "Strategy", "LTBI det.", "ATB det.", "Deaths", "Deaths avrt"))
cat(strrep("-", 82), "\n")

for (i in seq_len(nrow(frontier_rows))) {
  strat_nm  <- frontier_rows$strategy[i]
  strat_key <- names(strategies)[sapply(names(strategies), function(k) strategies[[k]]$name == strat_nm)]
  if (length(strat_key) == 0) next
  sr_idx    <- which(sapply(strategy_results, function(x) x$strategy_name == strat_nm))
  if (length(sr_idx) == 0) next
  sr        <- strategy_results[[sr_idx]]
  deaths_i  <- sr$state_membership[nrow(sr$state_membership), "Dead"]
  ltbi_det  <- unname(strategies[[strat_key[1]]]$init["LatentDiagnosed"])
  atb_det   <- unname(strategies[[strat_key[1]]]$init["ActiveDiagnosed"])
  cat(sprintf("  %-34s  %10.0f  %10.0f  %10.0f  %12.0f\n",
              strat_nm,
              ltbi_det, atb_det, deaths_i,
              no_screen_deaths - deaths_i))
}
cat("\n")
cat("  Notes: LTBI detected = persons in LatentDiagnosed state at t=0 (entry screening)\n")
cat("         ATB detected  = true positives (TP) from Zenner 2025 decision tree\n")
cat("         Deaths        = undiscounted cumulative deaths over 55-year horizon\n\n")

# Save event counts for all strategies
event_counts_df <- bind_rows(lapply(names(strategies), function(k) {
  sr_idx <- which(sapply(strategy_results, function(x) x$strategy_name == strategies[[k]]$name))
  if (length(sr_idx) == 0) return(NULL)
  sr <- strategy_results[[sr_idx]]
  tibble(
    strategy        = strategies[[k]]$name,
    ltbi_detected   = unname(strategies[[k]]$init["LatentDiagnosed"]),
    active_detected = unname(strategies[[k]]$init["ActiveDiagnosed"]),
    deaths_55yr     = sr$state_membership[nrow(sr$state_membership), "Dead"],
    deaths_averted  = no_screen_deaths - sr$state_membership[nrow(sr$state_membership), "Dead"]
  )
}))
write.csv(event_counts_df, "output/csv/event_counts_basecase.csv", row.names = FALSE)
cat("Saved: output/csv/event_counts_basecase.csv\n\n")

# Save base case ICER table as CSV (all 43 strategies), with detection columns added
write.csv(icer_table %>%
  mutate(cost_per_person = cost / n_c, qaly_per_person = qaly / n_c) %>%
  left_join(event_counts_df %>% select(strategy, ltbi_detected, active_detected),
            by = "strategy"),
  "output/csv/icer_nhs_perspective.csv", row.names = FALSE)
cat("Saved: output/csv/icer_nhs_perspective.csv\n\n")

# Non-dominated strategy names (ref + efficient frontier) — used to focus labels
# and filter PSA plots. 43 strategies make fully-labelled/coloured plots unreadable.
non_dominated_names <- icer_table$strategy[
  icer_table$dominance %in% c("ref", "non-dominated")
]

# -------------------- Base case cost-effectiveness scatter plot ----------------
# Each point is coloured by dominance status:
#   dark rose = reference (Passive case finding) or efficient frontier
#   light pink = dominated (simply or extendedly)
# No WTP line — this is total cost vs total QALYs (not incremental); WTP line
# appears in the incremental CE plane (plot_ce_plane) and ICER forest plot.
icer_table$ce_status <- dplyr::case_when(
  icer_table$dominance == "ref"                  ~ "Reference",
  icer_table$dominance == "non-dominated"        ~ "Efficient frontier",
  icer_table$dominance == "extendedly dominated" ~ "Extendedly dominated",
  TRUE                                           ~ "Dominated"
)
# Shape aesthetic: triangle for parallel IGRA strategies, circle for all others
icer_table$strategy_type <- ifelse(
  grepl("^Parallel", icer_table$strategy, ignore.case = FALSE),
  "Parallel IGRA", "Sequential"
)

p_icer <- ggplot(icer_table, aes(x = qaly / n_c, y = cost / n_c,
                                   colour = ce_status, shape = strategy_type)) +
  geom_point(aes(size = strategy_type)) +
  ggrepel::geom_text_repel(
    data   = icer_table %>% filter(dominance %in% c("ref", "non-dominated")),
    aes(label = strategy,
        fontface = ifelse(strategy_type == "Parallel IGRA", "bold", "plain")),
    size = 3, max.overlaps = Inf, colour = "grey20",
    box.padding = 0.8, point.padding = 0.5, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Reference" = "#ad002a", "Efficient frontier" = "#ed0000",
               "Extendedly dominated" = "#6b7b8d", "Dominated" = "#c8cfd6"),
    breaks = c("Reference", "Efficient frontier", "Extendedly dominated", "Dominated"),
    name   = NULL
  ) +
  scale_shape_manual(
    values = c("Parallel IGRA" = 17L, "Sequential" = 16L),
    name   = NULL
  ) +
  scale_size_manual(
    values = c("Parallel IGRA" = 4.5, "Sequential" = 3.5),
    guide  = "none"
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "QALYs per person", y = "Cost per person (\u00a3)",
       title = NULL,
       subtitle = NULL) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")


# =============================================================================
# STRATEGY SPACE VISUALISATION
#
# Plots all strategies from Zenner et al. 2025 Eur Respir J decision tree to contextualise the
# strategies selected for the Markov analysis.  Uses excel_ic already loaded
# at the top of the script.
#
# Plot 1 — scatter: TP detected (x) vs diagnostic cost per person (y),
#           coloured by test family; selected strategies labelled.
# Plot 2 — horizontal bar: TP detected, ranked descending.
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("         STRATEGY SPACE VISUALISATION (all Excel strategies)                   \n")
cat("================================================================================\n\n")

# Build data frame from all rows in excel_ic ----------------------------------------
strategy_space <- data.frame(
  name      = as.character(excel_ic[[1]]),
  TP        = suppressWarnings(as.numeric(excel_ic[[4]])),
  Tot_costs = suppressWarnings(as.numeric(excel_ic[[6]])),
  stringsAsFactors = FALSE
)
strategy_space <- strategy_space[
  !is.na(strategy_space$TP) & !is.na(strategy_space$Tot_costs) &
  strategy_space$name != "" & !is.na(strategy_space$name), ]
strategy_space$cost_per_person <- strategy_space$Tot_costs / 100000

strategy_space$family <- sapply(strategy_space$name, assign_family)

# All feasible strategies (non-NA costs) are included in the Markov analysis
strategy_space$selected <- !is.na(strategy_space$Tot_costs)

strategy_space$display_name <- sapply(strategy_space$name, clean_excel_name)

# --- Plot 1: Scatter — TP detected vs diagnostic cost per person -----------------
# Points coloured by primary screening test family.
family_pal  <- c(
  "CXR / Symptom screen"    = "#fdaf91",
  "QFT-GIT"                 = "#0099b4",
  "T-SPOT.TB"               = "#00468b",
  "TST Mantoux"             = "#ad002a",
  "Molecular (Xpert/Ultra)" = "#ad002a"
)

p_space_scatter <- ggplot(
  strategy_space,
  aes(x = TP, y = cost_per_person, colour = family)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = display_name), size = 2.5,
                  max.overlaps = 20, seed = 42, colour = "#333333",
                  show.legend = FALSE) +
  scale_colour_manual(values = family_pal, name = "Primary screening test") +
  theme_minimal(base_size = 13) +
  labs(
    x     = "True positives detected (active TB, per 100,000 screened)",
    y     = "Diagnostic cost per person (\u00a3)",
    title = NULL
  ) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )


# --- Plot 2: Two-panel bar — Active TB detection (left) + LTBI detection (right)
# Uses event_counts_df (all 43 strategies); colour-coded by pathway type.
# Parallel IGRA strategies: medium purple (QFT) / dark purple (T-SPOT) — consistent with CE scatter.
# TST is NOT an IGRA — only QFT-GIT and T-SPOT.TB are IGRAs.
# TST strategies grouped with non-IGRA (CXR/symptom screen).
classify_igra <- function(nm) {
  if (grepl("^Parallel", nm, ignore.case = TRUE)) return("Parallel IGRA (QFT-GIT / T-SPOT.TB)")
  if (grepl("QFT-GIT|T-SPOT", nm, ignore.case = TRUE)) return("IGRA-inclusive (sequential pathway)")
  return("Non-IGRA (CXR / symptom screen / TST)")
}

igra_pal_bar <- c(
  "Non-IGRA (CXR / symptom screen / TST)"  = "#adb6b6",
  "IGRA-inclusive (sequential pathway)"     = "#0099b4",
  "Parallel IGRA (QFT-GIT / T-SPOT.TB)"    = "#00468b"
)

combined_bar_df <- event_counts_df %>%
  filter(strategy != "Passive case finding") %>%
  mutate(
    igra_type = factor(sapply(strategy, classify_igra),
      levels = c("Non-IGRA (CXR / symptom screen / TST)",
                 "IGRA-inclusive (sequential pathway)",
                 "Parallel IGRA (QFT-GIT / T-SPOT.TB)"))
  ) %>%
  arrange(igra_type, active_detected) %>%
  mutate(strategy_fct = factor(strategy, levels = strategy))

shared_bar_fill <- scale_fill_manual(values = igra_pal_bar, name = NULL)

p_active_bar <- ggplot(combined_bar_df,
    aes(x = active_detected / n_c * 100, y = strategy_fct, fill = igra_type)) +
  geom_col(width = 0.7, colour = NA) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08)),
    labels = function(x) paste0(x, "%")
  ) +
  shared_bar_fill +
  theme_minimal(base_size = 9.5) +
  labs(
    x     = "Active TB cases detected (% of cohort)",
    y     = NULL,
    title = NULL
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 7)
  )

p_ltbi_bar <- ggplot(combined_bar_df,
    aes(x = ltbi_detected / n_c * 100, y = strategy_fct, fill = igra_type)) +
  geom_col(width = 0.7, colour = NA) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = function(x) paste0(x, "%")
  ) +
  shared_bar_fill +
  theme_minimal(base_size = 9.5) +
  labs(
    x     = "LTBI cases detected (% of cohort)",
    y     = NULL,
    title = NULL
  ) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank()
  )

p_space_bar <- (p_active_bar + p_ltbi_bar) +
  plot_layout(guides = "collect", widths = c(1.15, 1)) &
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.8, "lines"))

p_space_bar <- p_space_bar +
  plot_annotation(
    subtitle = "42 active screening strategies shown (Passive case finding reference excluded) \u00b7 grouped by test type \u00b7 sorted by active TB detection within each group \u00b7 percentage of cohort (n = 100,000)",
    theme = theme(
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )


# =============================================================================
# LTBI DETECTION GAP FIGURE
#
# Horizontal bar chart showing LTBI detected at entry per 100,000 across all
# 43 strategies, coloured by pathway type. This is the central structural
# finding: sequential strategies detect ~72–116 LTBI/100k (driven by the
# decision tree's sequential pathway design — IGRA applied only to
# symptom-screen positives); parallel IGRA strategies detect 14,774–15,664/100k
# (IGRA applied to all migrants). The ~100× gap is the mechanism behind the
# frontier composition.
# =============================================================================

# Classify each strategy by pathway type for colouring
classify_pathway <- function(nm) {
  if (nm == "Passive case finding")                          return("Passive case finding")
  if (grepl("^Parallel.*QFT",  nm, ignore.case = TRUE)) return("Parallel IGRA (QFT-Plus)")
  if (grepl("^Parallel.*T-SPOT|^Parallel.*TSPT|^Parallel.*tspt", nm, ignore.case = TRUE)) return("Parallel IGRA (T-SPOT.TB)")
  if (grepl("QFT-GIT|T-SPOT", nm, ignore.case = TRUE)) return("Sequential IGRA")
  return("Symptom screen / CXR")  # TST is a skin test, NOT an IGRA
}

validate_model_correctness <- function(strategies, event_counts_df, icer_table) {
  errors <- character(0)

  # --- 1. IGRA classification invariants ---
  # TST is NOT an IGRA. Only QFT-GIT and T-SPOT.TB are IGRAs.
  all_names <- event_counts_df$strategy

  tst_strategies <- all_names[grepl("TST", all_names, ignore.case = FALSE)]
  for (nm in tst_strategies) {
    grp <- classify_igra(nm)
    if (grp != "Non-IGRA (CXR / symptom screen / TST)")
      errors <- c(errors, sprintf("classify_igra: TST strategy '%s' -> '%s' (should be Non-IGRA)", nm, grp))
    grp2 <- classify_pathway(nm)
    if (grp2 %in% c("Sequential IGRA", "Parallel IGRA (QFT-Plus)", "Parallel IGRA (T-SPOT.TB)"))
      errors <- c(errors, sprintf("classify_pathway: TST strategy '%s' -> '%s' (TST is not an IGRA)", nm, grp2))
  }

  # All Parallel strategies must be in the Parallel IGRA group
  parallel_names <- all_names[grepl("^Parallel", all_names)]
  for (nm in parallel_names) {
    grp <- classify_igra(nm)
    if (!grepl("Parallel IGRA", grp))
      errors <- c(errors, sprintf("classify_igra: Parallel strategy '%s' -> '%s' (should be Parallel IGRA)", nm, grp))
    grp2 <- classify_pathway(nm)
    if (!grepl("Parallel IGRA", grp2))
      errors <- c(errors, sprintf("classify_pathway: Parallel strategy '%s' -> '%s' (should be Parallel IGRA)", nm, grp2))
  }

  # Xpert and Ultra (without QFT/T-SPOT/Parallel) must be Non-IGRA
  mol_only <- all_names[grepl("(Xpert|Ultra)$", all_names) &
                         !grepl("QFT|T-SPOT|^Parallel", all_names)]
  for (nm in mol_only) {
    grp <- classify_igra(nm)
    if (grp != "Non-IGRA (CXR / symptom screen / TST)")
      errors <- c(errors, sprintf("classify_igra: Molecular-only strategy '%s' -> '%s' (should be Non-IGRA)", nm, grp))
  }

  # --- 2. Colour palette coverage ---
  # Every named strategy in a palette must exist in the data; no silent grey points
  palettes_to_check <- list(
    ceac_strat_cols    = c("Passive case finding", "Cough+CXR (TB sx)",
                           "Symptom screen+CXR", "Parallel Sx+QFT (Ultra)"),
    ce_ellipse_cols    = c("Cough+CXR (TB sx)", "Symptom screen+CXR",
                           "Parallel Sx+QFT (Ultra)"),
    strat_cols_ov      = c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                           "Parallel Cough+T-SPOT (Ultra)", "Parallel Sx+T-SPOT (Ultra)"),
    strat_cols_uptake  = c("Parallel Cough+QFT (Ultra)", "Parallel Sx+QFT (Ultra)",
                           "Parallel Cough+T-SPOT (Ultra)", "Parallel Sx+T-SPOT (Ultra)")
  )
  for (pal_name in names(palettes_to_check)) {
    for (strat in palettes_to_check[[pal_name]]) {
      if (!strat %in% all_names)
        errors <- c(errors, sprintf("Palette '%s' key '%s' not found in strategy names", pal_name, strat))
    }
  }

  # --- 3. Frontier ICERs are cost-effective at NICE threshold ---
  frontier <- icer_table[icer_table$dominance %in% c("ref", "non-dominated"), ]
  for (i in seq_len(nrow(frontier))) {
    icer <- frontier$icer[i]
    if (!is.na(icer) && icer > 25000)
      errors <- c(errors, sprintf("Frontier strategy '%s' has ICER L%s > L25,000/QALY",
                                  frontier$strategy[i], format(round(icer), big.mark=",")))
  }

  # --- 4. Passive case finding must be the lowest-cost strategy ---
  ns_cost <- icer_table$cost[icer_table$strategy == "Passive case finding"]
  if (length(ns_cost) == 0)
    errors <- c(errors, "'Passive case finding' strategy not found in icer_table")
  else if (any(icer_table$cost < ns_cost - 1))
    errors <- c(errors, sprintf("Strategy cheaper than Passive case finding found: %s (L%.0f vs L%.0f)",
                                icer_table$strategy[which.min(icer_table$cost)],
                                min(icer_table$cost), ns_cost))

  # --- 5. Dominance values are valid ---
  valid_dom <- c("ref", "non-dominated", "simply dominated", "extendedly dominated")
  bad_dom <- setdiff(unique(icer_table$dominance), valid_dom)
  if (length(bad_dom) > 0)
    errors <- c(errors, paste("Invalid dominance values:", paste(bad_dom, collapse=", ")))

  # --- 6. Strategy count ---
  n_strat <- nrow(icer_table)
  if (n_strat != 43)
    errors <- c(errors, sprintf("Expected 43 strategies, found %d", n_strat))

  # --- Report ---
  if (length(errors) > 0) {
    cat("\n")
    cat("================================================================================\n")
    cat("  MODEL CORRECTNESS VALIDATION FAILED\n")
    cat("================================================================================\n")
    for (e in errors) cat(sprintf("  x %s\n", e))
    cat("\n")
    stop("Validation failed -- fix errors above before proceeding.", call. = FALSE)
  }
  cat("- Model correctness validation passed (", length(all_names), " strategies checked)\n", sep="")
}

validate_model_correctness(strategies, event_counts_df, icer_table)

ltbi_gap_df <- event_counts_df %>%
  mutate(
    pathway_type = sapply(strategy, classify_pathway),
    pathway_type = factor(pathway_type,
      levels = c("Passive case finding", "Symptom screen / CXR",
                 "Sequential IGRA", "Parallel IGRA (QFT-Plus)", "Parallel IGRA (T-SPOT.TB)"))
  )

pathway_pal <- c(
  "Passive case finding"       = "#adb6b6",  # gray
  "Symptom screen / CXR"      = "#c0cad2",  # light steel gray
  "Sequential IGRA"            = "#0099b4",  # teal
  "Parallel IGRA (QFT-Plus)"   = "#00468b",  # navy
  "Parallel IGRA (T-SPOT.TB)"  = "#3a5a7a"   # steel navy
)

# Summary per pathway type: median bar + individual strategy dots
ltbi_gap_summary <- ltbi_gap_df %>%
  group_by(pathway_type) %>%
  summarise(median_ltbi = median(ltbi_detected), .groups = "drop")

p_ltbi_gap <- ggplot() +
  geom_col(data = ltbi_gap_summary,
           aes(x = median_ltbi, y = pathway_type, fill = pathway_type),
           width = 0.55, colour = NA) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::comma,
    breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000)
  ) +
  scale_fill_manual(values = pathway_pal, guide = "none") +
  theme_minimal(base_size = 12) +
  labs(
    x = "LTBI cases detected at entry (per 100,000 migrants screened)",
    y = NULL
  ) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 12)
  )

ggsave("output/supplementary/ltbi_detection_gap_supplementary.tiff", p_ltbi_gap,
       width = 10, height = 5, dpi = 300)
cat("Saved: output/supplementary/ltbi_detection_gap_supplementary.tiff\n")

# =============================================================================

source("scripts/run_sa.R")
# =============================================================================
# PSA SUMMARY — EFFICIENT FRONTIER STRATEGIES
#
# Extracts PSA credible intervals for the five efficient-frontier strategies
# (No Screening reference + four non-dominated strategies). Placed after the
# full 43-strategy table to provide a reader-friendly focal summary.
#
# Key contrast highlighted here:
#   Cough+CXR (TB sx):      ~100% P(CE at £25k) — ICER driven by active TB detection
#   Symptom screen+CXR:     ~100% P(CE at £25k) — similar driver
#   Parallel IGRA strategies: lower P(CE) — driven by LTBI treatment cascade
#     (initiation, completion, prevention of reactivation — all sampled in PSA)
# =============================================================================
cat("\n================================================================================\n")
cat("          PSA SUMMARY — EFFICIENT FRONTIER STRATEGIES ONLY                      \n")
cat("================================================================================\n\n")
cat("  (Full 43-strategy table above; this block highlights the narrative contrast.)\n\n")

frontier_psa <- icer_ci %>%
  filter(strategy %in% non_dominated_names) %>%
  arrange(mean_inc_cost)

cat(sprintf("  %-34s  %12s  %26s  %10s\n",
            "Strategy", "Median ICER", "95% CrI", "P(CE £25k)"))
cat(strrep("-", 88), "\n")
for (i in seq_len(nrow(frontier_psa))) {
  med_str <- if (is.na(frontier_psa$median_icer[i])) "Reference"
             else paste0("\u00a3", format(round(frontier_psa$median_icer[i]), big.mark = ","))
  lo_str  <- format(round(frontier_psa$icer_lo[i]), big.mark = ",")
  hi_str  <- format(round(frontier_psa$icer_hi[i]), big.mark = ",")
  cri_str <- paste0("(\u00a3", lo_str, " to \u00a3", hi_str, ")")
  pce_str <- if (is.na(frontier_psa$prob_ce_25k[i])) "—"
             else sprintf("%.1f%%", 100 * frontier_psa$prob_ce_25k[i])
  cat(sprintf("  %-34s  %12s  %26s  %10s\n",
              frontier_psa$strategy[i], med_str, cri_str, pce_str))
}

# Derive summary stats dynamically from data for interpretation block
.igra_row   <- frontier_psa[frontier_psa$strategy == "Parallel Sx+QFT (Ultra)", ]
.igra_pce   <- if (nrow(.igra_row) > 0) sprintf("%.1f%%", 100 * .igra_row$prob_ce_25k) else "—"
.igra_pdom  <- if (nrow(.igra_row) > 0) sprintf("%.1f%%", 100 * .igra_row$prob_dominated) else "—"

cat("\n")
cat("  Interpretation:\n")
cat("  - Non-IGRA frontier strategies (Cough+CXR, Symptom screen+CXR) have near-certain\n")
cat("    CE — ICER driven by active TB detection, which is stable across PSA draws.\n")
cat(sprintf("  - Parallel Sx+QFT (Ultra) P(CE £25k) = %s; P(dominated vs No Screening) = %s.\n",
            .igra_pce, .igra_pdom))
cat("    Wider uncertainty because QALY advantage depends on LTBI treatment cascade\n")
cat("    (initiation, completion, prevention of reactivation), each sampled independently.\n")
cat("\n================================================================================\n")

# Save ICER credible interval table for results reporting
icer_ci_out <- icer_ci %>%
  mutate(
    mean_inc_cost_pp = mean_inc_cost / n_c,
    mean_inc_qaly_pp = mean_inc_qaly / n_c
  )
write.csv(icer_ci_out, "output/csv/icer_confidence_intervals.csv", row.names = FALSE)
cat("ICER CI table saved to output/csv/icer_confidence_intervals.csv\n")

# -------------------- ICER forest plot with 95% credible intervals ------------
# ICERs are capped at ±£200,000 for axis readability; strategies with very
# high uncertainty (e.g. dominated in most simulations) can produce extreme
# ICER values that would otherwise compress the scale. Capping is standard
# practice in health economic visualisation and does not affect interpretation.
icer_ci_plot <- icer_ci %>%
  mutate(
    median_icer_cap = pmax(pmin(median_icer, 200000), -200000),
    icer_lo_cap     = pmax(pmin(icer_lo, 200000), -200000),
    icer_hi_cap     = pmax(pmin(icer_hi, 200000), -200000),
    strategy        = reorder(strategy, median_icer)
  )

p_icer_forest <- ggplot(icer_ci_plot, aes(x = median_icer_cap, y = strategy)) +
  geom_vline(xintercept = 25000, linetype = "dashed", color = "#ad002a", linewidth = 0.8) +
  geom_vline(xintercept = 35000, linetype = "dashed", color = "#ed0000", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = icer_lo_cap, xmax = icer_hi_cap),
                width = 0.3, linewidth = 0.8, color = "grey40", orientation = "y") +
  geom_point(size = 4, color = "#ed0000") +
  annotate("label", x = 25000, y = Inf, vjust = 1.3, label = "NICE £25k\nWTP threshold",
           color = "#ad002a", fill = "white", alpha = 0.85,
           size = 3, hjust = 1.05, fontface = "bold", linewidth = 0.3) +
  annotate("label", x = 35000, y = Inf, vjust = 1.3, label = "NICE £35k\nWTP threshold",
           color = "#ed0000", fill = "white", alpha = 0.85,
           size = 3, hjust = -0.05, fontface = "bold", linewidth = 0.3) +
  scale_x_continuous(labels = function(x) paste0("£", trimws(format(x, big.mark = ",", scientific = FALSE)))) +
  theme_minimal(base_size = 14) +
  labs(x = "ICER (£/QALY)", y = "",
       title = NULL,
       subtitle = "Median ICER from 1,000 PSA simulations vs No Screening") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9))


# -------------------- Probability of cost-effectiveness stacked bar chart -----
# Each bar shows the proportion of PSA simulations in which the strategy was
# cost-effective (ICER < £25,000/QALY), not cost-effective, or dominated.
pce_data <- icer_ci %>%
  select(strategy, prob_ce_25k, prob_dominated) %>%
  mutate(
    prob_not_ce = 1 - prob_ce_25k - prob_dominated,
    strategy = reorder(strategy, prob_ce_25k)
  ) %>%
  pivot_longer(cols = c(prob_ce_25k, prob_dominated, prob_not_ce),
               names_to = "outcome", values_to = "prob") %>%
  mutate(outcome = factor(outcome,
    levels = c("prob_ce_25k", "prob_not_ce", "prob_dominated"),
    labels = c("Cost-effective (ICER < £25k)", "Not cost-effective", "Dominated")))

p_pce_bar <- ggplot(pce_data, aes(x = strategy, y = prob, fill = outcome)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Cost-effective (ICER < £25k)" = "#ed0000",
                                "Not cost-effective" = "#fdaf91",
                                "Dominated" = "#adb6b6")) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Proportion of PSA simulations",
       title = NULL,
       subtitle = "At £25,000/QALY WTP threshold | 1,000 PSA simulations",
       fill = "") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size = 9))


# -------------------- Cost-effectiveness plane with 95% confidence ellipses ---
# Each cloud of points represents one strategy across 1,000 PSA simulations.
# The 95% ellipse captures joint cost-QALY uncertainty. Diamond markers show
# the mean position. The dotted line is the £25,000/QALY WTP frontier.
icer_psa_pp <- icer_psa %>%
  mutate(inc_cost_pp = inc_cost / n_c, inc_qaly_pp = inc_qaly / n_c)

icer_means <- icer_psa_pp %>%
  group_by(strategy) %>%
  summarise(inc_cost_pp = mean(inc_cost_pp), inc_qaly_pp = mean(inc_qaly_pp), .groups = "drop")

all_strat_names  <- unique(psa_df$strategy)
ce_plane_colours <- build_strategy_colours(all_strat_names)

# For ellipse / CE plane / CEAC: filter to non-dominated strategies only.
# 35 strategies produce unreadable overlapping clouds and legends.
# Standard HTA practice: show efficient frontier + reference on these plots;
# all 43 are shown in the strategy space scatter and forest/bar plots.
nd_strats_psa    <- non_dominated_names[non_dominated_names != "Passive case finding"]
icer_psa_nd      <- icer_psa_pp %>% filter(strategy %in% nd_strats_psa)
icer_means_nd    <- icer_means   %>% filter(strategy %in% nd_strats_psa)
ce_ellipse_cols  <- build_strategy_colours(nd_strats_psa)

# Distinct colours for CE plane — match CEAC palette for consistency
ce_ellipse_named_cols <- c(
  "Cough+CXR (TB sx)"          = "#ed0000",
  "Symptom screen+CXR"         = "#0099b4",
  "Parallel Sx+QFT (Ultra)"    = "#00468b"
)
extra_ellipse <- setdiff(nd_strats_psa, names(ce_ellipse_named_cols))
if (length(extra_ellipse) > 0) {
  extra_ce_cols <- setNames(project_pal(length(extra_ellipse)), extra_ellipse)
  ce_ellipse_named_cols <- c(ce_ellipse_named_cols, extra_ce_cols)
}

p_ce_ellipse <- ggplot(icer_psa_nd, aes(x = inc_qaly_pp, y = inc_cost_pp, color = strategy)) +
  geom_point(alpha = 0.15, size = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(data = icer_means_nd, size = 4, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(intercept = 0, slope = 25000, linetype = "dotted", color = "black", linewidth = 1) +
  # Annotate the Cough+CXR cluster which sits at the origin and is too small
  # to see at this scale (inc. QALYs ~0.006 pp; inc. cost ~£7 pp).
  annotate("segment", x = 0.005, xend = 0.001, y = 30, yend = 5,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           colour = "#ed0000", linewidth = 0.8) +
  annotate("text", x = 0.006, y = 38,
           label = "Cough+CXR (TB sx)\n(inc. cost ~\u00a37, inc. QALY ~0.006;\nsee zoomed plot)",
           colour = "#ed0000", size = 3, hjust = 0) +
  scale_colour_manual(values = ce_ellipse_named_cols) +
  scale_x_continuous(labels = function(x) format(round(x, 4), nsmall = 4, scientific = FALSE)) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), "\u00a3", scales::comma(abs(x)))) +
  theme_minimal(base_size = 14) +
  labs(x = "Incremental QALYs per person",
       y = "Incremental cost per person (\u00a3)",
       title = NULL,
       subtitle = "Efficient frontier strategies only | vs No Screening | Dotted = \u00a325,000/QALY | 1,000 PSA simulations",
       color = "Strategy") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "right")


# Zoomed CE plane — shows the Cough+CXR (TB sx) cluster near the origin
# which is invisible in the full-scale ellipses plot above.
p_ce_ellipse_zoom <- ggplot(icer_psa_nd, aes(x = inc_qaly_pp, y = inc_cost_pp, color = strategy)) +
  geom_point(alpha = 0.25, size = 1) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  geom_point(data = icer_means_nd, size = 4, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(intercept = 0, slope = 25000, linetype = "dotted", color = "black", linewidth = 1) +
  coord_cartesian(xlim = c(-0.005, 0.020), ylim = c(-20, 60)) +
  scale_colour_manual(values = ce_ellipse_named_cols) +
  scale_x_continuous(labels = function(x) format(round(x, 4), nsmall = 4, scientific = FALSE)) +
  scale_y_continuous(labels = function(x) paste0(ifelse(x < 0, "-", ""), "\u00a3", scales::comma(abs(x)))) +
  theme_minimal(base_size = 14) +
  labs(x = "Incremental QALYs per person",
       y = "Incremental cost per person (\u00a3)",
       title = NULL,
       subtitle = "Zoomed to origin region | Parallel IGRA ellipses extend far right (see full-scale plot)",
       color = "Strategy") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "grey40"),
        panel.grid.minor = element_blank(),
        legend.position = "right")

cat("ICER visualisations saved to output/ folder\n")

# CE plane and CEAC: filter PSA data to non-dominated strategies for readability.
# "Passive case finding" (reference) must be included for incremental calculations.
psa_df_nd <- psa_df %>% filter(strategy %in% non_dominated_names)

# Cost-effectiveness plane: scatter of incremental cost vs incremental QALY
# across all 1,000 simulations for efficient frontier strategies vs no screening.
p_ce_plane <- plot_ce_plane(psa_df_nd, reference = "Passive case finding", wtp = 25000)

# NHS CEAC omitted: body Figure 2 (ceac_panel.tiff) shows the NHS perspective
# as a dashed overlay on the societal CEAC.

cat("\nPSA plots saved to output/ folder\n")

# =============================================================================
# SUPPLEMENTARY EXCEL FILES — S4 and S5
# =============================================================================
# S4: All 43 strategies × all 3 conditions (base case + 2 scenarios)
# S5: All sensitivity analyses (PSA, one-way DSA) × all scenarios
#
# Both files are written to output/supplementary/ folder.
# All data come from CSVs already written to output/csv/ in this run.
# =============================================================================

library(openxlsx)

.today_str <- format(Sys.Date(), "%Y-%m-%d")

# Helper: read a CSV and write it as a formatted table into a workbook sheet.
# Returns the data frame (invisibly) so callers can inspect it if needed.
.write_supp_sheet <- function(wb, sheet_name, csv_path, description) {
  df <- tryCatch(
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    cat(sprintf("  [SKIP] %s not found\n", csv_path))
    return(invisible(NULL))
  }
  addWorksheet(wb, sheetName = sheet_name)

  # Row 1: description note
  writeData(wb, sheet_name, description,
            startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name,
           createStyle(fontSize = 9, fontColour = "#555555", wrapText = TRUE),
           rows = 1, cols = 1, gridExpand = FALSE)
  setRowHeights(wb, sheet_name, rows = 1, heights = 28)

  # Row 2+: data table with auto-filter
  writeDataTable(wb, sheet_name, df,
                 startRow = 2, startCol = 1,
                 tableStyle = "TableStyleLight9",
                 withFilter  = TRUE)
  setColWidths(wb, sheet_name,
               cols   = seq_len(ncol(df)),
               widths = "auto")
  invisible(df)
}

# ---------------------------------------------------------------------------
# S4 — All strategies, all conditions
# ---------------------------------------------------------------------------
.s4_path <- file.path("output/supplementary", "CE_Results_supplementary.xlsx")
.wb4 <- createWorkbook()
modifyBaseFont(.wb4, fontSize = 10, fontName = "Calibri")

.write_supp_sheet(
  .wb4, "Base case",
  "output/csv/icer_nhs_perspective.csv",
  paste0("S4a · Base case — all 43 strategies · 55-year lifetime horizon · ",
         "n = 100,000 cohort · 2023/24 GBP · ", .today_str,
         " · Columns: cost and QALYs per cohort; cost/QALY per person; ",
         "ICER vs No Screening; sequential ICER; dominance classification; ",
         "LTBI detected; active TB detected.")
)

.write_supp_sheet(
  .wb4, "Primary — Transmission",
  "output/csv/icer_societal_perspective.csv",
  paste0("S4b · Primary analysis (Transmission prevention, societal perspective) — all 43 strategies · ",
         "5 contacts/case/year (conservative base) × 4.1% annual transmission risk = ",
         "β = 0.205 infections/year · £6,055/secondary case prevented · ",
         "All other parameters at base-case values · 2023/24 GBP · ", .today_str)
)

saveWorkbook(.wb4, .s4_path, overwrite = TRUE)
cat(sprintf("Saved: %s\n", .s4_path))

# ---------------------------------------------------------------------------
# S5 — All sensitivity analyses, DSA and PSA
# ---------------------------------------------------------------------------
.s5_path <- file.path("output/supplementary", "Sensitivity_Analyses_supplementary.xlsx")
.wb5 <- createWorkbook()
modifyBaseFont(.wb5, fontSize = 10, fontName = "Calibri")

# --- PSA ---
.write_supp_sheet(
  .wb5, "PSA — Base case summary",
  "output/csv/icer_confidence_intervals.csv",
  paste0("S5a · PSA base case — ICER 95% credible intervals and P(CE) for all 43 strategies · ",
         "1,000 simulations · Beta distributions for proportions/utilities; ",
         "Gamma distributions for costs · WTP thresholds £25,000 and £35,000/QALY · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "PSA — Raw simulations",
  "output/csv/psa_results.csv",
  paste0("S5b · PSA raw results — cost and QALYs per cohort (100,000) for each of ",
         "1,000 simulations × 43 strategies (43,000 rows) · ",
         "sim = simulation index; cost = total cohort cost (£); qaly = total cohort QALYs · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "PSA — Transmission scenario",
  "output/csv/transmission_psa_summary.csv",
  paste0("S5c · Primary analysis PSA summary (transmission prevention, societal perspective) — frontier strategies · ",
         "Contacts per case/year drawn from Uniform(5, 15) across 1,000 simulations · ",
         "β = contacts × 4.1% annual transmission risk · £6,055/secondary case · ",
         "2023/24 GBP · ", .today_str)
)

# --- DSA ---
.write_supp_sheet(
  .wb5, "DSA — Cough+CXR (TB sx)",
  "output/csv/dsa_results_coughcxr.csv",
  paste0("S5d · One-way DSA — Cough+CXR (TB sx) vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Symptom screen+CXR",
  "output/csv/dsa_results_symscrCXR.csv",
  paste0("S5e · One-way DSA — Symptom screen+CXR vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Parallel Sx+QFT (Ultra)",
  "output/csv/dsa_results.csv",
  paste0("S5f · One-way DSA — Parallel Sx+QFT (Ultra) vs No Screening · ",
         "NMB at £25,000/QALY WTP · Each parameter varied independently between ",
         "its lower and upper bound; all others held at base-case values · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "DSA — Transmission contacts",
  "output/csv/transmission_dsa_results.csv",
  paste0("S5g · Primary analysis DSA — close contacts per active TB case · ",
         "Values: 3, 5 (base), 7, 10, 15 contacts/case/year · ",
         "β = contacts × 4.1% per year · £6,055/secondary case · ",
         "2023/24 GBP · ", .today_str)
)

# --- Structural SAs ---
.write_supp_sheet(
  .wb5, "SA — LTBI completion",
  "output/csv/ltbi_completion_sensitivity.csv",
  paste0("S5h · Structural SA — LTBI treatment completion · ",
         "Scenarios: trial 76.7% (base; UKHSA 2024 TB Prevention report) vs ",
         "real-world 55.5% (UKHSA 2024) · 2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — LTBI prevalence",
  "output/csv/ltbi_prevalence_sensitivity.csv",
  paste0("S5i · Structural SA — LTBI prevalence at entry · ",
         "Scenarios: base 17.8% (UKHSA 2021 migrant cohort) vs ",
         "15.1% (UKHSA 2025 updated estimate) · 2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — Active TB prevalence",
  "output/csv/active_tb_prev_sensitivity.csv",
  paste0("S5j · Structural SA — Active TB prevalence at entry · ",
         "Scenarios: base 1.0%, high-risk 0.44%, pooled low-burden 0.215% · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA uptake",
  "output/csv/igra_uptake_sensitivity.csv",
  paste0("S5k · Structural SA — IGRA programme uptake rate · ",
         "Varied 70% to 100% · Applies to all parallel IGRA strategies · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA specificity",
  "output/csv/igra_specificity_sensitivity.csv",
  paste0("S5l · Structural SA — IGRA test specificity · ",
         "QFT-GIT base 0.96; T-SPOT.TB base 0.93; varied 0.90 to 1.00 · ",
         "2023/24 GBP · ", .today_str)
)

.write_supp_sheet(
  .wb5, "SA — IGRA programme cost",
  "output/csv/igra_programme_cost_sensitivity.csv",
  paste0("S5m · Structural SA — IGRA programme delivery cost per person screened · ",
         "Varied across plausible range (£4 to £60) · ",
         "Cost added to initial screening cost only; does not affect Markov transitions · ",
         "2023/24 GBP · ", .today_str)
)

saveWorkbook(.wb5, .s5_path, overwrite = TRUE)
cat(sprintf("Saved: %s\n", .s5_path))

cat("\n================================================================================\n")
cat("                    GENERATING FIGURES AND TABLES                               \n")
cat("================================================================================\n")

source("scripts/run_plots_and_tables.R")

cat("\n================================================================================\n")
cat("                    ALL ANALYSES COMPLETE                                       \n")
cat("================================================================================\n")
