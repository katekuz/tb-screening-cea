# Model Assumptions — TB Screening CEA
# Model Assumptions — TB Screening CEA

All assumptions recorded here with justification and limitation flag where relevant.

---

## 1. Model Structure

| # | Assumption | Justification | Limitation? |
|---|---|---|---|
| S1 | Closed cohort — no new entrants after t=0 | Standard for CEA of a single screened cohort | Yes — ignores future migration waves |
| S2 | Monthly cycles (1 month per cycle) | Balances computational precision with disease timescale | Minor |
| S3 | Lifetime time horizon (55 years, 660 months; entry age ~25 → age ~80; NICE reference case) | NICE DSU guidance: horizon must capture all important differences in costs and outcomes — effectively lifetime for LTBI screening. Dale et al. 2022 (AJE) uses lifetime. Age-varying ONS background mortality implemented. | Minor — ONS rates held constant within each 5-yr band; individual emigration not modelled |
| S4 | Cohort size = 100,000 | Matches decision tree denominator | — |
| S5 | Markov memorylessness — transition probabilities constant, no history dependence | Standard Markov assumption; simplifies computation | Yes — reactivation risk actually declines over time post-infection |
| S6 | Half-cycle correction applied | Standard practice; corrects for overcount in cycle-start assumption | — |
| S7 | No dynamic transmission in base case — force of infection fixed at entry; Scenario 2 applies a static secondary-case multiplier (β = 0.205; conservative base: 5 contacts/case/yr × 4.1% per-contact transmission probability) | Base case: standard convention (Jit 2011, Dale 2022, Green 2025). Scenario 2: per-contact probability anchored to Brooks-Pollock et al. 2020 UK β estimate; 5 contacts/yr chosen as conservative lower bound for UK contact investigation; DSA: 5/10/15 contacts/yr. | Yes — Scenario 2 is a static approximation, not a dynamic transmission model; ignores time-varying force of infection and herd immunity effects |
| S8 | Sequential incremental ICER analysis with extended dominance detection (NICE methods) | Implemented 28 Feb 2026: pairwise vs No Screening retained for reference + sequential ICERs computed on efficient frontier; simply and extendedly dominated strategies identified and flagged | — |
| S9 | DSA run for three comparisons: (1) Parallel Sx+QFT (Ultra) vs No Screening — top of efficient frontier; (2) Symptom screen+CXR vs No Screening — second strategy on the frontier (seq. ICER £3,822/QALY); (3) Cough+CXR (TB sx) vs No Screening — entry-level non-IGRA strategy (updated 19 Mar 2026 from Parallel Cough+QFT (Ultra), which is now extendedly dominated) | Three tornado diagrams characterise the efficient frontier strategies | Yes — extendedly dominated strategies not individually characterised |
| S10 | PSA: if sampled parameters cause row sums > 1, outgoing probabilities rescaled proportionally | Standard practice for PSA in Markov models to maintain valid probability matrices | Minor — very rare in practice with well-specified distributions |

---

## 2. Population & Epidemiology

| # | Assumption | Value | Source | Limitation? |
|---|---|---|---|---|
| E1 | LTBI prevalence at cohort entry | 17.8% | Berrocal-Almanza et al. 2022 Lancet Public Health (17.8% IGRA positivity, n=37,268, age 16-35, England 2011-2018); confirmed by UKHSA TB in England 2025 (Prevention 2024): 15.1% (6,097/40,490 tested); Pareek 2011: 20% | May vary by country of origin and year of arrival |
| E2 | Active TB prevalence at cohort entry | 1% (1,000/100,000) | Zenner et al. 2025 Eur Respir J decision tree (Initial_conditions.xlsx col prevalence = 0.01). All TP/FN values in the Excel are internally consistent with this denominator and cannot be changed independently. | Yes — higher than published pooled estimates for screened migrant populations: pooled meta-analysis 215/100,000 (95% CI 112–350; 36 studies, 40M+ migrants; Chen et al. 2025 Global Health Research and Policy); highest-risk subgroups (refugees, countries >300/100k) reach ~440–492/100,000. The Zenner decision tree value likely represents a specifically high-risk target population (recent arrivals from highest-burden countries, >500/100k TB rate). Direction of bias: if true prevalence is lower (~200–500/100k), the model overstates the absolute benefit of active TB detection at entry; relative comparisons between strategies are unaffected (TP/FN proportions are consistent). Parameter retained to maintain consistency with the decision tree denominator. |
| E3 | Uninfected proportion at entry | 81.2% | Derived (1 − 0.178 − 0.01) | — |
| E4 | No MDR-TB modelled | — | UK has low drug-resistance context; MDR-TB <2% of UK cases (UKHSA 2024) | Yes — may underestimate costs if MDR-TB proportion rises |
| E5 | Rate of progression from latent to active TB is constant (linear monthly probability) | Lu→Iu: 0.0912%/month (1.09%/yr) | Berrocal-Almanza et al. 2022 Lancet Public Health (UK IGRA+ migrants); monthly conversion 1-(1-0.0109)^(1/12). Note: this rate (1.09%/yr) is at the HIGH end of the range used in published LTBI CEA studies (0.02%–1.25%/yr; Oxlade et al. 2013 PLoS ONE SR of 13 studies). It is appropriate here as Berrocal-Almanza 2022 used a UK IGRA+ migrant cohort (higher-risk than general population); most other studies used lower general-population rates. Reactivation rate is the 2nd most influential parameter in LTBI CEA (potential influence 92–115%; Oxlade 2013). | Yes — true risk is highest in first 2 years post-infection; semi-Markov would be more accurate |
| E5b | **Two-phase reactivation (implemented 19 Mar 2026):** Phase 1 (cycles 1–60, years 0–5): 0.000912/mo (1.09%/yr) — Berrocal-Almanza et al. 2022 Lancet Public Health; short-term follow-up rate in newly arrived UK IGRA+ migrants. Phase 2 (cycles 61–660, years 5–55): 0.000333/mo (0.40%/yr) — Horsburgh 2004 NEJM PMID15496613; declining reactivation risk for remote LTBI; consistent with Dale et al. 2022 two-phase approach. Rationale: constant application of 1.09%/yr over 55 years gives ~45% lifetime reactivation, exceeding the literature consensus of 10–22% for untreated LTBI; the two-phase model corrects this to ~22% lifetime reactivation. Config parameter: p_react_phase2 = 0.000333/mo. To revert to single-phase: set p_react_phase2 = p_LatentUndiagnosed_ActiveUndiagnosed (0.000912). | Horsburgh CR Jr. NEJM. 2004;350(20):2060–7. PMID15496613. Berrocal-Almanza et al. 2022. Dale et al. 2022 AJE. | Yes — phase boundary (5 years) is approximate; transition to lower rate is abrupt rather than continuous |
| E6 | Within-cohort new LTBI acquisition: Uninfected → LatentUndiagnosed | 0.0001/month (~0.12%/yr) | Low rate reflecting within-migrant-community transmission only; most LTBI is prevalent at entry (UKHSA 2024: UK-acquired transmission ~15–20% of cases) | Yes — rate is a model assumption, not directly measured; UK-acquired LTBI is likely underestimated |
| E7 | Background mortality varies by 5-year age band across the 55-year lifetime horizon (IMPLEMENTED 10 Mar 2026; extended to 11 bands 11 Mar 2026) | ONS age bands 25-29 (49/100k/yr), 30-34 (70), 35-39 (101), 40-44 (155), 45-49 (244), 50-54 (388), 55-59 (618), 60-64 (989), 65-69 (1460), 70-74 (2430), 75-79 (4260). Previous fixed rate 0.0000417/mo (ONS 16-35) replaced. | Working-age migrant cohort, ONS National Life Tables 2020-2022 | Minor — ONS rates held constant within each 5-yr band; does not model emigration or individual-level heterogeneity. |
| E8 | Active TB state mortality rates are fixed across the 55-year horizon (not age-stratified) | Active TB-specific mortality rates (p_ActiveUndiagnosed_Dead = 0.010/mo, p_ActiveTreated_Dead = 0.0012/mo, etc.) are held constant at values calibrated to a young UK migrant cohort (Pedrazzoli 2019, age 15–44; Tiemersma 2011) and are not updated with age bands. Background mortality age-banding (E7) does not apply to active TB states. | Quantitative impact is expected to be minimal: (1) active TB mortality (0.010–0.015/mo untreated) substantially exceeds background mortality even at age 75–79 (0.00355/mo), so the relative contribution of background ageing is small; (2) most active TB events in the model occur in the first 10–20 cycles, when the cohort is young and the discrepancy is negligible; (3) discounting down-weights late-horizon events. Acknowledged as a simplifying assumption. |

---

## 3. Screening & Initial Conditions

| # | Assumption | Justification | Limitation? |
|---|---|---|---|
| SC1 | Active TB TP/FN at entry taken directly from Zenner et al. 2025 Eur Respir J decision tree (Excel col 4: TP, col 9: FN) | Excel accounts for full pathway (symptom screen → test → culture confirmation). Source paper: Zenner D, Haghparast-Bidgoli H, Chaudhry T, et al. How to diagnose TB in migrants? Eur Respir J 2025;66:2402000. Umbrella review + decision tree analysis of 34 two- and three-test combinations in 100,000 migrants; provides pooled test sensitivities/specificities from 32 systematic reviews (437 original studies). The paper explicitly does not address cost-effectiveness: *"Costs of TBI tests vary widely and will need to be considered alongside local logistical, operational and resource considerations"* — the present Markov analysis fills this gap. | — |
| SC2 | LTBI detection taken directly from Zenner et al. 2025 Eur Respir J decision tree (Excel col 5). LatentUndiagnosed = 17,800 (Berrocal-Almanza 2022: 17.8% per 100k) − Excel detected. | Excel LTBI column used directly; consistent with decision tree design. The decision tree embeds its own lower LTBI prevalence assumption (~0.1–0.15%). | The Excel-implied LTBI prevalence (~0.1–0.15%) is much lower than background prevalence (17.8% Berrocal-Almanza 2022); discrepancy likely reflects the decision tree cohort including all migrants (not only high-risk), or sequential pathway filtering — should be discussed as a limitation |
| SC3 | False positives (FP) → Uninfected state; not treated; zero treatment cost | UK practice: culture confirmation required before LTBI/TB treatment starts → FPs identified and discharged without treatment. | Yes — in LMIC settings FPs may be treated empirically; model would underestimate costs there |
| SC4 | All screening occurs at t=0 (entry to cohort) | One-off screening programme at arrival | Ignores follow-up screening or opportunistic detection over time (partially captured via Lu→Ld transition) |
| SC5 | Upfront diagnostic cost = Tot_costs / 100,000 per person | From Zenner et al. 2025 Eur Respir J decision tree (Excel col 6) | Includes all test costs for the full pathway; may not reflect real-world wastage |
| SC6 | **Parallel IGRA strategies — Option B LTBI override:** For 8 parallel IGRA strategies (parallel allsx/cough × QFT/T-SPOT × Xpert/Ultra), Excel LTBI detection values (~72–116/100k) are overridden with population-based estimates accounting for both sensitivity and specificity: LTBI detected = n_ltbi × sensitivity − (1−specificity) × n_uninfected. Using Pai et al. 2008 BCG-vaccinated pooled estimates (QFT spec 0.96; T-SPOT spec 0.93): QFT = 17,800×0.83 − 0.04×81,200 = **11,526/100k**; T-SPOT = 17,800×0.88 − 0.07×81,200 = **9,980/100k**. Active TB TP/FN retained from Excel rows unchanged. Upfront cost = Excel row cost + £500k/100k programme delivery cost (all-sx screen) or £300k/100k (cough-only screen). **IGRA uptake: base case assumes 100% of eligible migrants accept the IGRA test.** | Rationale: Excel LTBI values reflect SEQUENTIAL pathway yield (~100/100k). Parallel pathways apply IGRA to ALL migrants, so detection = prevalence × sensitivity − FP rate × uninfected. Specificity corrected from implicit 1.00 to literature values (Pai M et al. Ann Intern Med. 2008; BCG-vaccinated meta-analysis: QFT 0.96, 95% CI 0.94–0.98; T-SPOT 0.93, 95% CI 0.86–1.00). T-SPOT strategies are SIMPLY DOMINATED at base (17.8% LTBI): lower specificity → more FP → higher cost + fewer QALYs than QFT-equivalent. 100% uptake: real-world typically 60–80% (UKHSA 2024); 100% is upper-bound. | Yes — two remaining limitations: (1) Programme delivery cost (£300k–£500k/100k) is a model estimate. DSA swept £0.50–£60/pp; breakeven ≈ £1,300/pp (far above plausible range); `output/igra_programme_cost_sensitivity.csv`. (2) 100% uptake overstates LTBI detection benefit. Uptake SA (50%–100%, 5% steps): ICERs remain <£25k across full range (£2,684–£3,740/QALY for QFT at 70–100% uptake); `output/igra_uptake_sensitivity.csv`. IGRA specificity SA: `output/igra_specificity_sensitivity.csv` (sweep 0.90/0.93/0.96/1.00). |

---

## 4. LTBI Test Sensitivities — FOR REFERENCE ONLY (superseded by SC2)

> **Note:** LTBI detection at entry is taken directly from the Zenner et al. 2025 Eur Respir J decision tree (Excel col 5). These sensitivity values are **no longer used in model calculations** — they are retained here for reference and methods documentation only.
>
> **Literature context (Oxlade et al. 2013 PLoS ONE 8(3): e56044):** A systematic review of 13 CEA studies of IGRA vs TST (high-income countries; migrants, contacts, HCWs) found that IGRA vs TST **effectiveness differences are clinically negligible** — typically <1 day of QALYs gained; for recent immigrants specifically (Linas 2011), gain = **0 days**. *"Cost is the main determinant of results, not effectiveness."* The most influential parameters in LTBI diagnostic CEA are: LTBI prevalence (147%/97% potential influence), reactivation rate (115%/92%), active TB treatment cost (106%/85%). Average Drummond quality score: 72%. This review validates the model finding that IGRA-based strategies are dominated by symptom screening: incremental effectiveness differences are negligible and cost drives the ranking.

| Test | Sensitivity | Source |
|---|---|---|
| QFT-GIT | 80% | NICE NG33 (2016); Pai et al. 2014 Lancet Infect Dis |
| T-SPOT.TB | 85% | NICE NG33 (2016); Pai et al. 2014 Lancet Infect Dis |
| TST Mantoux | 75% | Lower due to BCG cross-reactivity in migrant cohort |
| TST then QFT (serial) | 60% (= 0.75 × 0.80) | Combined serial sensitivity |
| CXR, Xpert, symptom screen | 0% | Cannot detect latent TB |

---

## 5. Treatment

| # | Assumption | Value/detail | Source |
|---|---|---|---|
| T1 | LTBI treatment regimen: 3HR (3 months isoniazid + rifampicin) as base case | £64/month (£192.67 total / 3 months: drugs £79.88 + tests £4.58 + staff £108.21) | Snoad et al. 2025 medRxiv (TB_costs_UK.xlsx Treatment sheet); NICE NG33 first-line recommendation |
| T2 | Active TB treatment | Standard 2HRZE/4HR (6 months) | NICE NG33; BTS guidelines |
| T3 | No treatment for LTBI in "not treated" state | Declined or not offered | UK care cascade data |
| T4 | LTBI treatment completion → LatentCompleted (monthly rate 0.253/month; Surey et al. 2021); then LatentCompleted → Uninfected | Treatment 60–90% efficacious; modelled as monthly transition rate over ~3 months (3HR regimen); real-world completion 55.5% contacts (UKHSA 2025 Prevention 2024) | Yes — structural assumption: efficacy modelled as monthly transition; does not capture exact time profile of immunological protection |
| T5 | No LTBI re-treatment after discontinuation (re-engages at LatentDiagnosed) | Conservative assumption | — |
| T6 | Active TB treatment completion → ActiveCompleted → Uninfected (0.10/month); utility 0.92 applied in ActiveCompleted state | Post-treatment return to population norms per Jit et al. 2011 (London UK); norm value 0.92 per Kind et al. 1998 (age 16–35) | Yes — permanent post-TB sequelae not captured; a dedicated "CuredTB" absorbing state with reduced utility (0.88–0.91; Ranzani 2021 SR) would be more accurate |
| T7 | LTBI treatment completion: monthly rate 0.253/mo (Surey et al. 2021, 76%/3mo 3HR arm). LTFU rate 0.057/mo (Surey 2021, 4/25=16%/3mo → geometric monthly). Effective model completion = 0.253/(0.253+0.02+0.057+0.000045) = **76.7%** ≈ 76% trial rate (self-calibrating). UKHSA 2025 real-world 55.5% used as SA lower bound | The base case correctly matches the Surey 2021 trial completion rate (76%). Conditional completion formula: p_comp/(p_comp+p_disc+p_ltfu+p_dead) | Yes — standard Markov approximation; impact of approximation now negligible as base case self-calibrates to 76.7% |
| T7b | Active TB LTFU rate: **0.005/mo** (UKHSA TB in England 2025, treatment outcomes 2023: ~3% LTFU over 6-month course; monthly geometric rate = 1−(0.97)^(1/6) ≈ 0.005/mo). PSA: beta(3.98, 791). | UKHSA TB in England 2025, Section 4 (treatment outcomes 2023). England active TB LTFU ~3% in 2023. | Minor — small compartment; impact on ICERs negligible |
| T8 | LatentCompleted has no outflow to latent or active TB states — treatment failures not modelled. Virtually all treated LTBI patients eventually become Uninfected (~99% over 40 years) | NICE NG33: 60–90% treatment efficacy implies 10–40% residual risk. A treatment-failure pathway (LatentCompleted → LatentUndiagnosed) would be more accurate | Yes — overestimates eventual LTBI treatment efficacy. Impact small given tiny detected LTBI numbers |

---

## 6. Costs

| # | Assumption | Justification | Limitation? |
|---|---|---|---|
| C1 | NHS & PSS perspective | NICE reference case for UK health technology assessment | — |
| C2 | Cost year: **2023/24 GBP** (pounds sterling). Older costs inflated using NHS Healthcare Cost Inflation Index (HCHS, PSSRU) where available (2008–2016), extended with ONS CPIH (2016–2023). See C8 for specific methodology. | NICE reference case recommends NHS cost inflation index for inflating historical NHS costs; HCHS is the standard PSSRU source | Inflation factors apply to cost components with known base years only; components already at 2023/24 prices (Snoad 2025, NHS tariffs) are used directly |
| C3 | No productivity losses or societal costs | Narrow NHS perspective; consistent with NICE TA approach | Yes — underestimates full societal cost of TB |
| C4 | Diagnostic costs are one-time at entry, not recurring | Screening programme is single-entry | — |
| C5 | Cost of undiagnosed LTBI (LatentUndiagnosed) = £0/month | No healthcare contact if LTBI undetected | Minor — small background primary care contacts not captured |
| C6 | Cost of undiagnosed active TB (ActiveUndiagnosed) = £0/month | No NHS cost until diagnosis | Yes — in reality, undiagnosed active TB patients make unplanned healthcare contacts (emergency, primary care); a model estimate would be needed; flagged as ⚠ in config.csv |
| C7 | Cost of LTBI LTFU = £0/month | No NHS costs if completely lost to follow-up | Yes — migrants may re-present; some sporadic contacts uncosted |
| C8a | **Cost_ActiveTreated DSA range cross-check:** Base case £770/month (NICE CG117 2011 ×1.63 inflation). DSA uses gamma(shape=16, scale=48.125) → 97.5th percentile ≈ £1,161 (chi-square approximation: 48.125/2 × qchisq(0.975,32) ≈ £1,161). White & Jit 2015 cross-check (NICE NG33 Appendix I2, Table 2.5) gives £1,183/month in 2023 GBP. The existing DSA 97.5th percentile (£1,161) approaches the White & Jit upper bound (2% gap) — existing distribution is adequate. No config change required. | White & Jit 2015 NICE NG33 Appendix I2; NICE CG117 | Minor — DSA range slightly undershoots published upper bound by ~2% |
| C8 | **Inflation of historical cost components:** (a) Cost_ActiveDiagnosed: Pareek et al. 2011 (Lancet Infect Dis, PMID 21514236) CXR £39.85 + culture £14.23 are in **2010 GBP** (confirmed from paper). Inflated using HCHS (PSSRU NHS Cost Inflation Index, annual % change rates) 2010→2016 cumulative factor = ×1.186, then ONS CPIH 2016→2023 = ×1.273; combined factor = **×1.510**; inflated subtotal = £81.64. Remaining components (Xpert Ultra £38.80; staff tariffs £23.25) are at 2023/24 NHS tariff prices. Final Cost_ActiveDiagnosed = **£144** (2023/24 GBP). HCHS data sourced via `inflatably` R package (HealthEconomicsHackathon/inflatably v0.1). (b) Cost_ActiveTreated: NICE CG117 Appendix L (2011), base year **2008/09 GBP** (NHS financial year). Inflated using NHS cost inflation ×1.63 = £4,622/episode ÷ 6 months = £770/month. HCHS+CPIH combined gives ×1.55 (£733/month); ×1.63 likely reflects NHS-specific labour cost inflation running above general CPI in 2009–2016. Conservative relative to White & Jit 2015 cross-check (£1,183/month in 2023 GBP). | PSSRU HCHS annual rates 2008–2016; ONS CPIH index (base 2015=100): 2016=101.0, 2023=128.6; Pareek 2011 PMID 21514236; NICE CG117 NBK97875 | Yes (minor) — HCHS available only to 2016; CPIH used for 2016–2023 period as proxy; CPIH grows slightly slower than NHS-specific inflation, so 2016–2023 component may be marginally underestimated |

---

## 7. Health Utilities (QALYs)

| State | Value | Source |
|---|---|---|
| Uninfected / asymptomatic LTBI (all latent states except on treatment) | 0.92 | Kind et al. 1998 BMJ (UK EQ-5D population norms, age-specific value for 16–35 yr cohort). Note: overall population mean in Kind 1998 is ~0.86; the age-specific value for 16–35 yr olds is ~0.92, reflecting higher utility in younger adults. Using the overall mean (0.90) would underestimate the healthy baseline for this cohort. |
| LTBI on treatment | 0.92 | Bauer et al. 2015 Qual Life Res; Dale et al. 2022 AJE — no measurable disutility from LTBI treatment; baseline utility retained. Standard UK CEA practice (consistent with Pareek 2011, Zenner 2022). SA with annual decrement 0.0133 (Bauer 2015) applied: Parallel Sx+QFT (Ultra) ICER rises to £14,626/QALY (frontier unchanged). Dale 2022 found all Australian strategies dominated at this decrement — different outcome here due to shorter UK treatment regimen (3HR ~3 months vs 9H), higher LTBI prevalence, and favourable WTP. CSV: output/csv/ltbi_tx_utility_sensitivity.csv |
| Active TB diagnosed (pre-treatment) | 0.68 | Jit et al. 2011 London TB study; diagnosis alone does not improve QoL before treatment starts |
| Active TB undiagnosed/untreated/discontinued/LTFU | 0.68 | Jit et al. 2011 London TB study (PMC3273731); directly measured in UK population |
| Active TB on treatment | 0.81 | Jit et al. 2011 London TB study (PMC3273731) |
| Post-TB completed (ActiveCompleted) | 0.92 | Jit et al. 2011 (London UK, PMC3273731): *"We assumed cases without active tuberculosis to have the same utility score as population norms."* Population norm age 16–35: 0.92 (Kind et al. 1998 BMJ). SA range 0.80–0.95 captures sequelae uncertainty (Ranzani 2021 SR: 0.88–0.91 in other studies). |
| Dead | 0 | Convention |

Notes:
- Previous versions used Mugwagwa et al. 2021 (South Africa) for active TB utilities. Replaced with Jit et al. 2011 UK-measured values as more appropriate for UK migrant population.
- **Limitation:** Full recovery (0.92) assumed in ActiveCompleted state. Permanent post-TB sequelae (lung damage, fatigue) are not captured. A dedicated "CuredTB" absorbing state with reduced utility would better reflect lifelong sequelae (see T6).

---

## 8. Discounting

| # | Assumption | Value | Source |
|---|---|---|---|
| D1 | Annual discount rate | 3.5% | NICE reference case (applied to both costs and QALYs) |
| D2 | Monthly discount rate | (1.035)^(1/12) − 1 ≈ 0.287%/month | Converted from annual |

---

## 9. Probabilistic Sensitivity Analysis

| # | Assumption |
|---|---|
| P1 | 1,000 Monte Carlo simulations |
| P2 | Beta distributions for probabilities and utilities (bounded 0–1) |
| P3 | Gamma distributions for costs (bounded >0, right-skewed) |
| P4 | Initial conditions (from Excel/prevalence) are NOT varied in PSA — treated as fixed; uncertainty addressed in DSA |
| P5 | WTP threshold: £25,000/QALY (NICE lower bound, updated 1 Dec 2025, effective April 2026); £35,000 shown in CEAC | Source: https://www.nice.org.uk/news/articles/changes-to-nice-s-cost-effectiveness-thresholds-confirmed |
| P6 | If PSA draw causes transition matrix row sum > 1, outgoing probabilities rescaled proportionally to restore validity; diagonal set to 0 |
| P7 | 8 diagnostic test cost parameters (Cost_TST, Cost_QFT, Cost_TSPOT, Cost_CXR, Cost_Culture, Cost_CoughScreen, Cost_Xpert, Cost_XpertUltra) removed from config.csv on 8 Mar 2026 — they were never referenced in MasterTBModel.R. Screening costs are embedded in the Excel decision tree (ic$Tot_costs / excel_denom) and flow into each strategy's test_cost directly. Removing these orphaned params has zero effect on model results. |

---

## 10. Background Mortality

| # | Assumption | Value | Source |
|---|---|---|---|
| M1 | Background monthly mortality (working-age cohort) | 0.0000417/month (≈ 50/100,000/yr) | ONS UK Mortality Statistics 2024 age-specific (16–35): ~50/100k/yr |
| M2 | LTBI states (all): background mortality only | 0.0000417/month | LTBI itself does not increase mortality |
| M3 | LTBI treated: slightly elevated | **0.0000450/month** | Background (0.0000417) + minimal hepatotoxicity increment. NICE NG33: fatal hepatitis <0.001%/course → ~0.0000017/month excess. |
| M4 | LTBI LTFU: background mortality only | **0.0000417/month** | LTFU patients are not on treatment → no hepatotoxicity risk. |
| M5 | Untreated active TB mortality (undiagnosed + not treated states) | **0.010/month** (undiagnosed); **0.015/month** (not treated) | Tiemersma et al. 2011 (PLOS ONE): ~50% 5-year mortality pre-chemotherapy. SA range: **0.008–0.025/month** (38–78% 5-year). PSA: beta(3.95, 391) and beta(8.85, 581). |
| M6 | Treated active TB mortality | **0.002/month** | UKHSA 2025: 3.8% 12-month all-cause mortality all ages → 0.0032/month. Cohort is 16–35 yr old migrants; age-adjusted lower. 0.002/month ≈ 1.2% per 6-month course — plausible but no age-specific UK source. ⚠ Model estimate. SA range: 0.0006–0.0034. |
| M7 | Background mortality is age-varying | ONS National Life Tables 2020-2022: 5-year bands updated every 60 cycles. Active TB mortality states unchanged — TB-specific mortality dominates background at working ages. | Minor residual limitation: within-band rates constant, no emigration modelling. |
| M8 | ActiveCompleted state uses age-varying background mortality | ActiveCompleted mortality updates with ONS age bands so it cannot fall below background at older ages. | — |


---

## 11. Calibration

Two parameters were calibrated to UK screening programme data (values in config.csv):

| Parameter | Value | Calibration target | Result | Source |
|---|---|---|---|---|
| p_LatentUndiagnosed_LatentDiagnosed | 0.001/month | UK LTBI programme: 17.2% of eligible migrants tested in 2024 (40,490/234,814; UKHSA 2025 Prevention 2024) | Model gives 18.4% diagnosed over 20-yr horizon — slight overshoot, acceptable | UKHSA TB in England 2025 (Prevention: England 2024) |
| p_LatentDiagnosed_LatentTreated | **0.151/month** | 74.8% of diagnosed LTBI initiate treatment (UKHSA 2024 TB Prevention, England 2023). Conditional rate: p_treat/(p_treat+0.051)=0.748 → p_treat=0.151 | 74.8% → LatentTreated; 24.7% → LatentNotreated | UKHSA TB Prevention England 2023 (TB in England 2024 Report) |

Note: Both parameters slightly overshoot their calibration targets. Further calibration would improve accuracy but is not critical for current analysis given the small absolute differences.

---

## Open / Pending Assumptions

| Issue | Status |
|---|---|
| Excel LTBI column interpretation | **RESOLVED** — Excel col 5 used directly; decision tree embeds its own ~0.1–0.15% LTBI prevalence assumption |
| Time-varying rate of progression from latent to active TB | Not implemented (flat monthly observed-incidence rate used); flagged as limitation (E5, S5) |
| Cost of inaction for undetected LTBI/active TB | **RESOLVED** — captured implicitly via LatentUndiagnosed→ActiveUndiagnosed reactivation; explicit Scenario 1 (Cost of Unmanaged TB) also implemented |
| Post-TB sequelae QALY after ActiveCompleted→Uninfected transition | No dedicated CuredTB state; qaly_ActiveCompleted = 0.92 (Kind 1998, age 16–35 population norm); full recovery assumption. |
| Age-dependent mortality | **RESOLVED** — ONS 5-year age-band rates implemented. Residual limitation: within-band rates constant, no emigration. |
| LTBI treatment completion (Markov approximation) | **RESOLVED** — effective completion 76.7%, matching Surey 2021 trial rate (76%). |
| LTBI treatment 100% eventual efficacy | Documented T8 — no treatment-failure pathway in LatentCompleted; small impact |

---

## 12. Discussion Limitations

These limitation paragraphs are intended for inclusion in the Discussion section of the thesis.

### Limitation 1 — Static (time-invariant) rate of progression from latent to active TB

> The model applies a constant monthly rate of progression from latent to active TB (0.000912/month, 1.09%/yr), derived from the observed incidence of active TB in untreated IGRA-positive UK migrants (Berrocal-Almanza et al. 2022, Table 4). This observed incidence conflates true endogenous reactivation with any residual new infection, and assumes the rate remains constant over the 55-year lifetime horizon. In reality, risk is highest in the first 2 years post-infection and declines over time (Horsburgh 2004). A semi-Markov or time-inhomogeneous model would better capture this dynamic — the current approach likely overestimates late-horizon progression and underestimates early risk. The net effect on ICERs is uncertain: if early progression drives most TB cases, strategies detecting LTBI at entry would appear even more cost-effective; if late progression predominates, the benefit of LTBI screening would be attenuated.

**References:** Horsburgh CR Jr. Priorities for the treatment of latent tuberculosis infection in the United States. NEJM. 2004;350(20):2060–7. | Berrocal-Almanza LC et al. Epidemiology & Infection. 2022.

### Limitation 2 — No onward TB transmission modelled in base case (addressed in Scenario 2)

> The base case does not capture onward TB transmission beyond the screened cohort. By detecting and treating active TB and LTBI, screening prevents secondary cases in the community that are not accounted for in a closed static cohort. Base case ICERs are therefore conservative (less favourable for screening than reality). **Scenario 2 (Secondary Cases Prevention) addresses this directly:** applying a secondary case multiplier (β = 0.205; conservative base: 5 contacts/case/yr × 4.1% per-contact probability; anchored to Brooks-Pollock et al. 2020) and cost per secondary case (£6,055; Green et al. 2025). At the conservative base, Parallel Cough+QFT (Ultra) enters the efficient frontier at £2,684/QALY; DSA sweeps 5/10/15 contacts/yr. See `output/scenario_frontier_comparison.png`.

**References:** Pareek M et al. Thorax. 2011. | Aldridge RW et al. Lancet Infect Dis. 2016. | Brooks-Pollock E et al. PLOS Comput Biol. 2020;16(4):e1007795 (PMID 32218567). | Green N et al. ERJ Open Res. 2025 (PMC12183743).

### Limitation 3 — Background mortality within-band rates held constant (minor)

> The model implements age-varying background mortality using 11 ONS bands (ages 25–29 to 75–79; rates 49, 70, 101, 155, 244, 388, 618, 989, 1460, 2430, 4260 per 100,000/yr), updated every 60 cycles (5 years). Rates are held constant within each 5-year band — the true continuous increase in mortality with age is approximated by a step function. With an 11-band, 660-cycle model the step increments are small relative to the scale of mortality differences between strategies. Impact is further attenuated by discounting (discount factor at year 55 ≈ 6.6, so terminal years contribute ~15% of their nominal value). Direction of bias: marginally favours all strategies equally, so effect on ICERs is minimal.

### Limitation 4 — PSA does not resample initial state distributions or detection rates

> Initial conditions (TP, FN, LTBI detected per 100,000 from the decision tree) and effective detection rates are treated as fixed in the probabilistic sensitivity analysis. Their uncertainty is addressed in the one-way DSA (tornado diagram). A fully probabilistic treatment would resample detection rates (eff_active_* parameters, which have beta distributions in config.csv) for each PSA simulation and rebuild strategy-specific initial distributions accordingly. This would widen credible intervals, particularly for CXR+Ultra (PSA P(CE at £25k) = 64%), but is unlikely to change conclusions for the clearly cost-effective frontier strategies (Cough+CXR TB sx, Parallel Cough+QFT (Ultra), Parallel Sx+QFT (Ultra), all P(CE) ≈ 100%).

### Limitation 5 — LTBI detection negligible for sequential strategies; corrected for parallel strategies

> For the 35 sequential strategies drawn directly from the decision tree, LTBI detection numbers are negligible (~100/100,000), meaning ICERs are driven almost entirely by active TB detection (TP/FN). A sensitivity analysis varying cohort LTBI prevalence (10%/17.8%/30%) produced identical ICERs (£1,218/£3,946/£51,469) for the pre-parallel frontier — confirming that LTBI treatment benefit was structurally absent in sequential pathways where IGRA was applied only to symptom-screen positives (~600/100k). This is a direct consequence of the decision tree's sequential pathway design (ISSUE-03).
>
> **Correction for parallel IGRA strategies (Option B, 14 Mar 2026):** For the 8 added parallel strategies, LTBI detection is overridden to n_ltbi × IGRA sensitivity (14,774–15,664/100k), reflecting universal IGRA application. This restores the primary value proposition of IGRA — LTBI detection and prevention of reactivation — and is consistent with the NICE NG33 recommended pathway. The new frontier includes four parallel IGRA strategies (CE at £25k), all of which derive meaningful benefit from LTBI treatment prevention.
>
> Implication: reported ICERs for sequential strategies are driven by early active TB detection; ICERs for parallel strategies additionally reflect LTBI treatment benefit. The two pathway types are not directly comparable on LTBI detection grounds. This is consistent with Oxlade et al. 2013 (PLoS ONE 8(3): e56044): *"Cost is the main determinant of results, not effectiveness"* — cost differences between sequential and parallel pathways drive ranking, not small differences in active TB sensitivity.

**References:** Oxlade O et al. PLoS ONE. 2013;8(3):e56044. | Oxlade O, Menzies D. Int J Tuberc Lung Dis. 2007. | Zenner D et al. Eur Respir J. 2025;66:2402000 (companion paper: establishes diagnostic accuracy superiority of IGRA-containing algorithms — dOR up to 24,670 — but explicitly excludes cost-effectiveness considerations).

### Limitation 7 — Parallel IGRA upfront costs are model-estimated, not empirically derived (ISSUE-14 resolved)

> Eight parallel IGRA strategies (parallel allsx/cough with QFT/T-SPOT followed by Xpert or Ultra) were included in the analysis (added 14 Mar 2026; SC6). Their upfront cost was estimated as the corresponding sequential Excel row cost plus an programme delivery cost of £500,000/100k (all-symptom screens) or £300,000/100k (cough-only screens), representing the additional IGRA test costs for universal application. This programme delivery cost is a model estimate based on unit costs (QFT ~£30–40/test; T-SPOT ~£35–45/test; Snoad et al. 2025; NICE NG33), not directly available from the decision tree (which did not include cost data for these pathways). LTBI detection was overridden using population-based estimates (n_ltbi × IGRA sensitivity: QFT 14,774/100k; T-SPOT 15,664/100k; Option B, SC6) rather than the decision tree's sequential pathway yield (~72–116/100k).
>
> **Sensitivity:** The cost-effectiveness of parallel IGRA strategies is sensitive to this programme delivery cost assumption. If universal IGRA programme costs are higher than assumed (e.g., £600k–£800k/100k due to logistics, staffing, or test price), some parallel strategies may become dominated. A one-way DSA on the IGRA programme programme delivery cost cost would quantify this uncertainty and is recommended for future analyses.

**Reference:** Zenner D, Haghparast-Bidgoli H, Chaudhry T, Abubakar I, Cobelens F. How to diagnose TB in migrants? A systematic review of reviews and decision tree analytical modelling exercise to evaluate properties for single and combined tuberculosis screening tests. Eur Respir J. 2025;66:2402000. https://doi.org/10.1183/13993003.02000-2024

### Limitation 8 — IGRA specificity: base case uses literature values

| Limitation 8 | IGRA specificity applied from Pai et al. 2008 (BCG-vaccinated populations) | QFT specificity = 0.96 (95% CI 94–98%); T-SPOT specificity = 0.93 (95% CI 86–100%). Formula: LatentDiagnosed = n_ltbi × sens − (1−spec) × n_uninfected. Net LTBI detection: QFT = 11,526/100k; T-SPOT = 9,980/100k. SA sweeps spec=0.90/0.93/0.96/1.00 jointly. | Pai M et al. Ann Intern Med. 2008;149:177 |

### Limitation 6 — IGRA cost savings from avoided BCG false positives not modelled

> In BCG-vaccinated populations, IGRA avoids false-positive TST results, reducing unnecessary LTBI treatment costs. Oxlade et al. 2013 (PLoS ONE 8(3): e56044) note that IGRA is often **cheaper** than TST in BCG-vaccinated populations for this reason, and that 6 of 13 reviewed studies found IGRA cost-saving relative to TST. The current model uses costs embedded in the decision tree (Zenner et al. 2025), which already reflects appropriate cascade costs per pathway including FP handling. However, the model does not explicitly model the differential false-positive rate of TST vs IGRA in BCG-vaccinated subgroups, nor separately model the cost of treating false LTBI diagnoses. This is a minor limitation for a BCG-vaccinated migrant cohort (UK migrants from high-burden countries are often BCG-vaccinated), and may cause slight underestimation of IGRA strategy cost-efficiency relative to TST-based strategies.

**Reference:** Oxlade O et al. PLoS ONE. 2013;8(3):e56044.

---

## 13. Scenario Analyses

### S13.1 — Scenario: Introducing Costs for Unmanaged Active TB (base case: £0)

**Literature context:** A systematic literature review (11 Mar 2026) found that **no published UK or high-income-country TB CEA Markov model assigns a non-zero per-cycle NHS cost to undiagnosed or not-treated active TB states.** The universal convention is zero cost until treatment initiation (Jit 2011 BMJ PMC3273731; Dale 2022 AJE PMID34017976; Green 2025 ERJ PMC12183743; Söderholm 2021 PMC7954754). This scenario is therefore a **novel sensitivity check** testing the impact of this universal simplifying assumption.

| # | Assumption | Value | Source | Limitation? |
|---|---|---|---|---|
| SA1 | Base case: £0/month to ActiveUndiagnosed and ActiveNotreated | Consistent with all published UK TB CEA models | Jit 2011 PMC3273731; Dale 2022 PMID34017976; Green 2025 PMC12183743 | Yes — universal simplification; real NHS contact costs exist |
| SA2 | Scenario: £150/month to ActiveUndiagnosed | Derived from UK pre-diagnosis contact literature using TOTAL delay as denominator. Loutet et al. 2018 (PHE, England, n=22,422; PMID29923481 — NOT a Zenner paper): median TOTAL delay 2.8 months (patient 1.3mo + healthcare 14 days). Mawer et al. 2007 (UK, BJGP PMID17263928): 2–4 GP visits over 2.8 months = 0.7–1.4/month × £49 (PSSRU 2022) = £34–£69/month. Schwartzman et al. 2002 (Canada, PMID11936743): 47% have ≥1 ED visit pre-diagnosis, mean 2.2 visits over 2.8 months = 0.37/month × £180 (NHS NCC 2022/23) = £67/month. Range: £100–£135/month; £150 adopted as conservative central estimate. Methodological precedent: Miners et al. 2017 Lancet HIV (UK, PMC5614770) assigns non-zero pre-diagnosis NHS costs in a Markov model. DSA range: £0–£300/month. | Loutet 2018 PMID29923481; Mawer 2007 PMID17263928; Schwartzman 2002 PMID11936743; NHS NCC 2022/23; PSSRU 2022; Miners 2017 PMC5614770 | Moderate — range supported by UK delay data and unit costs; LOW CONFIDENCE: no UK study directly quantifies pre-diagnosis TB NHS contact costs |
| SA3 | Scenario: £275/month to ActiveNotreated | £150/month symptom-driven contacts (same derivation as ActiveUndiagnosed — known TB patient with progressing untreated disease; Loutet 2018; Mawer 2007; Schwartzman 2002) + £125/month ECM outreach (~3 weekly visits × £44/visit; Hayward/Holden et al. 2019 Lancet PMID30799062; Hanif et al. 2017 BMC Public Health PMID29141600: refusing patients mapped to ECM Level 1–2 = weekly/fortnightly visits; NICE NG33). Internal consistency: ActiveNotreated (£275) > ActiveUndiagnosed (£150) since known patients in the system receive active outreach in addition to symptom-driven contacts. LOW CONFIDENCE — no published data on contact frequency specifically for refusing patients; UKHSA records refusal rates (~1–2%) but not contact frequencies. Limitation: assumes refusing patients comply with ECM visits — may overestimate costs if patients disengage entirely. | Loutet 2018 PMID29923481; Hayward 2019 PMID30799062; Hanif 2017 PMID29141600; NICE NG33; PSSRU 2022 | Yes — low-confidence model estimate; no direct evidence; upward bias possible if refusing patients disengage from ECM |
| SA4 | Cost_LatentUndiagnosed remains £0/month in scenario | LTBI is asymptomatic; no plausible NHS contact before progression to active TB | Pareek 2011 PMC3108102 | Minor |
| SA5 | Scenario output: `output/icer_table_inaction.csv`, `output/scenario_inaction_compare.png` | Both generated by MasterTBModel.R scenario block | — | — |

**Interpretation:** Assigning costs to unmanaged active TB states increases total costs in strategies with more undetected/untreated cases (especially No Screening). ICERs for screening strategies narrow — all strategies become more cost-effective relative to the no-screening comparator. Magnitude of shift depends on the proportion remaining undiagnosed or untreated under each strategy.

### S13.2 — LTBI Prevalence Scenario (considered; not implemented)

Tested low/base/high LTBI prevalence (10%/17.8%/30%) but found ICERs identical across scenarios — a model artifact: because LTBI detection numbers from the decision tree are negligible (~100/100k), the incremental differences between strategies are driven entirely by active TB detection (TP/FN, fixed from Excel) and do not vary with LTBI prevalence. Cut from final analysis; noted as a model limitation instead.

---

### S13.3 — Scenario 2: Transmission Prevention

**Literature context:** The base case follows universal UK TB CEA convention in excluding onward transmission (Jit 2011; Dale 2022; Green 2025). This scenario quantifies the value of prevented secondary transmission using a static multiplier approach, consistent with Pareek et al. 2011 (Lancet ID) and Green et al. 2025 (ERJ Open Res).

| # | Assumption | Value | Source | Limitation? |
|---|---|---|---|---|
| TX1 | Close contacts per untreated active TB case per year | **5** (base/conservative); 10 (DSA mid); 15 (DSA high) | Conservative lower bound for UK contact investigation. Per-contact transmission probability = 0.041, anchored to Brooks-Pollock et al. 2020 UK β = 0.41 at 10 contacts (PLOS Comput Biol, PMID 32218567). Base β = 5 × 0.041 = 0.205. p_tx_per_contact held fixed; only contacts vary in DSA/PSA. | Yes — fitted to national cluster data; may not reflect migrant-specific transmission networks |
| TX2 | Cost per secondary TB case averted | **£6,055** (base); £12,110 (DSA upper) | Green et al. 2025 ERJ Open Res (PMC12183743): NHS direct costs, 2024 GBP, drug-sensitive TB. DSA range from same paper: base doubled (£12,110) as upper bound. | Yes — NHS-direct costs only; indirect/societal costs excluded (conservative) |
| TX3 | New active TB cases estimated as: person-months in reactivating latent states × reactivation rate, with two-phase rates applied | Standard Markov cohort incidence approximation (cycle length = 1 month). Phase 1 (cycles 1–60): LatentUndiagnosed, LatentDiagnosed, LatentNotreated, LatentLtfu at 0.000912/mo; LatentDiscontinued at 0.0008/mo. Phase 2 (cycles 61–660): phase 1 rates scaled by 0.000333/0.000912 = 0.365, consistent with the two-phase reactivation implemented in the Markov model (E5b). **Corrected 19 Mar 2026 (commit 879a9f4):** original code incorrectly applied 0.000912/mo to all 660 cycles, inflating No Screening new ATB from ~3,684 to 8,444 (47% vs correct 22.5% lifetime reactivation). Fix does not affect base-case or SA ICERs; Scenario 2 Parallel Sx+QFT ICER corrected from £5,972 → £6,960/QALY. | Yes — approximation; ignores time-varying risk and competing risks |
| TX4 | Active TB cases prevented vs No Screening = No_screening_atb − strategy_atb; secondary cases avoided = cases_prevented × β | Linear multiplier; ignores herd immunity threshold effects | Yes — conservative for high-screening-coverage strategies |
| TX5 | Transmission savings subtracted from total strategy cost (not added to QALYs) | Savings reflect downstream NHS costs of secondary cases averted; QALYs of secondary cases are excluded (conservative) | Yes — does not capture QALYs gained in secondary cases prevented; understates true CE of screening |
| TX6 | PSA: contacts ~ Uniform(5, 15) per untreated case/yr; implied β ~ Uniform(0.21, 0.62) via p_tx_per_contact = 0.041; cost_per_secondary held at £6,055 base across PSA simulations | Slightly wider than Brooks-Pollock 95% CrI (0.30–0.60); 1,000 simulations; results in `output/transmission_psa_summary.csv` | — |
| TX7 | Outputs: `output/icer_societal_perspective.csv`, `output/scenario_transmission_compare.png`, `output/transmission_dsa_results.csv`, `output/transmission_dsa_bar.png`, `output/transmission_psa_summary.csv`, `output/scenario_frontier_comparison.png` (three-panel comparison across all scenarios) | All generated by MasterTBModel.R and run_plots_and_tables.R | — |

**Key results (base β = 0.205, 5 contacts/case, £6,055/case) — CURRENT (commit 879a9f4):**
- No screening new ATB cases: **3,684/100k** (22.5% lifetime reactivation; consistent with Horsburgh 2004 NEJM)
- **Frontier composition UNCHANGED from base case** (same 4 strategies): No screening → Cough+CXR → Symptom screen+CXR → Parallel Sx+QFT
- Cough+CXR (TB sx): ICER essentially unchanged at **£1,179/QALY** — negligible LTBI detection, minimal transmission savings
- Symptom screen+CXR: seq. ICER **£3,822/QALY** — same as base case; modest LTBI detection (667/100k)
- Parallel Sx+QFT (Ultra): seq. ICER **£6,960/QALY** (base case £7,745; modest ~11% improvement from transmission savings)
- Parallel Cough+QFT (Ultra): **extendedly dominated in both base case and Scenario 2** (does not enter frontier)
- Sequential IGRAs (QFT/T-SPOT/TST rows): still simply dominated — too few LTBI detected (~83/100k)

**Structural implication:** The frontier composition is robust across base case and Scenario 2. Secondary transmission savings from IGRA strategies (11,526 LTBI detected) reduce the Parallel Sx+QFT ICER by ~£785/QALY but do not change strategy ordering. The LTBI detection gap between sequential strategies (~83/100k) and parallel IGRA (11,526/100k) drives this finding, which is robust across all DSA contacts values (5–15/yr; implied β 0.21–0.62).

**References:** Brooks-Pollock E, Danon L, Jombart T, Pellis L. Defining the reproduction number of COVID-19 in the UK. PLOS Comput Biol. 2020;16(4):e1007795. | Green N, et al. Latent tuberculosis infection screening of adult close contacts: a cost–utility analysis. ERJ Open Res. 2025. PMC12183743. | Pareek M, et al. Lancet Infect Dis. 2011.

---

## 14. Pre-submission Audit Findings (1 Apr 2026)

Pre-submission internal audit (1 Apr 2026): three lenses — internal consistency,
literature benchmarking, methodological framing. No structural errors or fatal
inconsistencies found. Findings below are documentation gaps and one assumption
requiring supervisor confirmation (E2/ISSUE-1).

### L1 — Cost_ActiveTreated: inflation methodology note

Base case £770/month derived from NICE CG117 (2008/09 GBP) inflated ×1.63 to
2023/24 GBP. The combined HCHS+CPIH factor for 2008–2023 is approximately ×1.55
(£733/month); the ×1.63 factor applied here likely reflects NHS-specific labour
cost inflation running above general CPIH between 2009 and 2016. The resulting
base case value (£770/month) is conservative relative to the White & Jit 2015
cross-check (NICE NG33 Appendix I2: £1,183/month in 2023 GBP). The existing DSA
97.5th percentile (≈£1,161) approaches but does not fully reach this upper bound
(2% gap). No parameter change required; note in methods: *"The active TB treatment
cost (£770/month) was derived from NICE CG117 inflated using NHS Hospital and
Community Health Services indices to 2016, then ONS CPIH to 2023/24. This may
understate NHS labour-cost inflation post-2016; the DSA upper bound (£1,161/month)
partially captures this uncertainty."*

### L2 — "Societal perspective" scope

The primary perspective labelled "societal" includes NHS costs plus savings from
prevented secondary TB and LTBI cases (£6,055 and £458 per case respectively;
total £15.73/pp at 75% IGRA uptake base case). It excludes productivity losses (~£87–175/pp, Oxlade 2013),
patient time costs, and informal care costs. Methods note to add: *"The societal
perspective presented here incorporates NHS costs and savings from secondary TB and
LTBI cases prevented (£6,055 and £458 per case respectively) but excludes
productivity losses and patient time costs, which were not available for this
migrant cohort."* See ISSUE-5.

### L3 — Sequential vs parallel LTBI detection on different evidential bases

Sequential LTBI detection (~100–150/100k) derives from the Zenner 2025 decision
tree, which embeds a cohort LTBI prevalence of ~0.1–0.15%. Parallel LTBI detection
(~7,485–8,644/100k at 75% IGRA uptake base case) is recalculated at the model's
background prevalence of 17.8% (Berrocal-Almanza 2022). The two pathway types therefore cannot be directly compared
on LTBI detection grounds. A structural SA scaling sequential IGRA strategies to
17.8% prevalence is implemented in `scripts/run_dsa.R` (output:
`output/csv/sequential_ltbi_sa.csv`); results inform Discussion limitation.
See also Limitation 5 (section 12).

### L4 — Post-TB sequelae utility: qaly_ActiveCompleted = 0.92 (full recovery)

Ranzani et al. 2021 SR reports residual utility 0.88–0.91 post-TB treatment
(sequelae: bronchiectasis, COPD, fatigue, reduced exercise capacity). The base
case assigns 0.92 (Kind 1998 population norm, age 16–35), equivalent to full
recovery. Impact: a reduction from 0.92 to 0.88 raises the Parallel
Sx+QFT ICER by approximately 4%. SA across the Ranzani 2021 range (0.88, 0.89,
0.91) implemented in `scripts/run_dsa.R` (output:
`output/csv/post_tb_sequelae_sensitivity.csv`). Conclusion: all frontier ICERs
remain within NICE range across full sequelae range. See T6.

### L5 — PSA: uniform SD for active TB effective detection rates

All 11 non-zero active TB effective detection rate parameters (eff_active_*) carry a
uniform assumed SD of 0.05 in config.csv. The Zenner 2025 decision tree does not report
confidence intervals per strategy; the SD is a proxy. PSA uncertainty bounds for
active TB detection are therefore conservative estimates only. Limitation note for
PSA methods: *"Effective active TB detection rates were drawn from Zenner et al. 2025
with a uniform assumed standard deviation of 0.05; strategy-specific empirical
uncertainty was not available from the source decision tree."*

### L6 — Cost_LatentTreated: Snoad 2025 is an unreviewed preprint

Cost_LatentTreated = £64/month (£192.67 total for 3HR; drugs £79.88 + tests £4.58 +
staff £108.21) derives from Snoad et al. 2025 medRxiv. As of April 2026 this paper
has not been peer-reviewed. Methods note: *"LTBI treatment cost was taken from Snoad
et al. (2025, medRxiv), which was available as a preprint at the time of analysis;
the value was cross-checked against NICE NG33 drug costs and found consistent."*
Re-check on publication.

### L7 — Active TB prevalence 1% vs pooled literature (ISSUE-1)

Base case: 1% (1,000/100k) from Zenner 2025 ERJ decision tree. Pooled literature
(Chen et al. 2025 Global Health Research and Policy: 36 studies, 40M+ migrants): 0.215%;
highest-risk subgroups (refugees, countries >300/100k): ~0.44–0.49%. The 1% likely
reflects the specifically high-risk target population in the Zenner decision tree
(recent arrivals from countries with TB incidence >500/100k). Methods sentence to
add once confirmed with supervisor: *"Active TB prevalence (1%) reflects the
high-risk subgroup definition in Zenner et al. 2025 (countries with TB incidence
>150/100k); this exceeds pooled estimates across all screened migrant groups
(0.215%, Chen et al. 2025). Sensitivity analyses at 0.215% and 0.44% are reported
in supplementary table 4."* See ISSUE-1.
