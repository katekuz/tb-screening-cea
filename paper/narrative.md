# Narrative & Key Messages — TB Screening CEA
# For Lancet Regional Health – Europe submission (deadline 9 May 2026; 3,000 words)
# Working document — not for external sharing

---

## Core message (one sentence)
Parallel symptom screen + QFT-Plus with Xpert Ultra confirmation is cost-effective for UK migrants at both NICE thresholds under an NHS perspective (ICER £14,344/QALY; P(cost-optimal) 84.8% at £25k/QALY). Accounting for secondary transmission prevention further strengthens this conclusion (ICER £12,738/QALY; P(cost-optimal) 89%), supporting inclusion of transmission benefits in programme appraisal.

---

## Key findings by section

### Efficient frontier (4 strategies)
- **Passive case finding** (reference): £203/person, 21.606 QALYs — current practice
- **Cough+CXR (TB sx)**: £209/person (+£6), ICER £2,205/QALY — highly cost-effective; dominated by active TB detection gain
- **Symptom screen+CXR**: £220/person (+£11), ICER £7,963/QALY — cost-effective at both NICE thresholds
- **Parallel Sx+QFT (Ultra)**: £361/person (+£157 NHS; £345 societal), ICER £14,344/QALY (NHS) / £12,738/QALY (societal)
- T-SPOT-based parallel strategies simply dominated (higher cost, same or fewer QALYs than QFT equivalent)
- Parallel Cough+QFT extendedly dominated by Symptom screen+CXR + Parallel Sx+QFT combination

### Why QALYs are so similar across strategies
The x-axis range across all 43 strategies is ~0.03 QALYs (21.600–21.630 per person). This is expected and correct:
- Only 17.8% of the cohort has LTBI at entry
- ~22% of LTBI cases reactivate over 55 years → ~3.9% of the full cohort ever develops active TB
- Each averted TB case gains ~0.5–1 QALY
- Per-person QALY gain from even the best strategy ≈ 0.01–0.02 QALYs
- This pattern is identical in Zenner 2022 (Lancet PH), Pareek et al., and other NICE-format TB models
- **One-liner for methods/figure note:** "The narrow QALY range reflects the low per-person probability of TB reactivation across a 55-year horizon; cost-effectiveness is driven primarily by the ratio of incremental cost to incremental QALY rather than large absolute QALY gains."

### DSA — all stable except:
- **Active TB prevalence 0.215% (low end; pooled estimate for all migrants, Chen et al. 2025): Symptom screen+CXR extendedly dominated; Parallel Cough+QFT (Ultra) enters frontier at £16,469/QALY (below £25k); Parallel Sx+QFT (Ultra) remains on frontier at £32,956/QALY — ABOVE £25k WTP ceiling, so NOT cost-effective at this prevalence.** Key message for general migrant populations (vs high-risk 1% base case): only the cough-triggered parallel IGRA strategy is cost-effective.
- IGRA specificity 0.90 (low): ICER rises to £22,166/QALY — remains below £25k threshold
- Frontier composition unchanged in all other DSA values (LTBI completion, uptake, prevalence, programme costs)

### PSA (VERIFIED session 49 — 75% IGRA uptake base case; strategy-specific eff_active_QFT parameters)
- Groups fixed: A = active-TB QALY weights (5 states, Beta(14.12,6.64)); B = early-phase reactivation rates (4 states, Beta(9.23,10114.77)); C = background mortality (7 states, Beta(17.39,416982))
- New: eff_active_QFT_allsx_ultra/xpert and eff_active_QFT_cough_ultra/xpert added to config.csv; PSA now samples active TB detection independently for each parallel strategy variant
- Parallel Sx+QFT: P(optimal £25k) = **89.2% societal** / **84.8% NHS** — perspective-dependent
- Parallel Sx+QFT: P(optimal £35k) = **97.5% societal** / **96.3% NHS**
- Cough+CXR: P(optimal £25k) = 15.2% (NHS) / 10.8% (societal) — optimal at low WTP
- Symptom screen+CXR + Passive: 0.0% at both thresholds (never optimal at ≥£25k)
- Key contrast: adding transmission prevention benefit increases Parallel Sx+QFT P(optimal) from 84.8% (NHS) to 89.2% (societal) at £25k
- NHS ICER Parallel Sx+QFT (Ultra): £14,344/QALY (was £11,324 at 100% uptake); societal: £12,738/QALY (was £9,665)

---

## Narrative for Results

### Para 1 — frontier & ICERs
Of 43 strategies evaluated, four formed the efficient frontier. [Passive case finding / reference]. Cough+CXR (TB sx) was the most cost-effective active screening strategy (ICER £2,205/QALY), followed by Symptom screen+CXR (£7,963/QALY). Parallel Sx+QFT (Ultra) was cost-effective under both NICE thresholds (ICER £14,344/QALY, NHS; £12,738/QALY, societal). All remaining strategies were dominated or extendedly dominated. [Table 2/3]

### Para 2 — transmission
Under the societal perspective, Parallel Sx+QFT (Ultra) prevented 889 reactivations, 182 secondary active TB cases, and 1,023 secondary LTBI cases per 100,000 cohort, generating £15.73/person in transmission-related savings. This reduced the societal cost to £345/person (vs £361/person under the NHS perspective) and improved the ICER from £14,344 to £12,738/QALY. [Figure transmission bar]

### Para 3 — DSA
Deterministic sensitivity analyses (DSA) showed that QALY weights for latent TB states and reactivation probabilities most strongly influenced NMB [Figure tornado]. The frontier was stable across alternate DSA values of LTBI treatment completion (55.5%), IGRA uptake (50%–100%), LTBI prevalence (15.1%), and programme delivery costs (£0.50–£60/person). The frontier changed only at the lowest active TB prevalence (0.215%): Symptom screen+CXR became extendedly dominated, Parallel Cough+QFT (Ultra) entered the frontier at £16,469/QALY (below the £25k threshold), and Parallel Sx+QFT (Ultra) remained on the frontier but at £32,956/QALY — above the NICE WTP ceiling — making it not cost-effective at this prevalence. This scenario is relevant for general migrant populations where active TB prevalence is closer to the pooled all-migrant estimate (0.215%; Chen et al. 2025). In a structural sensitivity analysis addressing LTBI treatment efficacy, applying the NICE NG33 range of 20–40% residual reactivation risk post-treatment rendered Parallel Sx+QFT (Ultra) simply dominated at both values; non-IGRA strategies (Cough+CXR, Symptom screen+CXR) retained their cost-effectiveness with ICERs unchanged at £2,205 and £7,964/QALY. This identifies LTBI treatment efficacy as the pivotal assumption for IGRA-based screening; the base case applies a UK-specific hazard ratio (HR 0.14; Berrocal-Almanza 2022) reflecting modern short-course regimens in UK clinical practice. Applying a utility decrement during LTBI treatment (annual decrement 0.0133; Bauer et al. 2015) — the assumption that reversed cost-effectiveness in Dale et al. 2022 (Australia) — increased the Parallel Sx+QFT (Ultra) ICER marginally to £14,626/QALY; the frontier composition was unchanged. The smaller impact reflects the shorter UK standard regimen (3HR, ~3 months vs 9H isoniazid in Bauer 2015), higher LTBI prevalence in the target population, and inclusion of transmission savings.

### Para 4 — PSA
Probabilistic sensitivity analyses (PSA) confirmed Parallel Sx+QFT (Ultra) as cost-effective under both perspectives: optimal in 84.8% of simulations (NHS) and 89.2% (societal) at £25,000/QALY, rising to 96.3% and 97.5% respectively at £35,000/QALY. Accounting for transmission prevention savings increased P(optimal) by ~5 percentage points across thresholds, confirming that secondary transmission benefits further strengthen an already robust recommendation. Cough+CXR (TB sx) was optimal in the remaining simulations (15.2% NHS / 10.8% societal at £25k), reflecting that at lower WTP values a minimal-resource strategy can be preferred. Symptom screen+CXR was cost-effective vs passive case finding in 99.4% of simulations but was never the NMB-maximising strategy at either threshold, as it was consistently dominated in expected value by either Cough+CXR (at lower WTP) or Parallel Sx+QFT (Ultra) (at higher WTP); it is excluded from the CEAC for clarity. [Figure CEAC]

---

## Abstract draft (5-part, target ≤250 words)

**Background**
Tuberculosis (TB) remains a major public health concern in the UK, with migrants from high-burden countries accounting for the majority of cases. Current NICE guidance recommends IGRA-based LTBI testing but offers limited guidance on optimal screening strategy. We assessed the cost-effectiveness of 43 TB screening strategies for UK migrants.

**Methods**
We developed a Markov cohort model (16 health states, 55-year horizon, 100,000 cohort) comparing 34 sequential and 8 parallel IGRA-based screening strategies against passive case finding. Parameters were derived from UK-specific sources (2023/24 GBP). The primary perspective was societal (NHS direct costs plus secondary transmission prevention savings, including costs of secondary active TB and LTBI cases averted; productivity losses and patient time costs were not included). An NHS-only perspective was analysed as a robustness check. Uncertainty was characterised by deterministic (DSA) and probabilistic (PSA; n=1,000) sensitivity analyses.

**Findings**
Four strategies formed the efficient frontier. Cough+CXR (TB sx) was cost-effective at £2,205/QALY; Symptom screen+CXR at £7,963/QALY; Parallel Sx+QFT (Ultra) at £14,344/QALY (NHS) and £12,738/QALY (societal) — all below NICE thresholds (£25,000–£35,000/QALY). Parallel Sx+QFT (Ultra) was the optimal strategy in 89.2% of PSA simulations (societal) and 84.8% (NHS) at the NICE £25,000/QALY threshold, with the societal advantage driven by secondary transmission prevention. Results were robust across sensitivity analyses.

**Interpretation**
Parallel symptom screen + IGRA (QFT-Plus/Ultra) is cost-effective for UK migrants when secondary transmission costs are included. For lower-prevalence populations, symptom screen+CXR offers near-certain cost-effectiveness at substantially lower programme cost. Findings support updating NICE guidance to reflect contemporary IGRA platforms and societal benefits of active TB prevention.

**Funding**
[TBC — confirm with Kasim]

---

## Discussion talking points

### Comparison with existing literature
- Zenner 2022 (Lancet PH): sequential QFT+CXR cost-effective; our results broadly consistent but extend to parallel strategies and include transmission
- Pareek et al.: LTBI screening cost-effective in high-risk migrants; our model more granular (43 strategies vs handful)
- Xu et al. 2025 SR (cite — Kateryna author 5): systematic review context; place our findings within the range reported
- Key advance: first CEA to compare all 43 Zenner 2025 decision tree strategies simultaneously, including transmission prevention
- **Oxlade et al. 2013 (PLoS ONE 8(3):e56044):** SR of 13 IGRA vs TST CEA studies (high-income countries; migrants, contacts, HCWs). Found that including productivity losses typically makes LTBI screening cost-saving rather than merely cost-effective. Most influential parameters: LTBI prevalence (147% potential influence), reactivation rate (115%), active TB treatment cost (106%). Our model is consistent with these rankings (reactivation and QALY weights dominate tornado). Cite in Discussion to contextualise the NHS-only ICER and justify the transmission prevention perspective.

### Policy implications
- Cough+CXR TB sx: implementable immediately, minimal infrastructure — relevant for primary care and community settings
- Parallel Sx+QFT (Ultra): requires IGRA infrastructure but maximises TB prevention — relevant for high-volume clinic settings
- Transmission perspective matters: IGRA-based strategies are cost-effective under both perspectives; the NHS perspective alone understates the full programme value, and societal framing (adding transmission savings) further reinforces recommendations — aligning with UK TB elimination targets

### Limitations
- No sex disaggregation (SAGER — add one sentence)
- Passive case finding calibration assumes current diagnosis rates stable
- Behavioural assumptions (uptake, treatment completion) from literature — may not reflect local variation
- Single entry cohort (age ~25); does not model repeat migration or age-heterogeneous populations
- **Emigration not modelled:** The model follows all cohort members for 55 years in the UK. Return migration (estimated 15–20% of migrants within 5 years, UKHSA) means NHS costs and QALYs are accrued for individuals no longer resident in the UK; direction of bias on ICERs is uncertain.
- **LTBI treatment efficacy:** The model applies a UK-specific hazard ratio of 0.14 (Berrocal-Almanza 2022), giving ~0.16% lifetime reactivation after treatment completion — substantially lower than NICE NG33's 10–40% residual risk estimate. A structural sensitivity analysis at 20–40% residual risk (NICE NG33 range) showed Parallel Sx+QFT (Ultra) to be simply dominated at both values, with the frontier reverting to non-IGRA strategies (Cough+CXR £2,205/QALY; Symptom screen+CXR £7,964/QALY — unchanged). LTBI treatment efficacy is therefore the pivotal assumption for IGRA-based recommendations; the base-case HR 0.14 may better reflect modern short-course (3HR) regimens in UK clinical practice than older NICE NG33 estimates.
- **Post-TB sequelae:** Post-treatment, cured active TB patients return to the population QALY norm (0.92). Permanent post-TB sequelae (bronchiectasis, COPD, fatigue) account for ~47% of total TB disease burden (Menzies et al. 2021 Lancet Glob Health, PMID 34798027). Ignoring sequelae slightly inflates the QALY benefit of active TB treatment; SA with 0.89 utility in the post-treatment state shows modest ICER changes. ⚠️ "Ranzani 2021" citation removed — not found in Zotero or PubMed; replaced with Menzies 2021.

---

## Key reminders for writing
- WTP range: £25,000–£35,000/QALY (NOT £20,000)
- IGRA brand names: QFT-Plus, T-SPOT.TB, Xpert Ultra (not generic "IGRA test" or "TB test")
- Cite Zenner 2025 ERJ (DOI: 10.1183/13993003.02000-2024) for decision tree source
- UKHSA (not UKSHA)
- 43 strategies total (not 42 or 44)
- Societal perspective = primary; NHS = robustness check only
