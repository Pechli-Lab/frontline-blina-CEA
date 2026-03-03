#' Economic Evaluation: Blinatumomab vs Best Supportive Care
#'
#' @description
#' Compares blinatumomab + chemotherapy vs chemotherapy only under different
#' relapse pathway assumptions (1st-line vs 2nd-line usage). Script:
#'   1. Selects model output directories based on logical flags
#'   2. Loads aggregated total costs (tc) and effectiveness (te) RDS files
#'   3. Applies inflation adjustment
#'   4. Removes cost outliers using 1st-99th percentile filtering
#'   5. Calculates cost-effectiveness metrics (ICER, net benefit)
#'   6. Creates PSA object and generates:
#'        - Incremental cost-effectiveness scatter plot
#'        - Cost-effectiveness acceptability curve (CEAC)
#' @note
#' - Need to have run the blina and bsc models first to generate the results
#' - `output.path` needs to be defined
#'
#' @inputs
#'   - RDS files (per chosen scenario):
#'       * tc.i_new_blin.rds
#'       * tc.i_new_bsc.rds
#'       * te.i_new_blin.rds
#'       * te.i_new_bsc.rds
#'
#' @outputs
#'   - ce_summary.csv              : Summary table from calculate_ce_summary()
#'   - psa_blina_summary.txt       : Summary() of PSA object (cost & effect distributions)
#'   - ce-plot.pdf                 : Incremental cost-effectiveness plot
#'   - cea_table.txt               : ICER table
#'   - blina_ceac.csv              : CEAC data (% probability cost-effective by WTP)
#'   - ceac-plot.pdf               : CEAC plot (probability vs WTP)
#'

# Total discounted and undiscounted costs difference between bsc and blin ------

cat("\nCalculating cost-effectiveness results for blina vs bsc...\n")

# Specify where to get the output files from the microsimulation analyses
output_path <- file.path(getwd(), path.out2)

# Save path in the report/ folder
# - get the last folder name from path.out2 to create a save path
FOLDER <- basename(path.out2)
SAVE_PATH <- file.path(getwd(), "report", FOLDER)

if (!dir.exists(SAVE_PATH)) {
  dir.create(SAVE_PATH, recursive = TRUE)
}

bsc_path <- paste0(output_path, "/bsc/")
blin_path <- paste0(output_path, "/blin/")

# Base case analysis: discounted and undiscounted ------------------------------

res1_base <- process_cea(
  root = output_path,
  tc_blin = "tc.i_new_blin.rds",
  tc_bsc = "tc.i_new_bsc.rds",
  te_blin = "te.i_new_blin.rds",
  te_bsc = "te.i_new_bsc.rds"
)

res2_undiscounted <- process_cea(
  root = output_path,
  tc_blin = "tc.i.ud_new_blin.rds",
  tc_bsc = "tc.i.ud_new_bsc.rds",
  te_blin = "te.i.ud_new_blin.rds",
  te_bsc = "te.i.ud_new_bsc.rds",
  discounted = FALSE
)

# Scenario analysis ------------------------------------------------------------

# 1) +/- 10% on blina drug cost
res3_scenario_low <- process_cea(
  root = output_path,
  tc_blin = "tc.i_new_blin_m10.rds",
  tc_bsc = "tc.i_new_bsc_m10.rds",
  te_blin = "te.i_new_blin.rds", # effect unchanged since only cost changes
  te_bsc = "te.i_new_bsc.rds",   # effect unchanged since only cost changes
  scenario = "blina_cost_low"
)

res4_scenario_high <- process_cea(
  root = output_path,
  tc_blin = "tc.i_new_blin_p10.rds",
  tc_bsc = "tc.i_new_bsc_p10.rds",
  te_blin = "te.i_new_blin.rds", # effect unchanged since only cost changes
  te_bsc = "te.i_new_bsc.rds",   # effect unchanged since only cost changes
  scenario = "blina_cost_high"
)

# 2) IVIG scenario analysis
res5_scenario_ivig <- process_cea(
  root = output_path,
  tc_blin = "tc.i_new_blin_ivig.rds",
  tc_bsc = "tc.i_new_bsc_ivig.rds",
  te_blin = "te.i_ivig_new_blin.rds",
  te_bsc = "te.i_ivig_new_bsc.rds",
  scenario = "ivig"
)

# Combine all results ----------------------------------------------------------

ec_table_out <- rbind(
  res1_base$ec_table,
  res2_undiscounted$ec_table,
  res3_scenario_low$ec_table,
  res4_scenario_high$ec_table,
  res5_scenario_ivig$ec_table
)

psa_blina_out <- rbind(
  res1_base$psa_df,
  res2_undiscounted$psa_df,
  res3_scenario_low$psa_df,
  res4_scenario_high$psa_df,
  res5_scenario_ivig$psa_df
)

cea_table_out <- rbind(
  res1_base$cea_tbl,
  res2_undiscounted$cea_tbl,
  res3_scenario_low$cea_tbl,
  res4_scenario_high$cea_tbl,
  res5_scenario_ivig$cea_tbl
)

ceac_df_out <- rbind(
  res1_base$ceac_df,
  res2_undiscounted$ceac_df,
  res3_scenario_low$ceac_df,
  res4_scenario_high$ceac_df,
  res5_scenario_ivig$ceac_df
)


# Save the combined results to CSV files
write.csv(ec_table_out, file = paste0(output_path, "/ce_summary.csv"), row.names = FALSE, na = "")
write.csv(psa_blina_out, file = paste0(output_path, "/psa_blina_summary.csv"), row.names = FALSE, na = "")
write.csv(cea_table_out, file = paste0(output_path, "/cea_table.csv"), row.names = FALSE, na = "")
write.csv(ceac_df_out, file = paste0(output_path, "/ceac_data.csv"), row.names = FALSE, na = "")


################################################################################
###
### Bring all results together for reporting
###
################################################################################

# Life expectancy  -------------------------------------------------------------
# - conditional le from time of dx

le <- calc_le_metrics(bsc_path = bsc_path, blin_path = blin_path, rds = "les.rds", name = "le_i", out_var = "le", description = "Life Expectancy (years)")
le_bmt <- calc_le_metrics(bsc_path = bsc_path, blin_path = blin_path, rds = "lesBMT.rds", name = "le_iBMT", out_var = "le_bmt", description = "Life Expectancy post-BMT (years)")
le_rel <- calc_le_metrics(bsc_path = bsc_path, blin_path = blin_path, rds = "lesREL.rds", name = "le_iREL", out_var = "le_rel", description = "Life Expectancy relapse (years)")

le_summary <- dplyr::select(le, arm, le) |>
  left_join(le_bmt %>% dplyr::select(arm, le_bmt), by = "arm") |>
  left_join(le_rel %>% dplyr::select(arm, le_rel), by = "arm")

write.csv(
  le_summary,
  file = paste0(SAVE_PATH, "/life_expectancy_summary.csv"),
  row.names = FALSE,
  na = ""
)

# OS and PFS estimates ---------------------------------------------------------

# Define the time points (in years) at which to extract survival estimates
years <- c(1:10)

# Initialize an empty data frame to store survival estimates
survival_df <- data.frame()

# Loop over overall survival (os) and progression-free survival (pfs)
for (surv in c("os", "pfs")) {
  # Load the combined survival data
  survs <- read.csv(paste0(output_path, "/", surv, "_combined.csv")) |>
    dplyr::filter((time / 30 / 12) %in% years) |>
    dplyr::mutate(
      years = time / 30 / 12,
      arm = ifelse(arm != "BSC", "Blina", arm),
    ) |>
    dplyr::select(
      est = surv_median,
      lcl = surv_lcl,
      ucl = surv_ucl,
      years,
      arm,
      outcome
    ) |>
    dplyr::arrange(years, arm) |>
    mutate(
      est_chr = paste0(
        round(est * 100, 2),
        " (", round(lcl * 100, 2), ", ", round(ucl * 100, 2), ")"
      ),
    )

  # Add to overall dataframe
  survival_df <- rbind(survival_df, survs)
}

# surv_df <- survival_df |>
#   select(est_chr, years, arm, outcome) |>
#   arrange(desc(outcome), years, arm) |>
#   mutate(years = paste0("y", years)) |>
#   pivot_wider(
#     names_from = c(outcome, years),
#     values_from = est_chr
#   )

write.csv(
  survival_df,
  file = paste0(SAVE_PATH, "/survival_estimates.csv"),
  row.names = FALSE,
  na = ""
)

# Event probabilities ----------------------------------------------------------

blin_events <- readRDS(paste0(output_path, "/blin/event_probs.rds"))[, -1]
bsc_events <- readRDS(paste0(output_path, "/bsc/event_probs.rds"))[, -1]
event_ratio <- blin_events / bsc_events

# Get the median and 95% CI for each column
blin_events <- apply(blin_events, 2, function(x) quantile(x, c(0.5, 0.025, 0.975))) |>
  t() |>
  as.data.frame()
colnames(blin_events) <- c("median", "lcl", "ucl")
blin_events$arm <- "Blina"
blin_events$est_chr <- paste0(
  round(blin_events$median * 100, 2),
  " (", round(blin_events$lcl * 100, 2), ", ", round(blin_events$ucl * 100, 2), ")"
)

bsc_events <- apply(bsc_events, 2, function(x) quantile(x, c(0.5, 0.025, 0.975))) |>
  t() |>
  as.data.frame()
colnames(bsc_events) <- c("median", "lcl", "ucl")
bsc_events$arm <- "BSC"
bsc_events$est_chr <- paste0(
  round(bsc_events$median * 100, 2),
  " (", round(bsc_events$lcl * 100, 2), ", ", round(bsc_events$ucl * 100, 2), ")"
)

event_ratio <- apply(event_ratio, 2, function(x) quantile(x, c(0.5, 0.025, 0.975))) |>
  t() |>
  as.data.frame()
colnames(event_ratio) <- c("median", "lcl", "ucl")
event_ratio$arm <- "Ratio (Blina / BSC)"
event_ratio$est_chr <- paste0(
  round(event_ratio$median, 2),
  " (", round(event_ratio$lcl, 2), ", ", round(event_ratio$ucl, 2), ")"
)

df_event <- bind_cols(
  blin_events %>% select(est_chr) %>% rename(blina = est_chr),
  bsc_events %>% select(est_chr) %>% rename(bsc = est_chr),
  event_ratio %>% select(est_chr) %>% rename(risk_ratio = est_chr)
)
df_event$time <- rownames(df_event)
rownames(df_event) <- NULL

write.csv(
  df_event,
  file = paste0(SAVE_PATH, "/event_probabilities.csv"),
  row.names = FALSE,
  na = ""
)


# Willingness to pay threshold -------------------------------------------------

WTP_THRESHOLD <- 50000

wtp_prop <- read.csv(paste0(output_path, "/ceac_data.csv")) |>
  dplyr::filter(
    WTP == WTP_THRESHOLD,
    str_detect(Strategy, "Blinatumomab"),
  ) |>
  select(WTP, prop = Proportion, discount, scenario)

ceac_out <- read.csv(paste0(output_path, "/ce_summary.csv")) |>
    left_join(wtp_prop, by = c("discount", "scenario"))

write.csv(
  ceac_out,
  file = paste0(SAVE_PATH, "/ce_summary.csv"),
  row.names = FALSE,
  na = ""
)

cat("\nOutput for cost-effectiveness analysis saved to:\n\t", output_path, "\n\t", SAVE_PATH, "\n\n")
