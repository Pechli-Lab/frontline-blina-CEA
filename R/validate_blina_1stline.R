#' -----------------------------------------------------------------------------
#' Perform survival analysis across simulated runs
#'
#' Description:
#'   This script validates the survival analysis for Blinatumomab in first-line
#'   treatment settings. It simulates survival for Overall Survival
#'   (OS) and Progression-Free Survival (PFS) for two treatment arms:
#'   - Best Supportive Care (BSC)
#'   - Blinatumomab + BSC
#'
#' Dependencies:
#'   • `calc_sim_surv_probs()`: function to simulate survival curves
#'   • n_sim: number of simulation iterations
#'   • times: vector of follow-up timepoints (days)
#'
#' Outputs:
#'   Files written under `model-out/<scenario>/` (via `path.out2`):
#'   - os_combined.csv
#'   - pfs_combined.csv
#'   - os-plot.pdf
#'   - os-plot-5yr.pdf
#'   - pfs-plot.pdf
#'   - pfs-plot-5yr.pdf
#'
#' Notes:
#'   Ensure that the `calc_sim_surv_probs` function is defined and
#'   that `n_sim` and `times` are set in the environment.
#' -----------------------------------------------------------------------------

# Source the setup and model setup scripts
source(paste0(getwd(), "/R/functions/calc_sim_surv_probs.R"))

# Overall survival -------------------------------------------------------------

outcome <- "DEAD"

os_bsc <- calc_sim_surv_probs(
  events = outcome,
  n_sim = n_sim,
  output.path = paste0(getwd(), path.out2, "/bsc/"),
  times = times
) |>
  mutate(arm = "BSC")

os_blin <- calc_sim_surv_probs(
  events = outcome,
  n_sim = n_sim,
  output.path = paste0(getwd(), path.out2, "/blin/"),
  times = times
) |>
  mutate(arm = "Blina + BSC")

# Combine results
os_combined <- bind_rows(os_bsc, os_blin) |>
  mutate(
    outcome = "OS",
    n_sim = n_sim
  )

# Save the combined results
write.csv(os_combined, file = paste0(getwd(), path.out2, "/os_combined.csv"), row.names = FALSE)

# Plot the survival curves
plot_os <- ggplot(os_combined, aes(x = time / 365.25, color = arm)) +
  geom_ribbon(aes(ymin = surv_lcl, ymax = surv_ucl, fill = arm), alpha = 0.2, color = NA) +
  geom_step(aes(y = surv_median)) +
  labs(
    title = "Overall Survival Curves",
    x = "Time (years)",
    y = "Survival Probability"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) +
  scale_color_manual(values = c("BSC" = "blue", "Blina + BSC" = "red")) +
  scale_fill_manual(values = c("BSC" = "blue", "Blina + BSC" = "red")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
  )

# Save the plot
ggsave(
  plot = plot_os,
  filename = paste0(getwd(), path.out2, "/os-plot.pdf"),
  width = 8,
  height = 6,
  dpi = 300,
  units = "in"
)

ggsave(
  plot = plot_os + scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 0.5)) + scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, by = 0.05)),
  filename = paste0(getwd(), path.out2, "/os-plot-5yr.pdf"),
  width = 8,
  height = 6,
  dpi = 300,
  units = "in"
)

# Progression-free survival ----------------------------------------------------

outcome <- c("DEAD", "REL1", "BMT1_REL")

pfs_bsc <- calc_sim_surv_probs(
  events = outcome,
  n_sim = n_sim,
  output.path = paste0(getwd(), path.out2, "/bsc/"),
  times = times
) |>
  mutate(arm = "BSC")

pfs_blin <- calc_sim_surv_probs(
  events = outcome,
  n_sim = n_sim,
  output.path = paste0(getwd(), path.out2, "/blin/"),
  times = times
) |>
  mutate(arm = "Blina + BSC")

# Combine results
pfs_combined <- bind_rows(pfs_bsc, pfs_blin) |>
  mutate(
    outcome = "PFS",
    n_sim = n_sim
  )

# Save the combined results
write.csv(pfs_combined, file = paste0(getwd(), path.out2, "/pfs_combined.csv"), row.names = FALSE)

# Plot the progression-free survival curves
plot_pfs <- ggplot(pfs_combined, aes(x = time / 365.25, color = arm)) +
  geom_ribbon(aes(ymin = surv_lcl, ymax = surv_ucl, fill = arm), alpha = 0.2, color = NA) +
  geom_step(aes(y = surv_median)) +
  labs(
    title = "Progression-Free Survival Curves",
    x = "Time (years)",
    y = "Survival Probability",
    color = "Arm",
    fill = "Arm"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) +
  scale_color_manual(values = c("BSC" = "blue", "Blina + BSC" = "red")) +
  scale_fill_manual(values = c("BSC" = "blue", "Blina + BSC" = "red")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
  )

# Save the plot
ggsave(
  plot = plot_pfs,
  filename = paste0(getwd(), path.out2, "/pfs-plot.pdf"),
  width = 8,
  height = 6,
  dpi = 300,
  units = "in"
)

ggsave(
  plot = plot_pfs + scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 0.5)) + scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, by = 0.05)),
  filename = paste0(getwd(), path.out2, "/pfs-plot-5yr.pdf"),
  width = 8,
  height = 6,
  dpi = 300,
  units = "in"
)
