#' Process Cost-Effectiveness Analysis for Discounted/Undiscounted Results
#'
#' @param root String. Root directory path for input/output files
#' @param tc_blin file name for total cost (blin)
#' @param tc_bsc file name for total cost (bsc)
#' @param te_blin file name for total effect (blin)
#' @param te_bsc file name for total effect (bsc)
#' @param discounted Logical. If TRUE, use discounted results; if FALSE, undiscounted
#' @param scenario String or NULL. Scenario name for labeling outputs
#' @param inflation_factor Numeric. Factor to inflate costs (default 1.21 for 2018->2024)
#'
process_cea <- function(
    root,
    tc_blin,
    tc_bsc,
    te_blin,
    te_bsc,
    discounted = TRUE,
    scenario = NULL,
    inflation_factor = 1.21) {
  # Determine discount string and filter files
  discount_str <- ifelse(discounted, "discounted", "undiscounted")

  # Load data from RDS files ---------------------------------------------------
  tc_hat_all_blin <- readRDS(paste0(root, "/blin/", tc_blin))
  tc_hat_all_bsc <- readRDS(paste0(root, "/bsc/", tc_bsc))
  te_hat_all_blin <- readRDS(paste0(root, "/blin/", te_blin))
  te_hat_all_bsc <- readRDS(paste0(root, "/bsc/", te_bsc))

  # Inflate costs (model output is in 2018 CAD) --------------------------------
  tc_hat_all_blin$mean <- tc_hat_all_blin$mean * inflation_factor
  tc_hat_all_bsc$mean <- tc_hat_all_bsc$mean * inflation_factor

  # Prepare data frames
  tcs <- data.frame(cbind(tc_hat_all_blin$mean, tc_hat_all_bsc$mean))
  tes <- data.frame(cbind(te_hat_all_blin$mean, te_hat_all_bsc$mean))

  # Remove outliers
  low.blin <- quantile(tc_hat_all_blin$mean, .01)
  hi.blin <- quantile(tc_hat_all_blin$mean, .99)
  low.bsc <- quantile(tc_hat_all_bsc$mean, .01)
  hi.bsc <- quantile(tc_hat_all_bsc$mean, .99)

  range_blin <- between(tcs[, 1], low.blin, hi.blin)
  range_bsc <- between(tcs[, 2], low.bsc, hi.bsc)

  tcs <- tcs %>% dplyr::filter(range_blin & range_bsc)
  tes <- tes %>% dplyr::filter(range_blin & range_bsc)

  # Calculate cost-effectiveness results ---------------------------------------
  ec_table <- calculate_ce_summary(df_tc = tcs, df_te = tes)
  ec_table$discount <- discount_str
  ec_table$scenario <- ifelse(is.null(scenario), "base_case", scenario)

  # Generate PSA results -------------------------------------------------------
  strategies <- c("Blinatumomab and chemotherapy", "Chemotherapy only")
  psa_blina <- dampack::make_psa_obj(cost = tcs, effectiveness = tes, strategies = strategies)
  psa_blina$strategies <- strategies
  colnames(psa_blina$effectiveness) <- colnames(psa_blina$cost) <- strategies
  psa_df <- summary(psa_blina)
  psa_df$discount <- discount_str
  psa_df$scenario <- ifelse(is.null(scenario), "base_case", scenario)

  # Plot and save CE plot
  plot_ce <- plot(psa_blina, alpha = 1)
  plot_ce[["theme"]][["legend.position"]] <- "bottom"

  # ICER table ----------------------------------------------------------------
  blina_cea <- calculate_icers(
    cost = colMeans(tcs),
    effect = colMeans(tes),
    strategies = c("Blinatumomab", "Standard of Care")
  )
  cea_tbl <- as.data.frame(darthtools::format_table_cea(blina_cea))
  cea_tbl$discount <- discount_str
  cea_tbl$scenario <- ifelse(is.null(scenario), "base_case", scenario)

  # CEAC -----------------------------------------------------------------------
  ceac_blina <- ceac(wtp = seq(0, 200000, 5000), psa = psa_blina)
  ceac_blina$discount <- discount_str
  ceac_blina$scenario <- ifelse(is.null(scenario), "base_case", scenario)

  # Plot and save CEAC
  plot_ceac <- plot(ceac_blina, frontier = F)
  plot_ceac[["theme"]][["legend.position"]] <- "bottom"
  plot_ceac[["labels"]][["y"]] <- "Probability of Cost-Effective"
  plot_ceac[["data"]][["Strategy"]] <- strategies

  # Save plots only for base case and ignore plots for scenarios ---------------
  if (is.null(scenario)) {
    # Save CE plot
    ggplot2::ggsave(
      filename = paste0(root, "/ce-plot-", discount_str, ".pdf"),
      plot = plot_ce,
      width = 10, height = 7
    )

    # Save CEAC plot
    ggplot2::ggsave(
      filename = paste0(root, "/ceac-plot-", discount_str, ".pdf"),
      plot = plot_ceac,
      width = 10, height = 7
    )
  }

  # Return all results
  return(list(
    ec_table = ec_table,
    psa_df = psa_df,
    cea_tbl = cea_tbl,
    ceac_df = ceac_blina,
    scenario = scenario
  ))
}
