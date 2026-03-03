#' Calculate Cost-Effectiveness Summary Table
#'
#' This function computes summary statistics for costs and effects of two
#' interventions (e.g., chemotherapy only vs. blinatumomab and chemotherapy),
#' including means and 95% confidence intervals for each group and their
#' differences. It also calculates the Incremental Cost-Effectiveness Ratio
#' (ICER).
#'
#' @param df_tc Data frame containing cost values for each simulation/sample.
#' Columns should correspond to intervention and control groups.
#' @param df_te Data frame containing effect values for each simulation/sample.
#' Columns should correspond to intervention and control groups.
#' @param col_intervention Integer specifying the column index for the
#' intervention group (default is 1).
#' @param col_control Integer specifying the column index for the control group
#' (default is 2).
#'
#' @return A data frame summarizing mean and 95% confidence intervals for costs,
#' effects, and their differences, as well as the ICER.
#'
calculate_ce_summary <- function(df_tc, df_te, col_intervention = 1, col_control = 2) {
  
  # Differences between intervention and control for costs and effects ---------
  df_tc$diff <- df_tc[, col_intervention] - df_tc[, col_control]
  df_te$diff <- df_te[, col_intervention] - df_te[, col_control]

  # Calculate costs ------------------------------------------------------------
  c_int_mean <- mean(df_tc[, col_intervention], na.rm = TRUE)
  c_int_ci <- quantile(df_tc[, col_intervention], c(0.025, 0.975), na.rm = TRUE)
  c_int_summary <- sprintf("%.2f (%.2f, %.2f)", c_int_mean, c_int_ci[1], c_int_ci[2])

  c_ctrl_mean <- mean(df_tc[, col_control], na.rm = TRUE)
  c_ctrl_ci <- quantile(df_tc[, col_control], c(0.025, 0.975), na.rm = TRUE)
  c_ctrl_summary <- sprintf("%.2f (%.2f, %.2f)", c_ctrl_mean, c_ctrl_ci[1], c_ctrl_ci[2])

  c_diff_mean <- mean(df_tc$diff, na.rm = TRUE)
  c_diff_ci <- quantile(df_tc$diff, c(0.025, 0.975), na.rm = TRUE)
  c_diff_summary <- sprintf("%.2f (%.2f, %.2f)", c_diff_mean, c_diff_ci[1], c_diff_ci[2])

  # Calculate effects ----------------------------------------------------------
  e_int_mean <- mean(df_te[, col_intervention], na.rm = TRUE)
  e_int_ci <- quantile(df_te[, col_intervention], c(0.025, 0.975), na.rm = TRUE)
  e_int_summary <- sprintf("%.2f (%.2f, %.2f)", e_int_mean, e_int_ci[1], e_int_ci[2])

  e_ctrl_mean <- mean(df_te[, col_control], na.rm = TRUE)
  e_ctrl_ci <- quantile(df_te[, col_control], c(0.025, 0.975), na.rm = TRUE)
  e_ctrl_summary <- sprintf("%.2f (%.2f, %.2f)", e_ctrl_mean, e_ctrl_ci[1], e_ctrl_ci[2])

  e_diff_mean <- mean(df_te$diff, na.rm = TRUE)
  e_diff_ci <- quantile(df_te$diff, c(0.025, 0.975), na.rm = TRUE)
  e_diff_summary <- sprintf("%.2f (%.2f, %.2f)", e_diff_mean, e_diff_ci[1], e_diff_ci[2])

  # Calculate ICER -------------------------------------------------------------
  icer <- df_tc$diff / df_te$diff
  # icer_mean <- mean(icer, na.rm = TRUE)
  icer_mean <- c_diff_mean / e_diff_mean # Standard way to calculate ICER
  icer_ci <- quantile(icer, c(0.025, 0.975), na.rm = TRUE)

  # Create summary table
  results <- rbind(
    cbind("Costs", c_ctrl_summary, c_int_summary, c_diff_summary),
    cbind("QALYs", e_ctrl_summary, e_int_summary, e_diff_summary),
    cbind("ICER", "", "", sprintf("%.2f", icer_mean))
  ) |>
    as.data.frame()

  colnames(results) <- c("Metric", "Chemotherapy only", "Blinatumomab and chemotherapy", "Incremental")

  return(results)
}
