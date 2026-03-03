#' Calculate survival probabilities across multiple simulations
#'
#' Reads simulation output files, extracts first event times, fits a KM
#' survival curve for each simulation, and returns summary statistics
#' (median, mean, 95% CI) at specified time points.
#'
#' @param events Character vector of event names to consider (default: "DEAD").
#' @param n_sim Integer; number of simulations to process.
#' @param output.path String; path prefix for simulation data files.
#' @param times Numeric vector of time points at which to evaluate survival.
#'
#' @return A data frame with columns:
#'   \item{time}{Time points.}
#'   \item{surv_median}{Median survival probability across simulations.}
#'   \item{surv_mean}{Mean survival probability.}
#'   \item{surv_lcl}{Lower 95% confidence limit.}
#'   \item{surv_ucl}{Upper 95% confidence limit.}
#'
#' @import survival furrr
#' @export
#'
calc_sim_surv_probs <- function(events = c("DEAD"), n_sim, output.path, times) {

  # Load necessary libraries
  require(survival)
  require(furrr)

  # Setup parallel processing to read in data
  plan(multisession)

  # Make sure that the times vector is the same length as the .rds file columns
  # since the rds files have time as the columns and individuals as rows
  if (length(times) != ncol(readRDS(paste0(output.path, "data.1.rds")))) {
    stop("Length of 'times' does not match the number of columns in the data files.")
  }

  # Print to console what is happening
  cat("Calculating survival probabilities for events:", paste(events, collapse = ", "), "\n")

  # Calculate survival probabilities for each simulation
  # - Iterate over each simulation
  # - Build data frame with columns being each simulation survival probabilities
  df_surv <- future_map_dfc(
    .x = 1:n_sim,
    .f = ~ {
      dat <- paste0(output.path, "data.", .x, ".rds") |> readRDS()

      # Extract the first occurrence of each event for each individual
      # - rows are individuals, columns are time points
      event_time <- apply(dat, 1, function(x) {
        idx <- which(x %in% events)
        if (length(idx) > 0) times[idx[1]] else NA # First occurrence of event
      })

      # Clean up memory
      dat <- NULL
      gc()

      # Create a data frame with event_time, event, and time
      df <- data.frame(
        event_time = event_time,
        event = ifelse(is.na(event_time), 0, 1), # 0 if no event
        time = ifelse(is.na(event_time), max(times), event_time)
      )

      # Survival fit of the data
      surv_fit <- survfit(Surv(time, event) ~ 1, data = df)

      # Extract the survival probabilities at each time point
      surv_probs <- summary(surv_fit, times = times)$surv

      # Check if the length of surv_probs matches the length of times, and if 
      # not, issue a warning
      if (length(surv_probs) != length(times)) {
        warning(
          paste(
            "Simulation", .x, "has inconsistent survival probabilities length:", 
            length(surv_probs), "vs expected", length(times), ". Filling in with last known value."
          )
        )
      }

      # Ensure consistent length - if surv_probs is shorter than times, extend it
      if (length(surv_probs) < length(times)) {
        # Pad with last known value to match length of times
        surv_probs <- c(surv_probs, rep(tail(surv_probs, 1), length(times) - length(surv_probs)))
      }

      # Create data frame with simulation number
      surv_df <- data.frame(x = surv_probs)
      colnames(surv_df) <- paste0("sim_", .x)

      return(surv_df)
    }, .progress = TRUE
  )

  # Calculate median, mean, and confidence intervals across simulations
  surv_med <- apply(df_surv, 1, median, na.rm = TRUE)
  surv_mean <- apply(df_surv, 1, mean, na.rm = TRUE)
  surv_lcl <- apply(df_surv, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
  surv_ucl <- apply(df_surv, 1, function(x) quantile(x, 0.975, na.rm = TRUE))

  # Put results into a final data frame
  df_surv_final <- data.frame(
    time = times,
    surv_median = surv_med,
    surv_mean = surv_mean,
    surv_lcl = surv_lcl,
    surv_ucl = surv_ucl
  )

  # Print completion message
  cat("\nSurvival probabilities calculated for", n_sim, "simulations.\n")

  # Return the final data frame
  return(df_surv_final)

}
