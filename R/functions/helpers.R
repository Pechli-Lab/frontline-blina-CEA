#' Check if a variable or value is TRUE (safe for missing names)
#'
#' This function checks whether:
#' - an unquoted variable name refers to an existing variable that is exactly TRUE, or
#' - a quoted name refers to an existing variable that is exactly TRUE, or
#' - a provided value itself is exactly TRUE.
#'
#' @param x Unquoted variable name, quoted variable name (character scalar), or a value.
#' @param envir Environment to look up variable names (default: parent.frame()).
#' @param inherits Whether to search parent environments (default: TRUE).
#'
#' @return TRUE if condition holds, otherwise FALSE.
#'
#' @examples
#' is_TRUE(TRUE) # TRUE
#' is_TRUE(FALSE) # FALSE
#' is_TRUE(xx) # FALSE if xx doesn't exist; TRUE only if xx exists and is TRUE
#' is_TRUE("xx") # same behavior using a quoted name
#' flag <- TRUE
#' is_TRUE(flag) # TRUE
#' flag <- FALSE
#' is_TRUE(flag) # FALSE
is_TRUE <- function(x, envir = parent.frame(), inherits = TRUE) {
  expr <- substitute(x)

  # Case 1: quoted variable name for `x`
  if (is.character(expr) && length(expr) == 1) {
    nm <- expr
    if (exists(nm, envir = envir, inherits = inherits)) {
      return(isTRUE(get(nm, envir = envir, inherits = inherits)))
    }
    return(FALSE)
  }

  # Case 2: unquoted variable name (symbol), do NOT evaluate it
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (exists(nm, envir = envir, inherits = inherits)) {
      return(isTRUE(get(nm, envir = envir, inherits = inherits)))
    }
    return(FALSE)
  }

  # Case 3: direct value or expression (evaluate safely)
  val <- tryCatch(eval(expr, envir = envir), error = function(...) FALSE)

  return(isTRUE(val))
}

#' Trim Outliers from Data by Setting Extreme Values to NA
#'
#' @description
#'   Removes outliers by setting values outside specified quantile thresholds
#'   to NA. Works with both vectors and data frames. When applied to a data
#'   frame, trims outliers independently for each column.
#'
#' @param x Numeric vector or data frame. Data to trim for outliers.
#' @param prob Numeric vector of length 2. Quantile probabilities defining the
#'   lower and upper bounds for outlier detection. Default is c(0.01, 0.99),
#'   which removes the bottom 1% and top 1% of values.
#'
#' @return
#'   Returns the same type as input (vector or data frame) with outlier values
#'   replaced by NA. For data frames, each column is trimmed independently.
#'
#' @details
#'   The function identifies outliers as values falling below the lower quantile
#'   threshold (prob[1]) or above the upper quantile threshold (prob[2]).
#'   Missing values (NA) are ignored when calculating quantiles.
#'
#'   When applied to a data frame, the function recursively applies itself to
#'   each column while preserving the data frame structure.
#'
#' @examples
#' # Trim outliers from a vector (default 1% and 99% quantiles)
#' x <- c(1:100, 1000)
#' trim_outliers(x)
#'
#' # Custom quantile thresholds (remove bottom 5% and top 5%)
#' trim_outliers(x, prob = c(0.05, 0.95))
#'
#' # Trim outliers from a data frame
#' df <- data.frame(a = c(1:100, 1000), b = c(-50, 1:99, 500))
#' trim_outliers(df)
#'
trim_outliers <- function(x, prob = c(0.01, 0.99)) {
  # x: data to trim (can be vector or data frame)
  
  if (is.data.frame(x)) {
    # If x is a data frame, apply trim_outliers to each column
    x[] <- lapply(x, trim_outliers, prob = prob) # keeps the data frame structure then applies the function to each column
    return(x)
  }

  # For vectors: calculate the lower and upper bounds
  bounds <- quantile(x, prob, na.rm = TRUE)

  # Replace values outside bounds with NA
  x[x < bounds[1] | x > bounds[2]] <- NA

  # Return the trimmed vector
  x
}

#' Summarize Total Cost and Effects Output from Simulation Results
#'
#' @description
#'   Calculates summary statistics (mean, standard error, and confidence
#'   intervals) for total cost and effects outputs from simulation results. This
#'   function is designed to work with vectors or lists.
#'
#' @param x Numeric vector or list. Total values from simulation iterations.
#'   If a list is provided, it will be unlisted before calculations.
#' @param prob Numeric vector of length 2. Probabilities for the lower and upper
#'  bounds of the confidence interval. Default is c(0.025, 0.975) for a 95% CI.
#'
#' @return
#'   A data frame with one row and four columns:
#'   \describe{
#'     \item{tc_hat}{Mean total cost across all simulations}
#'     \item{se}{Standard error of the total cost}
#'     \item{upper}{Upper bound of the confidence interval}
#'     \item{lower}{Lower bound of the confidence interval}
#'   }
#'
#' @details
#'   The function computes confidence intervals using specific quantiles
#'   percentiles of the cost distribution. Row names are removed from the output
#'   data frame for cleaner presentation.
#'
summarize_sim_output <- function(x, prob = c(0.025, 0.975)) {
  x_unlisted <- unlist(x)
  
  # Calculate summary statistics
  mean_tc <- mean(x_unlisted, na.rm = TRUE) # total cost mean
  se_tc <- sd(x_unlisted, na.rm = TRUE)     # total cost standard error
  upper_tc <- quantile(x_unlisted, prob[2], na.rm = TRUE) # upper CI
  lower_tc <- quantile(x_unlisted, prob[1], na.rm = TRUE) # lower CI
  
  # Create output data frame
  df <- data.frame(
    tc_hat = mean_tc,
    se = se_tc,
    upper = upper_tc,
    lower = lower_tc
  )
  rownames(df) <- NULL
  return(df)
}

#' Calculate Event Risk from Simulation Data
#'
#' @description
#'   Calculates the proportion of patients who experienced a specific event
#'   (e.g., relapse, death) within a specified time horizon from microsimulation
#'   results.
#'
#' @param data Matrix of health states. Rows represent individual patients,
#'   columns represent time points. Each cell contains the health state at that
#'   time point. Columns should correspond to time in days.
#' @param event Character string. The health state of interest (e.g., "REL1"
#'   for first relapse, "DEAD" for death).
#' @param event_time Numeric or NULL. Time horizon in years for evaluating
#'   event risk (e.g., 5, 10, 20 years). If NULL, considers the entire follow-up
#'   period. Default is NULL.
#'
#' @return Numeric. The proportion of patients who experienced the event,
#'   ranging from 0 to 1.
#'
#' @details
#'   The function operates on microsimulation output where time is tracked in
#'   days (30 days per month, 12 months per year = 360 days per year).
#'
#'   When `event_time` is specified, the data matrix is subset to include only
#'   time points up to that horizon before calculating event incidence.
#'
#'   The function requires a global variable `times` representing the time points
#'   (in days) corresponding to the columns of the data matrix.
#' 
event_risk <- function(data, event, event_time = NULL) {

  # data: matrix of states (rows: patients, columns: time points)
  # event: state of interest (e.g., "REL1")
  # event_time: time horizon in years (e.g., 10 years)
  
  # Trim the data to only be interested in times of interest
  if (!is.null(event_time)) {

    # Subset to columns up to event_time (in years but matrix is in days)
    max_day <- event_time * 30 * 12

    data <- data[, times <= max_day]
  }

  # Identify patients who experienced the event
  id_event <- apply(data, 1, function(x) {
    any(x %in% event)
  })

  # Calculate the proportion of patients who experienced the event
  n_events <- sum(id_event)
  prop_events <- n_events / nrow(data)

  return(prop_events)
}


#' Simulate Risk Differences Between Two Proportions
#'
#' Calculates the risk difference between two proportions using a two-sample
#' proportion test, estimates the standard error from the 95% confidence
#' interval, and simulates uncertainty in the risk difference using a normal
#' approximation.
#'
#' @param numerator Number of events in each group (intervention, control).
#' @param denominator Denominator (sample size) for each group.
#' @param n_sim Number of simulations for uncertainty.
#'
#' @return Numeric vector of simulated risk differences uncertainty.
#'
#' @examples
#' # Simulate risk differences for two groups
#' risk_diff_prop(c(10, 20), c(100, 100), n_sim = 1000)
#'
sim_risk_diff <- function(numerator, denominator = denom_AALL1731, n_sim = n_sim) {
  prop_test <- prop.test(x = numerator, n = denominator)
  rd <- unname(prop_test$estimate[1] - prop_test$estimate[2]) # risk difference
  se <- abs((rd - prop_test$conf.int[1]) / qnorm(0.975)) # standard error
  sim_rd <- rnorm(n_sim, rd, se) # rd uncertainty using normal approximation
  return(sim_rd) # return vector of simulated risk differences
}

#' Simulate Weighted Average Duration for Adverse Events
#'
#' This function simulates the weighted average duration of adverse events.
#' It takes in trial data that are CEAE grade 3+ (_g3p) and grade 4 (_g4),
#' then simulates durations using beta distributions. It can handle inputs
#' that differ in hospitalization length of stay for each of the grades, then
#' calculates a weighted average duration.
#'
#' @param n_g3p Numeric. The number of grade 3+ events. If it is a vector,
#'   it will be summed. Contains both grade 3 and grade 4 events.
#' @param dur_g3 Numeric. Hospitalization duration for grade 3 events.
#' @param n_g4 Numeric. The number of grade 4 events. If it is a vector,
#'   it will be summed.
#' @param dur_g4 Numeric. Hospitalization duration for grade 4 events.
#' @param n_sim Integer. The number of simulations to perform. Default is 500.
#'
#' @return A numeric vector of length \code{n_sim} containing the simulated
#'   weighted average durations.
#'
#' @details
#' - The function first checks if \code{n_g3p} and \code{n_g4} are vectors
#'   and sums them if necessary.
#' - It then simulates durations for grade 3+ and grade 4 events using beta
#'   distributions.
#' - Finally, it calculates the weighted average duration by combining the
#'   simulated durations.
#'
#' @examples
#' # Example usage
#' result <- sim_weighted_ae_day(n_g3p = 10, dur_g3 = 5, n_g4 = 8,
#'   dur_g4 = 3, n_sim = 100)
#' summary(result)
#' 
sim_weighted_ae_day <- function(n_g3p, dur_g3, n_g4, dur_g4, n_sim = 500, seed = SEED) {
  set.seed(seed)

  # Check if inputs are vectors > 1, if so, sum them to get totals
  if (length(n_g3p) > 1) {
    n_g3p <- sum(n_g3p)
  }
  if (length(n_g4) > 1) {
    n_g4 <- sum(n_g4)
  }

  # Simulate durations using beta distributions
  dur_g3_sim <- rbeta(n_sim, n_g3p - n_g4, n_g4          ) * dur_g3 # duration of grade 3 events
  dur_g4_sim <- rbeta(n_sim, n_g4        , (n_g3p - n_g4)) * dur_g4 # duration of grade 4 events

  # Calculate weighted average duration
  dur_weighted <- dur_g3_sim + dur_g4_sim

  return(dur_weighted)
}


#' Calculate Life Expectancy Metrics for Two Treatment Arms
#'
#' This function reads life expectancy (or similar) data from RDS files for two
#' treatment arms (BSC and Blina), calculates summary statistics, and returns a
#' formatted data frame for reporting.
#'
#' @param rds Character string specifying the RDS filename to read.
#'   Default is "les.rds".
#' @param name Character string specifying the object name within the RDS file
#'   to extract. Default is "le_i".
#' @param out_var Character string specifying the column name in the output data
#'   frame. Default is "le".
#' @param description Character string describing the metric. Default is
#'   "Life Expectancy (years)".
#'
#' @return A data frame with three columns:
#'   \item{arm}{Treatment arm (BSC or Blina)}
#'   \item{<out_var>}{Formatted string with median and 95% CI:
#'     "median (lcl, ucl)"}
#'   \item{metric}{Description of the metric}
#'
#' @details
#' The function assumes that global variables `output_path_bsc` and
#' `output_path_blin` are defined, which specify the paths to the directories
#' containing the RDS files for each treatment arm.
#'
#' @importFrom dplyr mutate select
#' @importFrom stats quantile median mean
#'
calc_le_metrics <- function(bsc_path, blin_path, rds = "les.rds", name = "le_i", out_var = "le", description = "Life Expectancy (years)") {
  
  # BSC
  le_i_v <- readRDS(paste0(bsc_path, "/", rds))[[name]]
  le_i_bsc <- data.frame(
    mean = mean(le_i_v),
    median = median(le_i_v),
    lcl = quantile(le_i_v, 0.025),
    ucl = quantile(le_i_v, 0.975)
  ) |>
    mutate(arm = "BSC")

  # Blina
  le_i_v <- readRDS(paste0(blin_path, "/", rds))[[name]]
  le_i_blin <- data.frame(
    mean = mean(le_i_v),
    median = median(le_i_v),
    lcl = quantile(le_i_v, 0.025),
    ucl = quantile(le_i_v, 0.975)
  ) |>
    mutate(arm = "Blina")

  # Summary table
  le_summary <- rbind(le_i_bsc, le_i_blin) |>
    dplyr::mutate(
      est_chr = paste0(
        round(median, 2),
        " (", round(lcl, 2), ", ", round(ucl, 2), ")"
      ),
      metric = description
    )

  rownames(le_summary) <- NULL

  df_out <- le_summary %>%
    dplyr::select(arm, !!out_var := est_chr, metric)

  return(df_out)
}


#' Format Confidence Intervals with Custom Formatting
#'
#' @description
#'   Formats confidence interval strings by extracting numeric values and
#'   applying specified formatting (currency, number, or percentage). Handles
#'   vectorized input and returns formatted strings in the pattern:
#'   "estimate (lower CI separator upper CI)".
#'
#' @param x Character vector. Strings containing confidence intervals in the
#'   format "estimate (lower, upper)" where estimate, lower, and upper are
#'   numeric values.
#' @param accuracy Numeric. Number of decimal places or rounding increment for
#'   formatting. Default is 1.
#' @param scale Numeric. Scaling factor for percentages. Use 1 if input is
#'   already in percentage form (0-100), or 100 if converting from proportion
#'   (0-1). Default is 1.
#' @param ci_sep Character string. Separator between lower and upper confidence
#'   bounds. Default is " to ".
#' @param type Character string. Type of formatting to apply. Options are:
#'   \describe{
#'     \item{"dollar"}{Currency formatting with $ symbol}
#'     \item{"number"}{Plain number formatting with commas}
#'     \item{"pct"}{Percentage formatting with % symbol}
#'   }
#'   Default is "dollar".
#'
#' @return Character vector of formatted confidence intervals.
#'
#' @details
#'   The function extracts three numeric values from each input string using
#'   regex pattern matching. It expects the format "number (number, number)".
#'   The extracted values are formatted according to the specified type and
#'   recombined into a readable string.
#'
#' @examples
#' # Format as currency (default)
#' format_ci("272967.34 (242316.97, 311862.27)")
#' #> "$272,967 ($242,317 to $311,862)"
#'
#' # Format as numbers with 1 decimal place
#' format_ci("3.07 (1.57, 5.98)", accuracy = 0.1, type = "number")
#' #> "3.1 (1.6 to 6.0)"
#'
#' # Format as percentages
#' format_ci("0.45 (0.32, 0.58)", accuracy = 0.1, scale = 100, type = "pct")
#' #> "45.0% (32.0% to 58.0%)"
#'
#' # Vectorized usage in a data frame
#' library(dplyr)
#' df <- data.frame(ci = c("100 (90, 110)", "200 (180, 220)"))
#' df |> mutate(formatted = format_ci(ci, type = "number"))
#'
#' @importFrom stringr str_extract_all
#' @importFrom scales dollar number percent
#'
format_ci <- function(
  x,
  accuracy = 1,
  scale = 1,
  ci_sep = " to ",
  type = "dollar"
) {
  # Wrap in sapply to handle vectorized input
  sapply(
    x,
    function(val) {
      # Extract numbers
      nums <- as.numeric(str_extract_all(val, "\\d+\\.?\\d*")[[1]])

      if (type == "dollar") {
        # Format numbers as currency
        formatted <- scales::dollar(nums, accuracy = accuracy)
      } else if (type == "number") {
        # Format numbers as regular numbers
        formatted <- scales::number(nums, accuracy = accuracy)
      } else if (type == "pct") {
        # Format numbers as percentages
        formatted <- scales::percent(nums, accuracy = accuracy, scale = scale)
      } else {
        stop("Unsupported type. Use 'dollar', 'number', or 'pct'.")
      }

      # Combine into desired format
      out <- sprintf(
        paste0("%s (%s", ci_sep, "%s)"),
        formatted[1],
        formatted[2],
        formatted[3]
      )
    },
    USE.NAMES = FALSE
  )
}
