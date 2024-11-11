library(tidyverse)
library(pscl)
library(splines)
library(here)
# source(here("src/Utilities.R"))
# source(here("src/bayes/model.R"))
# source(here("src/bayes/bayes_fun.R"))


prepare_covid_data <- function(data) {
  data %>%
    group_by(FIPS) %>%
    mutate(
      log_pop = log(population),
      date_num = as.numeric(date),
      tmmx_scaled = scale(tmmx),
      rmax_scaled = scale(rmax)
    ) %>%
    ungroup()
}

create_model_formula <- function(outcome, lag_names) {
  as.formula(paste(
    outcome,
    "~ ns(tmmx_scaled, df = 2) +",
    "ns(rmax_scaled, df = 2) +",
    "ns(date_num, df = 4) +",
    "dayofweek +",
    "relative_change_feb +",
    paste(lag_names, collapse = " + "),
    "| 1"
  ))
}

fit_zinb_model <- function(
    formula, 
    data, 
    lag_matrix) {
  lag_matrix <<- lag_matrix
  with(list(data = data, lag_matrix = lag_matrix), {
    zeroinfl(
      formula,
      dist = "negbin",
      data = data,
      offset = log_pop,
      control = zeroinfl.control(EM = FALSE)
    )
  })
}

calculate_cumulative_effect <- function(model, lag_names, pm_delta = 10) {
  coefs <- coef(model)[lag_names]
  vcov_subset <- vcov(model)[lag_names, lag_names]
  
  cumulative <- sum(coefs)
  std_error <- sqrt(sum(vcov_subset))
  
  list(
    percent_change = (exp(cumulative * pm_delta) - 1) * 100,
    ci_lower = (exp((cumulative - 1.96 * std_error) * pm_delta) - 1) * 100,
    ci_upper = (exp((cumulative + 1.96 * std_error) * pm_delta) - 1) * 100
  )
}

extract_lag_effects <- function(model, lag_names) {
  coefs <- coef(model)[lag_names]
  ses <- sqrt(diag(vcov(model)[lag_names, lag_names]))
  
  data.frame(
    lag = seq_along(lag_names),
    effect = coefs,
    se = ses
  )
}

run_covid_analysis <- 
  function(data, outcome = "deaths", max_lag = 14) {
  clean_data <- prepare_covid_data(data)
  
  lag_matrix <- as.matrix(create.lag.value(
    dff = clean_data,
    value = "pm25",
    group = "FIPS",
    lags = 1:max_lag
  ))
  
  formula <- 
    create_model_formula(outcome, 
                         lag_names = "lag_matrix")
  
  model <- tryCatch({
    fit_zinb_model(formula, clean_data, lag_matrix)
  }, error = function(e) {
    return(list(error = e$message))
  })
  
  if ("error" %in% names(model)) {
    return(model)
  }

  lag_names <- paste0("count_lag_matrix.l", 1:max_lag) 
  cumulative <- 
    calculate_cumulative_effect(model, lag_names)
  lag_effects <- 
    extract_lag_effects(model, lag_names)
  
  list(
    model = model,
    cumulative_effect = cumulative,
    lag_effects = lag_effects,
    formula = formula
  )
}

plot_lag_pattern <- function(lag_effects) {
  ggplot(lag_effects, aes(x = lag, y = effect)) +
    geom_point() +
    geom_errorbar(
      aes(ymin = effect - 1.96*se, ymax = effect + 1.96*se),
      width = 0.2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(
      x = "Lag (days)",
      y = "Effect estimate",
      title = "Lag-specific PM2.5 effects"
    )
}
