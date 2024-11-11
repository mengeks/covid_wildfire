library(here)
source(here("src/Utilities.R"))
source(here("src/xmeng-model-utils.R"))
dff <- load.data()
results <- 
  run_covid_analysis(
    data = dff, 
    outcome = "deaths", 
    max_lag = 14)
cumulative_effect <- results$cumulative_effect
print(paste0("Percent change (95% CI): ",
            round(cumulative_effect$percent_change, 1), "% (",
            round(cumulative_effect$ci_lower, 1), "%, ",
            round(cumulative_effect$ci_upper, 1), "%)"))
# The result: "Percent change (95% CI): 3.7% (1.4%, 6%)": give an interpretation 
plot_lag_pattern(lag_effects = results$lag_effects)
ggsave(here("output", "xmeng","plot_lag_pattern.png"))
