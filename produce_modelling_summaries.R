# This is a script to easily source a summary for Rain and Flow data.
# The rain function requires the data to contain columns in the following order:

# sars_per_person, Date, avg_rain, avg_rain.lag2, avg_rain.lag3, avg_rain.lag4, avg_rain.lag5, avg_rain.lag6, avg_rain.lag7

# , And Flow requires the data to be in the following form:

library(tidyverse)
library(lubridate)
library(mgcv)
library(scales)



# Rain

produce_rainsummary.gam = function(data, .date_knots_overide = 12){
  # This is a function which takes a dataset from the ones produced above, fits a model, produces summary statistics and a map.

  # Fit the model
  model = gam(data = data %>% mutate(Date = as.numeric(as.Date(Date))),
              # 12 knots to account for monthly variation due to outbreaks.
              sars_per_person ~ s(Date, k = .date_knots_overide) + avg_rain + avg_rain.lag1 +
                avg_rain.lag2 +
                avg_rain.lag3 +
                avg_rain.lag4 +
                avg_rain.lag5 +
                avg_rain.lag6 +
                avg_rain.lag7,
              family = 'quasipoisson')


  # Statistics

  # Variance
  var = vcov(model)[2:9,2:9] %>% sum()
  # Estimate (beta)
  beta = coef(model)[2:9] %>% sum()
  # T-Statistic
  t = beta / sqrt(var)
  # P-Value (2-sided)
  model.p = 2 * min(pnorm(beta / sqrt(var)), 1 - pnorm(beta / sqrt(var)))
  # Estimate (response)
  model.est = 1 - beta %>% exp()
  # Confidence Interval
  model.confint = 1 - (beta + c(-1, 1) * 2 * sqrt(var)) %>% exp()



  # Graphics

  # Get ratio for first vis conversion ()

  # Rain Time-Series Plot
  rain_to_sarspp = with(data, max(sars_per_person, na.rm = T) / max(avg_rain, na.rm = T))
  location = data$SampleLocation %>% unique()

  p1 = data %>%
    mutate(avg_rain_trans = avg_rain * rain_to_sarspp) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = sars_per_person)) +
    geom_segment(aes(xend = Date, yend = avg_rain_trans, y = 0), col = 'blue', alpha = 0.5) +
    scale_y_continuous(sec.axis = sec_axis(~ . / rain_to_sarspp, name = paste('Average Rain over', location, '(mm)'))) +
    labs(y='Sars-Cov2 Genome Copies per Litre-Person-Population', x = 'Date',
         title = 'SARS-COV2 Copies with Average Rainfall by Date',
         subtitle = paste(location, '(# genomes per person-population)')) +
    theme_minimal()


  # Hex plot

  p2 = data %>%
    ggplot(aes(x = avg_rain, y = sars_per_person)) +
    geom_hex() +
    geom_smooth() +
    theme_minimal() +
    labs(title = 'SARS-COV2 Copies per Litre-Person-Population against Average Rainfall',
         x = 'Average Rain', y = 'Sars-Cov2 Genome Copies per Litre-Person-Population',
         subtitle = location, fill = 'Count')


  # Dilution plot

  p3 = tibble(
    location = rep(location, 2),
    lower = rep(model.confint[1], 2),
    upper = rep(model.confint[2], 2),
    est_type = c('est', 'est_c'),
    est = c(model.est, 1-model.est)) %>%
    ggplot(aes(x = location, y = est, fill = est_type)) +
    geom_col(position = 'stack') +
    geom_errorbar(aes(ymin = 1-lower, ymax = 1-upper), alpha = 0.2) +
    theme_bw() +
    labs(title = 'Observed Concentration of SARS-COV2 Genome copies',
         subtitle = 'Per (mm) of Rainwater for each of the past 7 days', x = location,
         y = 'Proportion of Actual', fill = 'Concentration') +
    scale_fill_manual(values = hue_pal()(2), labels = c('Conc. lost', 'Observed Conc.')) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme(axis.text.x = element_blank())


  # Map

  ### Figure OUT MAP


  # Smooth Plot

  smooth_plots = plot.gam_mine(model)
  p4 = smooth_plots[[2]] +
    labs(title = 'Date Smoother',
         subtitle = 'To account for variation explained by time',
         x = 'Date', y = 'Smooth Coefficient')


  # Output

  list(
    model = model,
    site = location,
    var = var,
    beta_effect = model.est,
    t = t,
    p = model.p,
    confint_95 = model.confint,
    timeseries_plot = p1,
    rain_sars_hexplot = p2,
    concentration_plot = p3,
    smooth_plot = p4
  )
}


print_rainsummary.gam = function(model){
  print(tibble(
    `Concentration Loss %` = round(model$beta_effect * 100, 3),
    `97.5%` = round(model$confint_95[1]*100, 3),
    `2.5%` = round(model$confint_95[2]*100, 3),
    `P-Value` = model$p
  ) %>% knitr::kable(caption = model$site))

  print(model$timeseries_plot)
  print(model$smooth_plot)
  print(model$rain_sars_hexplot)
  print(model$concentration_plot)
}









# Flow



produce_flowsummary.gam = function(data, .date_knots_overide = 12){
  # This is a function which takes a dataset from the ones produced above, fits a flow model, produces summary statistics and a map.

  # Fit the model
  model = gam(data = data %>% mutate(Date = as.numeric(as.Date(Date))),
              # 12 knots to account for monthly variation due to outbreaks.
              sars_per_person ~ s(Date, k = .date_knots_overide) +
                daily_flow_rate_m3,
              family = 'quasipoisson')


  # Statistics

  flow.summary = summary(model)$p.table[rownames(summary(model)$p.table) == 'daily_flow_rate_m3',]

  # Variance
  var = flow.summary[2]^2
  # Estimate (beta)
  beta = flow.summary[1]
  # Beta confint
  beta.confint = flow.summary[1] + c(-1, 1) * 2 * sqrt(var)
  # P-Value (2-sided)
  model.p = flow.summary[4]
  # Estimate (response)
  model.est = 1 - beta %>% exp()
  # Confidence Interval
  model.confint = 1 - beta.confint %>% exp()



  # Graphics

  # Get ratio for first vis conversion ()

  # Flow Time-Series Plot
  flow_to_sarspp = with(data, max(sars_per_person, na.rm = T) / max(daily_flow_rate_m3, na.rm = T))
  location = data$SampleLocation %>% unique()

  p1 = data %>%
    mutate(daily_flow_rate_m3_trans = daily_flow_rate_m3 * flow_to_sarspp) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = sars_per_person)) +
    geom_segment(aes(xend = Date, yend = daily_flow_rate_m3_trans, y = 0), col = 'red', alpha = 0.3) +
    scale_y_continuous(sec.axis = sec_axis(~ . / flow_to_sarspp, name = paste('Daily Flow Rate for', location, '(m3)'))) +
    labs(y='Sars-Cov2 Genome Copies per Litre-Person-Population', x = 'Date',
         title = 'SARS-COV2 Copies with Daily Flow by Date',
         subtitle = paste(location, '(# genomes per person-population)')) +
    theme_minimal()


  # Hex plot

  p2 = data %>%
    ggplot(aes(x = daily_flow_rate_m3, y = sars_per_person)) +
    geom_hex() +
    geom_smooth(col = 'red') +
    theme_minimal() +
    labs(title = 'SARS-COV2 Copies per Litre-Person-Population against Daily Flow',
         x = 'Daily Flow (m3)', y = 'Sars-Cov2 Genome Copies per Litre-Person-Population',
         subtitle = location, fill = 'Count') +
    scale_fill_gradient(low = 'black', high = 'red')



  # Dilution plot

  average_flow = median(data$daily_flow_rate_m3, na.rm = T)
  average_flow_effect = beta * average_flow
  average_flow_effect.confint = beta.confint * average_flow

  p3 = tibble(
    location = rep(location, 2),
    lower = rep(exp(average_flow_effect.confint[1]), 2),
    upper = rep(exp(average_flow_effect.confint[2]), 2),
    est_type = c('est', 'est_c'),
    est = c(1 - exp(average_flow_effect), exp(average_flow_effect))) %>%
    ggplot(aes(x = location, y = est, fill = est_type)) +
    geom_col(position = 'stack') +
    geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    theme_bw() +
    labs(title = 'Observed Concentration of SARS-COV2 Genome copies',
         subtitle = 'For the average (median) amount of Daily Flow (m3)', x = location,
         y = 'Proportion of Actual', fill = 'Concentration') +
    scale_fill_manual(values = hue_pal()(2), labels = c('Conc. lost', 'Observed Conc.')) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    theme(axis.text.x = element_blank())


  # Map

  ### Figure OUT MAP


  # Smooth Plot

  smooth_plots = plot.gam_mine(model)
  p4 = smooth_plots[[2]] +
    labs(title = 'Date Smoother',
         subtitle = 'To account for variation explained by time',
         x = 'Date', y = 'Smooth Coefficient')

  # Output

  list(
    model = model,
    site = location,
    flow.lq = quantile(data$daily_flow_rate_m3, 0.25, na.rm = T),
    flow.uq = quantile(data$daily_flow_rate_m3, 0.75, na.rm = T),
    flow.median = median(data$daily_flow_rate_m3, na.rm = T),
    var = var,
    beta = beta,
    p = model.p,
    beta_confint_95 = beta.confint,
    timeseries_plot = p1,
    rain_flow_hexplot = p2,
    concentration_plot = p3,
    smooth_plot = p4
  )
}


print_flowsummary.gam = function(model){
  print(tibble(
    Quantile = c('Lower 25%', 'Median', 'Upper 75%'),
    `Flow Rate` = c(model$flow.lq, model$flow.median, model$flow.uq),
    `% Concentration Loss Explained` = round((1 - exp(model$beta * `Flow Rate`)) * 100, 6),
    `97.5%` = round((1 - exp(model$beta_confint_95[1] * `Flow Rate`)) * 100, 6),
    `2.5%` = round((1 - exp(model$beta_confint_95[2] * `Flow Rate`)) * 100, 6),
    `P-Value` = model$p
  ) %>% knitr::kable(caption = model$site))

  print(model$timeseries_plot)
  print(model$smooth_plot)
  print(model$rain_flow_hexplot)
  print(model$concentration_plot)
}
























# Flow/Rain
# Produces Summary of Using Flow over Rain.
# P-value usually the only thing interpretable.

produce_flowrainsummary.gam = function(data, .date_knots_overide = 12){
  # This is a function which takes a dataset from the ones produced above, fits a flow model, produces summary statistics and a map.

  # Fit the model
  model = gam(data = data %>% mutate(Date = as.numeric(as.Date(Date))),
              # 12 knots to account for monthly variation due to outbreaks.
              sars_per_person ~ s(Date, k = .date_knots_overide) +
                avg_rain +
                avg_rain.lag1 +
                avg_rain.lag2 +
                avg_rain.lag3 +
                avg_rain.lag4 +
                avg_rain.lag5 +
                avg_rain.lag6 +
                avg_rain.lag7 +
                daily_flow_rate_m3,
              family = 'quasipoisson')


  # Statistics

  flow.summary = summary(model)$p.table[rownames(summary(model)$p.table) == 'daily_flow_rate_m3',]

  # Variance
  var = flow.summary[2]^2
  # Estimate (beta)
  beta = flow.summary[1]
  # Beta confint
  beta.confint = flow.summary[1] + c(-1, 1) * 2 * sqrt(var)
  # P-Value (2-sided)
  model.p = flow.summary[4]
  # Estimate (response)
  model.est = 1 - beta %>% exp()
  # Confidence Interval
  model.confint = 1 - beta.confint %>% exp()



  # Graphics

  # Get ratio for first vis conversion ()

  # Rain/Flow Time-Series Plot
  rain_to_sarspp = with(data, max(sars_per_person, na.rm = T) / max(avg_rain, na.rm = T))
  flow_to_sarspp = with(data, max(sars_per_person, na.rm = T) / max(daily_flow_rate_m3, na.rm = T))
  location = data$SampleLocation %>% unique()

  p1 = data %>%
    mutate(avg_rain_trans = avg_rain * rain_to_sarspp,
           daily_flow_rate_m3_trans = daily_flow_rate_m3 * flow_to_sarspp) %>%
    ggplot(aes(x = Date)) +
    geom_line(aes(y = sars_per_person)) +
    geom_segment(aes(xend = Date, yend = avg_rain_trans, y = 0), col = 'blue', alpha = 0.5) +
    geom_segment(aes(xend = Date, yend = daily_flow_rate_m3_trans, y = 0), col = 'red', alpha = 0.3) +
    scale_y_continuous(sec.axis = sec_axis(~ . / flow_to_sarspp, name = paste('Daily Flow Rate for', location, '(m3)'))) +
    labs(y='Sars-Cov2 Genome Copies per Litre-Person-Population', x = 'Date',
         title = 'SARS-COV2 Copies with Average Rainfall and Daily Flow by Date',
         subtitle = paste(location, '(# genomes per person-population)')) +
    theme_minimal()


  # Hex rain-flow plot

  p2 = data %>%
    ggplot(aes(x = avg_rain, y = daily_flow_rate_m3)) +
    geom_hex() +
    theme_minimal() +
    labs(title = 'SARS-COV2 Copies per Litre-Person-Population against Average Rainfall',
         x = 'Average Rain', y = 'Daily Flow Rate (m3)',
         subtitle = location, fill = 'Count') +
    scale_fill_gradient(low = 'black', high = 'purple')



  # Output

  list(
    model = model,
    site = location,
    flow.lq = quantile(data$daily_flow_rate_m3, 0.25, na.rm = T),
    flow.uq = quantile(data$daily_flow_rate_m3, 0.75, na.rm = T),
    flow.median = median(data$daily_flow_rate_m3, na.rm = T),
    var = var,
    beta = beta,
    p = model.p,
    beta_confint_95 = beta.confint,
    timeseries_plot = p1,
    rain_flow_hexplot = p2
  )
}


print_flowrainsummary.gam = function(model){
  print(tibble(
    Quantile = c('Lower 25%', 'Median', 'Upper 75%'),
    `Flow Rate` = c(model$flow.lq, model$flow.median, model$flow.uq),
    `% Concentration Loss Explained` = round((1 - exp(model$beta * `Flow Rate`)) * 100, 6),
    `97.5%` = round((1 - exp(model$beta_confint_95[1] * `Flow Rate`)) * 100, 6),
    `2.5%` = round((1 - exp(model$beta_confint_95[2] * `Flow Rate`)) * 100, 6),
    `P-Value` = model$p
  ) %>% knitr::kable(caption = model$site))

  print(model$timeseries_plot)
  print(model$rain_flow_hexplot)
}








# GAM plot

# GAM Smooth function which does the plot(gam) function but converts the numerics on the bottom to dates.

plot.gam_mine = function(model, .exp = F){
  pd = plot.gam(model, select = 0)
  plot.data = tibble(
    y = pd[[1]]$fit %>% as.numeric(),
    y_lower = y - 2 * pd[[1]]$se,
    y_upper = y + 2 * pd[[1]]$se,
    x = as.Date(pd[[1]]$x, origin = '1970-01-01'))

  # Get breaks for vis
  breaks = seq(min(plot.data$x), max(plot.data$x), by = 'months')
  labels = paste(breaks %>% month(label = T), breaks %>% year())


  p1 = plot.data %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(col = 'red') +
    geom_line(aes(y = y_lower), col = 'black', lty = 2) +
    geom_line(aes(y = y_upper), col = 'black', lty = 2) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    theme_bw()

  p1_exp = plot.data %>%
    ggplot(aes(x = x, y = exp(y))) +
    geom_line(col = 'red') +
    geom_line(aes(y = exp(y_lower)), col = 'black', lty = 2) +
    geom_line(aes(y = exp(y_upper)), col = 'black', lty = 2) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    # Account for large exponentiated interval.
    coord_cartesian(ylim = c(0, max(plot.data$y %>% exp())*2)) +
    theme_bw() +
    labs(caption = '*Exponentiated Y axis.')

  list(p1, p1_exp)
}


