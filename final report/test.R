models = model.data.flow[-4] %>% map(produce_flowsummary.gam)

data = map(models, ~c(.x$site, .x$flow.median, 1 - ((.x$flow.median * .x$beta * -1 ) %>% exp()), .x$p)) %>% do.call(rbind, .) %>% as_tibble()

colnames(data) = c('Sample Location', 'Median Flow (m3)', 'Estimated Loss in Concentration', 'p-value')

data = mutate(data, across(.cols = 2:4, ~as.numeric(.x))) %>% arrange(`p-value`)

write_csv(data, 'final report/flow_results.csv')
