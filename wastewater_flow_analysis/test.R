models$AU_Eastern
print_flowrainsummary.gam(models$AU_Eastern)

x = map(models, ~c(.x$site, .x$flow.median, 1-exp(.x$flow.median * .x$beta), .x$p)) %>% do.call(rbind, .) %>% as.data.frame()
colnames(x) = c('Sample Location', 'Median Flow (m3)', 'Estimated Loss in Concentration', 'p-value')

x = x %>% arrange(`p-value`)

write_csv(x, 'final report/flowrain_results.csv')
