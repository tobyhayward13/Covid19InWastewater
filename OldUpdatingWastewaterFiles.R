# This is the code I used to run before I was confident pulling branches from the ESR data.
# This was because I was phobic of github before and knew no other way to get updated waste data, other than webscraping it.
# It is much, much better this way. Just pull the data and get BUSY.
# Don't use this code unless you desperately need to.

library(tidyverse)
library(rvest)

# Waste
update_wastewater_files = function(){
  # Get file names
  waste.files = str_split("cases_national.csv cases_regional.csv cases_site.csv sites.csv ww_data_all.csv ww_national.csv ww_regional.csv ww_site.csv", ' ') %>% unlist()
  # list.files('data.test/')

  # Site to read from
  site = 'https://raw.githubusercontent.com/ESR-NZ/covid_in_wastewater/main/data/'

  # Paste files to site name to read from internet
  sites = paste0(site, waste.files)

  # Convert html to text and store in memory
  sites.text = map(sites, read_html) %>% map_chr(html_text)

  # Write to file
  for (i in 1:length(waste.files)) write(sites.text[i], paste0('data/', str_sub(waste.files, 1, -5)[i], '.csv'))
}

update_wastewater_files()

