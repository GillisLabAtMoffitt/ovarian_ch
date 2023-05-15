# Import Library
library(tidyverse)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "Ovarian CH")

Demographic <- 
  readxl::read_xlsx(paste0(path, "/raw data/Report for ovarian patients .xlsx"),
                    sheet = "Demographics") %>% 
  janitor::clean_names()

ovarian_DNA <- 
  readxl::read_xlsx(paste0(path, "/raw data/Report for ovarian patients .xlsx"),
                    sheet = "Biobanking") %>% 
  janitor::clean_names()

ovarian_info <- 
  readxl::read_xlsx(paste0(path, "/raw data/Report for ovarian patients .xlsx"),
                    sheet = "Cancer Characteristics") %>% 
  janitor::clean_names()

Chemot <- 
  readxl::read_xlsx(paste0(path, "/raw data/Report for ovarian patients .xlsx"),
                    sheet = "Treatments") %>% 
  janitor::clean_names()


