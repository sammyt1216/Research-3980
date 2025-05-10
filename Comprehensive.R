library(dplyr)
library(matrixStats)
library(fixest)
library(modelsummary)
library(ggplot2)

a_line_DiD <- read.csv("A_Line_DiD.csv")
g_line_DiD <- read.csv("G_Line_DiD.csv")
SMART_DiD <- read.csv("SMART_DiD.csv")
sunrail_DiD <- read.csv("SunRail_DiD.csv")

a_line_DiD <- a_line_DiD %>%
  mutate(line = "A Line")
g_line_DiD <- g_line_DiD %>%
  mutate(line = "G Line")
SMART_DiD <- SMART_DiD %>%
  mutate(line = "SMART")
SunRail <- sunrail_DiD %>%
  mutate(line = "SunRail")

# Combine your datasets
all_data <- bind_rows(a_line_DiD, g_line_DiD,SMART_DiD,sunrail_DiD)

# Add a priority column and arrange by GISJOIN, year, and priority
all_data <- all_data %>%
  mutate(priority = ifelse(line == "A Line", 1, 2)) %>%
  arrange(GISJOIN, YEAR, priority)

# Remove duplicates, keeping the first occurrence
final_data <- all_data %>%
  group_by(GISJOIN, YEAR) %>%
  slice(1) %>%
  ungroup() %>%
  select(-priority)

# Regressions for all_data
alldata_DiD_TT <- feols(log(med_travel_time) ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_dens <- feols(density ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_car <- feols(pct_car ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_transit <- feols(pct_transit ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_mhi <- feols(log(mhi) ~ TREATMENT : ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_medvalhu <- feols(log(medvalhu) ~ TREATMENT : ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)
alldata_DiD_unem <- feols(pct_unem ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN)

# Print results
alldata_results <- etable(list(alldata_DiD_dens,alldata_DiD_medvalhu,alldata_DiD_mhi,alldata_DiD_unem,
                               alldata_DiD_car,alldata_DiD_transit,alldata_DiD_TT), tex = TRUE)
cat(alldata_results)
