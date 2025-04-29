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
alldata_DiD_TT <- feols(log(med_travel_time) ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_dens <- feols(density ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_vehicles <- feols(avg_vehicles ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_car <- feols(pct_car ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_transit <- feols(pct_transit ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_mhi <- feols(log(mhi) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_medvalhu <- feols(log(medvalhu) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)
alldata_DiD_unem <- feols(pct_unem ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = all_data_clean, cluster = ~GISJOIN + YEAR)

# Print results
alldata_results <- etable(list(alldata_DiD_dens,alldata_DiD_medvalhu,alldata_DiD_mhi,alldata_DiD_unem,
                               alldata_DiD_vehicles,alldata_DiD_car,alldata_DiD_transit,alldata_DiD_TT), tex = TRUE)
cat(alldata_results)

# Graphs
weighted_mean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

denver_data <- bind_rows(a_line_DiD, g_line_DiD)
# Add a priority column and arrange by GISJOIN, year, and priority
denver_data <- denver_data %>%
  mutate(priority = ifelse(line == "A Line", 1, 2)) %>%
  arrange(GISJOIN, YEAR, priority)

# Remove duplicates, keeping the first occurrence
denver_data_clean <- denver_data %>%
  group_by(GISJOIN, YEAR) %>%
  slice(1) %>%
  ungroup() %>%
  select(-priority)

denver_treat <- denver_data_clean %>%
  filter(TREATMENT == 1)
denver_control <- denver_data_clean %>%
  filter(TREATMENT == 0)

denver_ols_dens <- feols(density ~ YEAR, data = denver_data_clean)
denver_ols_dens_treat <- feols(density ~ YEAR, data = denver_treat)
denver_ols_dens_control <- feols(density ~ YEAR, data = denver_control)

denver_ols_medvalhu <- feols(medvalhu ~ YEAR, data = denver_data_clean)
denver_ols_medvalhu_treat <- feols(medvalhu ~ YEAR, data = denver_treat)
denver_ols_medvalhu_control <- feols(medvalhu ~ YEAR, data = denver_control)

denver_ols_mhi <- feols(mhi ~ YEAR, data = denver_data_clean)
denver_ols_mhi_treat <- feols(mhi ~ YEAR, data = denver_treat)
denver_ols_mhi_control <- feols(mhi ~ YEAR, data = denver_control)

denver_ols_unem <- feols(pct_unem ~ YEAR, data = denver_data_clean)
denver_ols_unem_treat <- feols(pct_unem ~ YEAR, data = denver_treat)
denver_ols_unem_control <- feols(pct_unem ~ YEAR, data = denver_control)

denver_ols_tt <- feols(med_travel_time ~ YEAR, data = denver_data_clean)
denver_ols_tt_treat <- feols(med_travel_time ~ YEAR, data = denver_treat)
denver_ols_tt_control <- feols(med_travel_time ~ YEAR, data = denver_control)

# Group by YEAR and treatment status (or control)
denver_weighted_averages <- denver_data_clean %>%
  mutate(group = case_when(
    TREATMENT == 1 ~ "Treatment",
    TREATMENT == 0 ~ "Control"
  )) %>%
  group_by(YEAR, group) %>%
  summarise(
    weighted_density = weighted_mean(density, pop), # assuming 'pop' is your weight
    weighted_medvalhu = weighted_mean(medvalhu,QX7E001), # QX7E001 is total housing units (occupied + vacant)
    weighted_mhi = weighted_mean(mhi,pop),
    weighted_unem = weighted_mean(pct_unem,QXSE003), # QXSE003 is civilian LF
    weighted_tt = weighted_mean(med_travel_time,QTHE001), # QTHE001 is workers who did not work from home
    .groups = 'drop'
  )

# Add total group if you want
denver_weighted_total <- denver_data_clean %>%
  group_by(YEAR) %>%
  summarise(
    weighted_density = weighted_mean(density, pop),
    group = "Total",
    .groups = 'drop'
  )

# Combine all
denver_weighted_averages_all <- bind_rows(denver_weighted_averages, denver_weighted_total)

# Add fitted values to your datasets
denver_data_clean <- denver_data_clean %>%
  mutate(dens_predicted = predict(denver_ols_dens), medvalhu_predicted = predict(denver_ols_medvalhu),
         mhi_predicted = predict(denver_ols_mhi), pct_unem_predicted = predict(denver_ols_unem),
         tt_predicted = predict(denver_ols_tt))

denver_treat <- denver_treat %>%
  mutate(dens_predicted = predict(denver_ols_dens_treat), medvalhu_predicted = predict(denver_ols_medvalhu_treat),
         mhi_predicted = predict(denver_ols_mhi_treat), pct_unem_predicted = predict(denver_ols_unem_treat),
         tt_predicted = predict(denver_ols_tt_treat))

denver_control <- denver_control %>%
  mutate(dens_predicted = predict(denver_ols_dens_control), medvalhu_predicted = predict(denver_ols_medvalhu_control),
         mhi_predicted = predict(denver_ols_mhi_control), pct_unem_predicted = predict(denver_ols_unem_control),
         tt_predicted = predict(denver_ols_tt_control))

denver_plot_data <- bind_rows(
  denver_treat %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                          mhi, mhi_predicted, pct_unem, pct_unem_predicted, med_travel_time, tt_predicted) %>% 
    mutate(group = "Treatment"),
  denver_control %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                            mhi, mhi_predicted, pct_unem, pct_unem_predicted, med_travel_time, tt_predicted) %>% 
    mutate(group = "Control"),
  denver_data_clean %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                               mhi, mhi_predicted, pct_unem, pct_unem_predicted, med_travel_time, tt_predicted) %>% 
    mutate(group = "Total")
)

# Plot with ggplot
ggplot() +
  # Solid lines for weighted averages
  geom_line(data = weighted_averages_all, aes(x = YEAR, y = weighted_density, color = group), size = 1) +
  
  # Dashed lines for pooled OLS trends
  geom_smooth(data = denver_plot_data, aes(x = YEAR, y = dens_predicted, color = group), 
              method = "lm", se = FALSE, linetype = "dashed", size = 0.8) +
  
  # Add vertical lines for A Line and G Line openings
  geom_vline(xintercept = 2016, linetype = "dotted", color = "black", size = 0.8) +
  geom_vline(xintercept = 2019, linetype = "dotted", color = "black", size = 0.8) +
  
  labs(
    x = "Year",
    y = "Density (people per sq meter)",
    color = "Group"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +  # whole years
  
  # Clean up theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
