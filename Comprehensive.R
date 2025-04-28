library(dplyr)
library(matrixStats)
library(fixest)
library(modelsummary)
library(ggplot2)

a_line_DiD <- read.csv("A_Line_DiD.csv")
g_line_DiD <- read.csv("G_Line_DiD.csv")
SMART_DiD <- read.csv("SMART_DiD.csv")

a_line_DiD <- a_line_DiD %>%
  mutate(line = "A Line")
g_line_DiD <- g_line_DiD %>%
  mutate(line = "G Line")
SMART_DiD <- SMART_DiD %>%
  mutate(line = "SMART")

# Combine your datasets
all_data <- bind_rows(a_line_DiD, g_line_DiD,SMART_DiD)

## I'd like some help with this filtering please, I'm worried treatment data got taken out as the graph I made is inconsistent with my DiD analysis

priority <- all_data %>%
  group_by(GISJOIN) %>%
  arrange(line) %>%  # A Line alphabetically before G Line
  slice(1) %>%
  ungroup() %>%
  select(GISJOIN, line)  # This shows which line to keep

# Join back onto the full dataset
all_data_clean <- all_data %>%
  left_join(priority, by = "GISJOIN", suffix = c("", ".chosen")) %>%
  filter(line == line.chosen) %>%
  select(-line.chosen)

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
denver_data <- bind_rows(a_line_DiD, g_line_DiD)
denpriority <- denver_data %>%
  group_by(GISJOIN) %>%
  arrange(line) %>%  # A Line alphabetically before G Line
  slice(1) %>%
  ungroup() %>%
  select(GISJOIN, line)  # This shows which line to keep

denver_data_clean <- denver_data %>%
  left_join(priority, by = "GISJOIN", suffix = c("", ".chosen")) %>%
  filter(line == line.chosen) %>%
  select(-line.chosen)
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

# Add fitted values to your datasets
denver_data_clean <- denver_data_clean %>%
  mutate(dens_predicted = predict(denver_ols_dens), medvalhu_predicted = predict(denver_ols_medvalhu),
         mhi_predicted = predict(denver_ols_mhi), pct_unem_predicted = predict(denver_ols_unem))

denver_treat <- denver_treat %>%
  mutate(dens_predicted = predict(denver_ols_dens_treat), medvalhu_predicted = predict(denver_ols_medvalhu_treat),
         mhi_predicted = predict(denver_ols_mhi_treat), pct_unem_predicted = predict(denver_ols_unem_treat))

denver_control <- denver_control %>%
  mutate(dens_predicted = predict(denver_ols_dens_control), medvalhu_predicted = predict(denver_ols_medvalhu_control),
         mhi_predicted = predict(denver_ols_mhi_control), pct_unem_predicted = predict(denver_ols_unem_control))

denver_plot_data <- bind_rows(
  denver_treat %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                          mhi, mhi_predicted, pct_unem, pct_unem_predicted) %>% 
    mutate(Group = "Treatment"),
  denver_control %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                            mhi, mhi_predicted, pct_unem, pct_unem_predicted) %>% 
    mutate(Group = "Control"),
  denver_data_clean %>% select(YEAR, density, dens_predicted, medvalhu, medvalhu_predicted,
                               mhi, mhi_predicted, pct_unem, pct_unem_predicted) %>% 
    mutate(Group = "Total")
)

# Plot with ggplot
ggplot(denver_plot_data, aes(x = YEAR, y = density, color = Group)) +
  geom_line(aes(y = ifelse(Group == "Treatment", dens_predicted, 
                           ifelse(Group == "Control", dens_predicted, dens_predicted))),
            size = 1.2) + # Regression lines for each group
  geom_vline(xintercept = 2016, linetype = "dashed", color = "black") + # A Line operational
  geom_vline(xintercept = 2019, linetype = "dashed", color = "black") + # G Line operational
  labs(title = "Pooled OLS: Density over Time by Group",
       x = "Year",
       y = "Population Density",
       color = "Group") +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  scale_color_manual(values = c("Treatment" = "blue", "Control" = "red", "Total" = "green"))
