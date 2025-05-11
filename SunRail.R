library(dplyr)
library(matrixStats)
library(fixest)
library(modelsummary)
library(ggplot2)

sunrail <- read.csv("SunRail.csv")
orlando_control <- read.csv("OrlandoControl.csv")

weighted_median_interpolated <- function(counts, times) {
  total <- sum(counts)
  if (total == 0) return(NA)
  
  cum_counts <- cumsum(counts)
  below <- which(cum_counts < total / 2)
  above <- which(cum_counts >= total / 2)[1]
  
  if (length(below) == 0) {
    return(times[1])  # median is in the first bin
  }
  
  # Bin where the median lies
  i <- above
  
  # Interpolation between bin i-1 and i
  count_below <- if (i > 1) cum_counts[i - 1] else 0
  bin_count <- counts[i]
  prop_in_bin <- (total / 2 - count_below) / bin_count
  
  # Width of the bin (time[i] is bin midpoint, so we approximate bin width)
  if (i == 1) {
    lower <- times[i]
    upper <- (times[i] + times[i + 1]) / 2
  } else if (i == length(times)) {
    lower <- (times[i - 1] + times[i]) / 2
    upper <- times[i]
  } else {
    lower <- (times[i - 1] + times[i]) / 2
    upper <- (times[i] + times[i + 1]) / 2
  }
  bin_width <- upper - lower
  
  return(lower + prop_in_bin * bin_width)
}

# Define the weights vector
TTweights <- c(2.5, 7, 12, 17, 22, 27, 32, 37, 42, 52, 74.5, 95)  # Row weights

# Subset the rows (QTHE002 to QTHE013)
TTsubset_cols <- c("QTHE002", "QTHE003", "QTHE004", "QTHE005", "QTHE006",
                   "QTHE007", "QTHE008", "QTHE009", "QTHE010", "QTHE011", 
                   "QTHE012", "QTHE013")

# Filter the dataframe for the 13 rows and relevant columns
sunrailTTdata <- as.matrix(sunrail[,TTsubset_cols])

# Compute weighted medians for each column
sunrail_TT_medians <- apply(sunrailTTdata, 1, function(counts) {
  weighted_median_interpolated(as.numeric(counts), TTweights)
})

# Input the results into the table
sunrail$med_travel_time <- sunrail_TT_medians

ocontrolTTdata <- as.matrix(orlando_control[,TTsubset_cols])
ocontrol_TT_medians <- apply(ocontrolTTdata, 1, function(counts) {
  weighted_median_interpolated(as.numeric(counts), TTweights)
})
orlando_control$med_travel_time <- ocontrol_TT_medians

sunrail$TREATMENT <- 1
orlando_control$ACTIVE <- 0
orlando_control$TREATMENT <- 0

orlando_control <- orlando_control %>%
  mutate(ACTIVE = case_when(
    YEAR == 2014 ~ 0.6939891,
    YEAR > 2014 ~ 1,
    TRUE ~ 0
  ))
sunrail_DiD <- bind_rows(sunrail,orlando_control)
# Travel time DiD
sunrail_DiD_TT <- feols(log(med_travel_time) ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# Density DiD
sunrail_DiD$density <- sunrail_DiD$pop / sunrail_DiD$ALAND

sunrail_DiD_dens <- feols(log(density) ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# avg vehicles DiD
sunrail_DiD <- sunrail_DiD %>% 
  mutate(avg_vehicles = (QZBE004 + QZBE011 + (QZBE005 + QZBE012)*2 + (QZBE006 + QZBE013)*3 + (QZBE007 + QZBE014)*4 + (QZBE008 + QZBE015)*5)/QZBE001)
sunrail_DiD_vehicles <- feols(avg_vehicles ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# drive alone DiD
sunrail_DiD <- sunrail_DiD %>%
  mutate(pct_car = QTFE003/QTHE001)
sunrail_DiD_car <- feols(pct_car ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# transit DiD
sunrail_DiD <- sunrail_DiD %>%
  mutate(pct_transit = QTFE010/QTHE001)
sunrail_DiD_transit <- feols(pct_transit ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# mhi DiD
sunrail_DiD <- sunrail_DiD %>%
  mutate(mhi = as.numeric(mhi)) %>%
  filter(mhi >= 0)
sunrail_DiD_mhi <- feols(log(mhi) ~ TREATMENT : ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# medvalhu DiD
sunrail_DiD <- sunrail_DiD %>%
  mutate(medvalhu = as.numeric(medvalhu)) %>%
  filter(medvalhu >= 0)
sunrail_DiD_medvalhu <- feols(log(medvalhu) ~ TREATMENT : ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# unemployment DiD
sunrail_DiD <- sunrail_DiD %>%
  mutate(pct_unem = QXSE005/QXSE003)
sunrail_DiD_unem <- feols(pct_unem ~ TREATMENT:ACTIVE | GISJOIN + YEAR, data = sunrail_DiD, cluster = ~GISJOIN)

# collect and print results
sunrail_results <- etable(list(sunrail_DiD_dens,sunrail_DiD_medvalhu,sunrail_DiD_mhi,sunrail_DiD_unem,
                          sunrail_DiD_car,sunrail_DiD_transit,sunrail_DiD_TT), tex = TRUE)
cat(sunrail_results)
write.csv(sunrail_DiD,"SunRail_DiD.csv")