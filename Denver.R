library(dplyr)
library(matrixStats)
library(fixest)

a_line <- read.csv("A_Line.csv")
g_line <- read.csv("G_Line.csv")
control <- read.csv("DenverControl.csv")

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
a_lineTTdata <- as.matrix(a_line[,TTsubset_cols])
g_lineTTdata <- as.matrix(g_line[,TTsubset_cols])

# Compute weighted medians for each column
a_line_TT_medians <- apply(a_lineTTdata, 1, function(counts) {
  weighted_median_interpolated(as.numeric(counts), TTweights)
})
g_line_TT_medians <- apply(g_lineTTdata, 1, function(counts) {
  weighted_median_interpolated(as.numeric(counts), TTweights)
})

# Input the results into the table
a_line$med_travel_time <- a_line_TT_medians
g_line$med_travel_time <- g_line_TT_medians

controlTTdata <- as.matrix(control[,TTsubset_cols])
control_TT_medians <- apply(controlTTdata, 1, function(counts) {
  weighted_median_interpolated(as.numeric(counts), TTweights)
})
control$med_travel_time <- control_TT_medians

a_line$TREATMENT <- 1
g_line$TREATMENT <- 1
control$ACTIVE <- 0
control$TREATMENT <- 0

control_a_line <- control[control$YEAR <= 2021,] %>%
  mutate(ACTIVE = case_when(
    YEAR == 2016 ~ 0.6939891,
    YEAR > 2016 ~ 1,
    TRUE ~ 0
  ))
control_g_line <- control[control$YEAR >= 2016,] %>%
  mutate(ACTIVE = case_when(
    YEAR == 2019 ~ 0.6849315,
    YEAR > 2019 ~ 1,
    TRUE ~ 0
  ))

a_line_DiD <- bind_rows(a_line, control_a_line)
g_line_DiD <- bind_rows(g_line, control_g_line)

# Travel time DiD
a_line_DiD_TT <- feols(log(med_travel_time) ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_TT <- feols(log(med_travel_time) ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# Density DiD
a_line_DiD$density <- a_line_DiD$pop / a_line_DiD$ALAND
g_line_DiD$density <- g_line_DiD$pop / g_line_DiD$ALAND

a_line_DiD_dens <- feols(density ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_dens <- feols(density ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# avg vehicles DiD
a_line_DiD <- a_line_DiD %>% 
  mutate(avg_vehicles = (QZBE004 + QZBE011 + (QZBE005 + QZBE012)*2 + (QZBE006 + QZBE013)*3 + (QZBE007 + QZBE014)*4 + (QZBE008 + QZBE015)*5)/QZBE001)
g_line_DiD <- g_line_DiD %>% 
  mutate(avg_vehicles = (QZBE004 + QZBE011 + (QZBE005 + QZBE012)*2 + (QZBE006 + QZBE013)*3 + (QZBE007 + QZBE014)*4 + (QZBE008 + QZBE015)*5)/QZBE001)
a_line_DiD_vehicles <- feols(avg_vehicles ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_vehicles <- feols(avg_vehicles ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# drive alone DiD
a_line_DiD <- a_line_DiD %>%
  mutate(pct_car = QTFE003/QTHE001)
g_line_DiD <- g_line_DiD %>%
  mutate(pct_car = QTFE003/QTHE001)
a_line_DiD_car <- feols(pct_car ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_car <- feols(pct_car ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# transit DiD
a_line_DiD <- a_line_DiD %>%
  mutate(pct_transit = QTFE010/QTHE001)
g_line_DiD <- g_line_DiD %>%
  mutate(pct_transit = QTFE010/QTHE001)
a_line_DiD_transit <- feols(pct_transit ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_transit <- feols(pct_transit ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# mhi DiD
a_line_DiD$mhi <- as.numeric(a_line_DiD$mhi)
g_line_DiD$mhi <- as.numeric(g_line_DiD$mhi)
a_line_DiD_mhi <- feols(log(mhi) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_mhi <- feols(log(mhi) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# medvalhu DiD
a_line_DiD$medvalhu <- as.numeric(a_line_DiD$medvalhu)
g_line_DiD$medvalhu <- as.numeric(g_line_DiD$medvalhu)
a_line_DiD_medvalhu <- feols(log(medvalhu) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_medvalhu <- feols(log(medvalhu) ~ TREATMENT * ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# unemployment DiD
a_line_DiD <- a_line_DiD %>%
  mutate(pct_unem = QXSE005/QXSE003)
g_line_DiD <- g_line_DiD %>%
  mutate(pct_unem = QXSE005/QXSE003)
a_line_DiD_unem <- feols(pct_unem ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = a_line_DiD)
g_line_DiD_unem <- feols(pct_unem ~ TREATMENT*ACTIVE | GISJOIN + YEAR, data = g_line_DiD)

# collect and print results
a_line_results <- etable(a_line_DiD_dens,a_line_DiD_medvalhu,a_line_DiD_mhi,a_line_DiD_unem,
                         a_line_DiD_vehicles,a_line_DiD_car,a_line_DiD_transit,a_line_DiD_TT)
write.csv(a_line_results,"A Line Regression Results.csv")