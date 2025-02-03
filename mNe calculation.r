# Install and load necessary libraries
if (!require(moments)) install.packages("moments", dependencies = TRUE)
library(moments)
# Load the data
data_filename <- "~/NGS_analyses/fastsimcoal/GemHumMon2401_10_boot/GemHumMon_10_Recent_boot"
data <- read.table(data_filename, header = TRUE)
calculate_confidence_intervals <- function(column_data, method = "auto") {
  # Skip processing if the data size is too small
  if (length(column_data) < 3) {
    return(list(lower = NA, upper = NA, method = "insufficient_data", mean = NA, median = NA))
  }
  
  # Remove NA values
  column_data <- na.omit(column_data)
  
  # Calculate skewness and kurtosis of the distribution
  skewness_value <- skewness(column_data)
  kurtosis_value <- kurtosis(column_data)
  
  # Automatically select method if not specified
  if (method == "auto") {
    if (!is.na(skewness_value) && !is.na(kurtosis_value) && abs(skewness_value) < 1 && abs(kurtosis_value - 3) < 1) {
      method <- "normal"
    } else {
      method <- "bootstrap"
    }
  }
  
  # Compute confidence intervals
  if (method == "normal") {
    # Confidence interval based on normal distribution
    mean_value <- mean(column_data)
    std_error <- sd(column_data) / sqrt(length(column_data))
    conf_interval <- mean_value + qnorm(c(0.025, 0.975)) * std_error
    method_used <- "normal"
  } else if (method == "bootstrap") {
    # Confidence interval using bootstrap method
    n_iterations <- 1000
    means <- numeric(n_iterations)
    for (i in 1:n_iterations) {
      sample <- sample(column_data, size = length(column_data), replace = TRUE)
      means[i] <- mean(sample)
    }
    conf_interval <- quantile(means, probs = c(0.025, 0.975))
    mean_value <- mean(column_data)
    method_used <- "bootstrap"
  } else if (method == "quantile") {
    # Confidence interval using quantile method
    conf_interval <- quantile(column_data, probs = c(0.025, 0.975))
    mean_value <- mean(column_data)
    method_used <- "quantile"
  } else {
    stop("Invalid method specified.")
  }
  
  median_value <- median(column_data)
  
  return(list(lower = conf_interval[1], upper = conf_interval[2], method = method_used, mean = mean_value, median = median_value))
}
# Creating new calculated columns
data$Calculation1_mNe_01_all <- data$NPOP0 * data$RSANC1_gem * data$MIG_01
data$Calculation2_mNe_02_all <- data$NPOP0 * data$RSANC1_gem * data$MIG_02
data$Calculation3_mNe_20_all <- data$NPOP2 * data$RSANC1_mon * data$MIG_20
data$Calculation4_mNe_21_all <- data$NPOP2 * data$RSANC1_mon * data$MIG_21
data$Calculation5_mNe_10_all <- data$NPOP1 * data$RSANC1_hum * data$MIG_10
data$Calculation6_mNe_12_all <- data$NPOP1 * data$RSANC1_hum * data$MIG_12
## 5 generations
data$Calculation1_mNe_01_5gene <- ifelse(data$CHANGM < 5, 
                            data$NPOP0 * data$RSANC1_gem * data$MIG_01, 
                            NA)

data$Calculation2_mNe_02_5gene <- ifelse(data$CHANGM <= 5, 
                            data$NPOP0 * data$RSANC1_gem * data$MIG_02, 
                            NA)
data$Calculation3_mNe_20_5gene <- ifelse(data$CHANGM <= 5, 
                            data$NPOP2 * data$RSANC1_mon * data$MIG_20, 
                            NA)
data$Calculation4_mNe_21_5gene <- ifelse(data$CHANGM <= 5, 
                            data$NPOP2 * data$RSANC1_mon * data$MIG_21, 
                            NA)
data$Calculation5_mNe_10_5gene <- ifelse(data$CHANGM <= 5, 
                            data$NPOP1 * data$RSANC1_hum * data$MIG_10, 
                            NA)
data$Calculation6_mNe_12_5gene <- ifelse(data$CHANGM <= 5, 
                            data$NPOP1 * data$RSANC1_hum * data$MIG_12, 
                            NA)
## 0 generations
data$Calculation1_mNe_01_0gene <- ifelse(data$CHANGM <= 0, 
                              data$NPOP0 * data$RSANC1_gem * data$MIG_01, 
                              NA)

data$Calculation2_mNe_02_0gene <- ifelse(data$CHANGM <= 0 , 
                              data$NPOP0 * data$RSANC1_gem * data$MIG_02, 
                              NA)
data$Calculation3_mNe_20_0gene <- ifelse(data$CHANGM <= 0, 
                              data$NPOP2 * data$RSANC1_mon * data$MIG_20, 
                              NA)
data$Calculation4_mNe_21_0gene <- ifelse(data$CHANGM <= 0, 
                              data$NPOP2 * data$RSANC1_mon * data$MIG_21, 
                              NA)
data$Calculation5_mNe_10_0gene <- ifelse(data$CHANGM <= 0, 
                              data$NPOP1 * data$RSANC1_hum * data$MIG_10, 
                              NA)
data$Calculation6_mNe_12_0gene <- ifelse(data$CHANGM <= 0, 
                              data$NPOP1 * data$RSANC1_hum * data$MIG_12, 
                              NA)
# Calculate confidence intervals for each column
columns_of_interest <- c("Calculation1_mNe_01_all", "Calculation2_mNe_02_all", "Calculation3_mNe_20_all", "Calculation4_mNe_21_all", "Calculation5_mNe_10_all", "Calculation6_mNe_12_all",
                         "Calculation1_mNe_01_5gene", "Calculation2_mNe_02_5gene", "Calculation3_mNe_20_5gene", "Calculation4_mNe_21_5gene", "Calculation5_mNe_10_5gene", "Calculation6_mNe_12_5gene",
                         "Calculation1_mNe_01_0gene", "Calculation2_mNe_02_0gene", "Calculation3_mNe_21_0gene", "Calculation4_mNe_21_0gene", "Calculation5_mNe_10_0gene", "Calculation6_mNe_12_0gene")

results <- data.frame(Column = character(), `2.5%` = numeric(), `97.5%` = numeric(), Method = character(), Mean = numeric(), Median = numeric(), stringsAsFactors = FALSE)

for (col_name in columns_of_interest) {
  if (col_name %in% colnames(data)) {
    column_data <- data[[col_name]]
    if (is.numeric(column_data)) {
      conf_interval <- calculate_confidence_intervals(column_data)
      results <- rbind(results, data.frame(Column = col_name, `2.5%` = conf_interval$lower, `97.5%` = conf_interval$upper, Method = conf_interval$method, Mean = conf_interval$mean, Median = conf_interval$median, stringsAsFactors = FALSE))
    }
  }
}
# Save results to CSV file
csv_filename <- "~/Dropbox/Publication/GemHumMon/Data/Recent_boot_confidence_intervals_results.csv"
write.csv(results, csv_filename, row.names = FALSE)

# Display results
print(results)
