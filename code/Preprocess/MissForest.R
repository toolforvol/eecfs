# Use MissForest to impute missing values
library(missForest)
library(data.table)

# File path configuration
data_paths <- list(
  train_input = "../final_feature/train_2362_RNA141bp_final_feature_new_with_label.csv",
  test_input  = "../final_feature/test_238_RNA141bp_final_feature_new_with_label.csv",
  train_output = "../final_feature/train_2362_RNA141bp_final_feature_filled_data_new.csv",
  test_output  = "../final_feature/test_238_RNA141bp_final_feature_filled_data_new.csv"
)

# Monitor memory usage
memory_usage <- function() {
  mem <- sum(gc(reset = TRUE)[,2])  # Get current R session memory usage (MB)
  return(paste0(round(mem, 2), " MB"))
}

# Read data and convert to data frame format
load_data <- function(file_path) {
  data <- fread(file_path, data.table = FALSE)  # Use fread to read CSV for faster loading
  return(data)
}

# Perform missing value imputation using MissForest
missforest_imputation <- function(train_data, test_data) {
  cat("\nStart using MissForest to fill missing values...\n")
  start_time <- Sys.time()
  
  # Select data starting from column 6 for imputation
  train_data_to_fill <- train_data[, 6:ncol(train_data)]  # From column 6 to the last column
  test_data_to_fill <- test_data[, 6:ncol(test_data)]    # From column 6 to the last column
  
  # Train MissForest model
  cat("Filling train data...\n")
  train_filled <- missForest(train_data_to_fill, ntree = 10, maxiter = 10, verbose = TRUE)
  
  cat("Train data filled | Time:", round(difftime(Sys.time(), start_time, units = "secs"), 2), "s | Memory:", memory_usage(), "\n")
  
  # Fill test set using the same approach
  cat("Filling test data...\n")
  test_filled <- missForest(test_data_to_fill, ntree = 10, maxiter = 10, verbose = TRUE)
  
  cat("Test data filled | Total time:", round(difftime(Sys.time(), start_time, units = "secs"), 2), "s\n")
  
  # Merge imputed data back with original dataset (keep first 5 columns)
  train_filled_data <- cbind(train_data[, 1:5], train_filled$ximp)
  test_filled_data <- cbind(test_data[, 1:5], test_filled$ximp)
  
  return(list(train = train_filled_data, test = test_filled_data))  # Return imputed data
}

# Main function
main <- function() {
  cat("===== Start filling program =====\n")
  start_total <- Sys.time()
  
  # 1. Load data
  cat("Loading data...\n")
  train_data <- load_data(data_paths$train_input)
  test_data <- load_data(data_paths$test_input)
  
  cat("Data loaded. Train shape:", dim(train_data), "| Test shape:", dim(test_data), "\n")
  cat("Current memory usage:", memory_usage(), "\n")
  
  # 2. Missing value statistics
  cat("\nMissing value statistics:\n")
  cat("Train - Total missing:", sum(is.na(train_data)), sprintf("(%.2f%%)", mean(is.na(train_data)) * 100), "\n")
  cat("Test - Total missing:", sum(is.na(test_data)), sprintf("(%.2f%%)", mean(is.na(test_data)) * 100), "\n")
  
  # 3. Perform MissForest imputation
  filled_data <- missforest_imputation(train_data, test_data)
  
  # 4. Save imputation results
  cat("\nSaving filled results...\n")
  fwrite(filled_data$train, file = data_paths$train_output, row.names = FALSE)
  fwrite(filled_data$test, file = data_paths$test_output, row.names = FALSE)
  
  cat("\nTrain filled data saved to:", data_paths$train_output, "\n")
  cat("Test filled data saved to:", data_paths$test_output, "\n")
  cat("Total time:", round(difftime(Sys.time(), start_total, units = "secs"), 2), "seconds\n")
  cat("===== Finish =====\n")
}

# Run main program
main()