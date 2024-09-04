# Load necessary libraries
library(tidyverse)
library(rstatix)
library(permuco)
library(emmeans)

# Set paths 
path <- file.path("/Users/sara/Desktop/Projects/ExternalCollaborations/1_Fingerprinting_AD/repo/outputs/") # change this accordingly to where your outputs path is located 
setwd(path)

process_dataset <- function(prefix) {
  # Read the dataset
  data_iself <- read.csv(paste0("Iself_", prefix, ".csv"))
  data_iothers <- read.csv(paste0("Iothers_", prefix, ".csv"))
  data_iself$Group <- as.factor(data_iself$Group)
  data_iothers$Group <- as.factor(data_iothers$Group)
  
  # Outliers for Iself
  iself_outliers <- data_iself %>%
    group_by(Group) %>%
    identify_outliers(Iself)
  
  # Check Normality assumption by Groups for Iself
  iself_normality <- data_iself %>%
    group_by(Group) %>%
    shapiro_test(Iself)
  
  # Permutation ANOVA for Iself
  aov_iself <- aovperm(Iself ~ Group + Age + Sex + YoE + delta_FD, data = data_iself, np = 5000)
  
  # Outliers for Iothers
  iothers_outliers <- data_iothers %>%
    group_by(Group) %>%
    identify_outliers(Iothers)
  
  # Check Normality assumption by Groups for Iothers
  iothers_normality <- data_iothers %>%
    group_by(Group) %>%
    shapiro_test(Iothers)
  
  # Conditional ANOVA formula for Iothers based on prefix
  if (prefix == "ADNI") {
    # Permutation ANOVA for Iothers including ScannerType for ADNI
    aov_iothers <- aovperm(Iothers ~ Group + Age + Sex + YoE + FD + Scanner, data = data_iothers, np = 5000)
  } else {
    # Permutation ANOVA for Iothers without ScannerType for non-ADNI datasets
    aov_iothers <- aovperm(Iothers ~ Group + Age + Sex + YoE + FD, data = data_iothers, np = 5000)
  }
  
  # Post-hoc pairwise comparisons for Iothers
  pwc <- data_iothers %>%
    emmeans_test(
      Iothers ~ Group, covariate = FD,
      p.adjust.method = "bonferroni"
    )
  
  return(list(
    iself_outliers = iself_outliers, 
    iself_normality_test = iself_normality, 
    aov_iself = aov_iself, 
    iothers_outliers = iothers_outliers, 
    iothers_normality_test = iothers_normality, 
    aov_iothers = aov_iothers, 
    pwc = pwc
  ))
}

# Dataset prefixes
datasets <- c("GENEVA", "ADNI")

# Loop through the datasets
results <- lapply(datasets, process_dataset)

# Extract results for each dataset
results_geneva <- results[[1]]
results_adni <- results[[2]]

# Print results for verification
print("Results for GENEVA:")
print(results_geneva)

print("Results for ADNI:")
print(results_adni)
