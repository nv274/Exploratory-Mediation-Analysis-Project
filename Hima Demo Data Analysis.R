library(broom)
library(ggplot2)
library(boot)
perform_bootstrap <- function(sourcefolder, x, y = x, func) {
  csv_files <- list.files(path = source_folder, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  boot_list <- list()
  boot_func <- function(data, indices) {
    return(func(data[indices]))
  }
  for (file in csv_files) {
    data <- read.csv(file)
    file_name <- tools::file_path_sans_ext(basename(file))
    if(y %in% colnames(data)) {
      print(file_name)
      boot_result = boot_tests(data, x, y, statistic = boot_func, R = 1000)
      boot_list[[file_name]] = boot_result
    }
  }
  return(boot_list)
}
boot_tests <- function(data, x, y, statistic, R = 1000) {
  results <- matrix(, nrow=length(unique(data[,x])), ncol = 2)
  rownames(results) <- unique(data[,x])
  colnames(results) <- c("2.5%", "97.5%")
  for(row in rownames(results)) {
    if(y == x) {boot_result <- boot(data[,row], statistic, R)}
    else {boot_result <- boot(data[data[,x] == row, y], statistic, R)}
    CI <- boot.ci(boot_result, type = "perc")
    results[row, "2.5%"] <- CI$percent[4]
    results[row, "97.5%"] <- CI$percent[5]
  }
  return(results)
}

kstests <- function(dat = E1classic) {
  ks <- matrix(, nrow=length(union(unique(dat$DominantEffect), "Direct")), ncol = 2)
  rownames(ks) <- union(unique(dat$DominantEffect), "Direct")
  colnames(ks) <- c("D_Statistic", "P_Value")
  for(row in rownames(ks)) {
    result <- ks.test(dat[,row], 'pnorm')
    ks[row, "D_Statistic"] <- result$statistic
    ks[row, "P_Value"] <- result$p.value
  }
  Direct <- ks.test(dat[dat$Treatment != 0, ]$Direct, 'pnorm')
  ks["Direct", "D_Statistic"] <- Direct$statistic
  ks["Direct", "P_Value"] <- Direct$p.value
  
  return(ks)
}
perform_ks_tests <- function(sourcefolder) {
  # List all CSV files in the source folder
  csv_files <- list.files(path = source_folder, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # Initialize a list to store the results
  ks_results_list <- list()
  # Loop through each CSV file
  for (file in csv_files) {
    # Read the CSV file
    data <- read.csv(file)
    
    # Extract the file name without the extension
    file_name <- tools::file_path_sans_ext(basename(file))
    print(file_name)
    # Create a matrix with the KS test results
    ks_matrix <- kstests(dat = data)
    # Add the matrix to the list
    ks_results_list[[file_name]] <- ks_matrix
  }
  # Combine all matrices into one large data frame
  # Return the data frame with KS test results
  return(ks_results_list)
}

perform_chisq_tests <- function(sourcefolder, x, y, exclude = "") {
  csv_files <- list.files(path = source_folder, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # Initialize a list to store the results
  chi_results_list <- list()
  # Loop through each CSV file
  for (file in csv_files) {
    # Read the CSV file
    data <- read.csv(file)
    
    # Extract the file name without the extension
    file_name <- tools::file_path_sans_ext(basename(file))
    if(length(unique(data[data[,x] != exclude,y])) < length(data[data[,x] != exclude,y]) & (y %in% colnames(data))) {
      print(file_name)
      chi_result <- chisq.test(data[data[,x] != exclude,y], data[data[,x] != exclude,x])
      chi_results_list[[file_name]] <- chi_result
    }
  }
  return(chi_results_list)
}

perform_aov_tests <- function(sourcefolder, x, y) {
  csv_files <- list.files(path = source_folder, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  # Initialize a list to store the results
  aov_results_list <- list()
  # Loop through each CSV file
  for (file in csv_files) {
    # Read the CSV file
    data <- read.csv(file)
    
    # Extract the file name without the extension
    file_name <- tools::file_path_sans_ext(basename(file))
    if(y %in% colnames(data)) {
      print(file_name)
      aov_result <- aov(y ~ x, data = data.frame(y = data[,y], x = data[,x]))
      aov_results_list[[file_name]] <- aov_result
    }
  }
  return(aov_results_list)
}

perform_tukey <- function(aov_tests) {
  tukey_list <- list()
  for(result in names(aov_tests)) {
    if(summary(aov_tests[[result]])[[1]][["Pr(>F)"]][1] <= 0.05) {
      print(result)
      tukey_list[[result]] <- TukeyHSD(aov_tests[[result]], conf.level = 0.95)
    } 
  }
  return(tukey_list)
}

save_tukeys <- function(tukey_list, value, group) {
  csv_files <- list.files(path = getwd(), pattern = "\\.csv$", recursive = TRUE, full.names = FALSE)
  if (length(tukey_list) > 0) {
    names <- names(tukey_list)
    for (file in csv_files) {
      file_name <- tools::file_path_sans_ext(basename(file))
      if(file_name %in% names) {
        plotname <- paste(file_name,value,"~",group,"Pairwise 95% CI Plot")
        filename <- paste0(dirname(file),"plots/",file_name,value,"~",group,"Tukey.png")
        tidy_tukey <- broom::tidy(tukey_list[[file_name]])
        pl <- ggplot(tidy_tukey, aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high)) +
          geom_pointrange() +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          coord_flip() +
          labs(
            title = plotname,
            x = paste(group,"Comparison"),
            y = "Difference in Means"
          ) +
          theme_minimal()
        ggsave(filename = filename, plot = pl,  width = 8, height = 8, dpi = 300)
      }
    }
  }
}

save_violins <- function(tukey_list, group, value) {
  csv_files <- list.files(path = getwd(), pattern = "\\.csv$", recursive = TRUE, full.names = FALSE)
  if (length(tukey_list) > 0) {
    names <- names(tukey_list)
    for (file in csv_files) {
      file_name <- tools::file_path_sans_ext(basename(file))
      if(file_name %in% names) {
        data <- read.csv(file)
        if(group %in% colnames(data)) {
          plotname <- paste(file_name,value,"~",group,"Violin Plot")
          filename <- paste0(dirname(file),"plots/",file_name,value,"~",group,"Violin.png")
          pl <- ggplot(data, aes_string(x = group, y = value, fill = group)) +
            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),trim = FALSE) +
            scale_fill_brewer(palette = "Set3") + # Use a predefined color palette
            labs(title = plotname, x = group, y = value) +
            theme_minimal()
          ggsave(filename = filename, plot = pl,  width = 8, height = 8, dpi = 300)
        }
      }
    }
  }
}

source_folder <- getwd()
# Example usage
# Assuming the CSV files are located in a folder named "data_folder"
ks_results <- perform_ks_tests(sourcefolder = source_folder)
chi_results_sex <- perform_chisq_tests(sourcefolder = source_folder, x = "DominantEffect", y = "Sex")
chi_results_treatment <- perform_chisq_tests(sourcefolder = source_folder
                                             , x = "DominantEffect"
                                             , y = "Treatment"
                                             , exclude = "Direct")

chi_results_E2 <- perform_chisq_tests(sourcefolder = source_folder
                                      , x = "DominantEffect"
                                      , y = "Disease")

aov_results_outcome <- perform_aov_tests(source_folder, x = "DominantEffect", y = "Outcome")
aov_results_E2 <- perform_aov_tests(source_folder, x = "DominantEffect", y = "totalLogOdds")
aov_results_E3 <- perform_aov_tests(source_folder, x = "DominantEffect", y = "totalHazardRatio")
aov_results_age <- perform_aov_tests(source_folder, x = "DominantEffect", y = "Age")
aov_results_time <- perform_aov_tests(source_folder, x = "DominantEffect", y = "Time")

tukey_outcome <- perform_tukey(aov_results_outcome)
tukey_age <- perform_tukey(aov_results_age)
tukey_e2 <- perform_tukey(aov_results_E2)
tukey_e3 <- perform_tukey(aov_results_E3)
tukey_time <- perform_tukey(aov_results_time)

save_violins(tukey_outcome, "DominantEffect", "Outcome")
save_violins(tukey_age, "DominantEffect", "Age")
save_violins(tukey_e2, "DominantEffect", "totalLogOdds")
save_violins(tukey_e3, "DominantEffect", "totalHazardRatio")
save_violins(tukey_time, "DominantEffect", "Time")

save_tukeys(tukey_outcome, "Outcome", "DominantEffect")
save_tukeys(tukey_age, "Age", "DominantEffect")
save_tukeys(tukey_e2, "totalLogOdds", "DominantEffect")
save_tukeys(tukey_e3, "totalHazardRatio", "DominantEffect")
save_tukeys(tukey_time, "Time", "DominantEffect")

outcome_CI <- perform_bootstrap(source_folder, "DominantEffect", "Outcome", mean)
age_CI <- perform_bootstrap(source_folder, "DominantEffect", "Age", mean)
logodds_CI <- perform_bootstrap(source_folder, "DominantEffect", "totalLogOdds", mean)
hazard_CI <- perform_bootstrap(source_folder, "DominantEffect", "totalHazardRatio", mean)
time_CI <- perform_bootstrap(source_folder, "DominantEffect", "Time", mean)
effects_CI <- perform_bootstrap(source_folder, x = "DominantEffect", func = mean)
