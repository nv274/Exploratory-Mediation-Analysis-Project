functions.pck <-
c("functions.pck", "survival", "ggplot2", "broom")
perform_bootstrap <- function(sourcefolder, x, y, func) {
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
perform_cor_tests <- function(data, x, y = x, alternative, method) {
	cor_results <- list()
	if (setequal(x,y)) {
		for(i in 1:(length(x)-1)) {
			for(j in (i+1):length(x)) {
				cor_results[[paste(x[i],"-",x[j])]] <- cor.test(data[,x[i]], data[,x[j]], alternative, method)
			}
		}
	}
	else {
		for(name in x) {
			cor_results[[paste(name,"-", y)]] <- cor.test(data[,name], data[,y], alternative, method)
		}
	}
	return(cor_results)
}
perform_kruskal <- function(data, x, y) {
	results <- list()
	for(name in y) {
		results[[name]] <- kruskal.test(y_val ~ x_val, data = data.frame(y_val = as.vector(data[,name]), x_val = as.vector(data[,x])))
		print(name)
	}
	return(results)
}
boot_tests <- function(data, x, y = x, statistic, R = 1000) {
  boot_statistic <- function(data, indices) {
    return(statistic(data[indices]))
  }
  results <- matrix(, nrow = length(unique(data[,x])), ncol = 2)
  rownames(results) <- unique(data[,x])
  colnames(results) <- c("2.5%", "97.5%")
  for(row in rownames(results)) {
    if(y == x) {boot_result <- boot(data[,row], boot_statistic, R)}
    else {boot_result <- boot(data[data[,x] == row, y], boot_statistic, R)}
    CI <- boot.ci(boot_result, type = "perc")
    results[row, "2.5%"] <- CI$percent[4]
    results[row, "97.5%"] <- CI$percent[5]
  }
  return(results)
}

kstests <- function(dat = E1classic, measured = colnames(dat), exclude) {
  ks <- matrix(, nrow=length(setdiff(measured, exclude)), ncol = 2)
  rownames(ks) <- setdiff(measured, exclude)
  print(rownames(ks))
  colnames(ks) <- c("D_Statistic", "P_Value")
  for(row in rownames(ks)) {
    result <- ks.test(dat[,row], 'pnorm')
    ks[row, "D_Statistic"] <- result$statistic
    ks[row, "P_Value"] <- result$p.value
  }
  
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
          tidy_tukey <- broom::tukey_list[[name]]
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
plot_violins <- function(dat, x, y, plotname) {
	pl <- ggplot(dat, aes_string(x = x, y = y, fill = x)) +
            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),trim = FALSE) +
            scale_fill_brewer(palette = "Set3") + # Use a predefined color palette
            labs(title = plotname, x = x, y = y) +
            theme_minimal()
	return(pl)

}
estimatedEffects <- function(fit = hima.fit,
                             betas = fit[,"beta_hat"],
                             IDEcoefs = fit[,"IDE"],
                             type = "continuous",
                             dat = himaDat$Example1,
                             X = himaDat$Example1$PhenoData$Treatment,
                             Y = himaDat$Example1$PhenoData$Outcome,
                             Cov = himaDat$Example1$PhenoData[, c("Sex", "Age")]) {
  
  # Start with base data frame
  result <- data.frame(dat$PhenoData)
  # Multiply mediator values by their respective coefficients
  mediator_indices <- fit[,1]
  result[, mediator_indices] <- sweep(dat$Mediator[, mediator_indices, drop = FALSE], 2, betas, FUN = `*`)
  IDEs <- (sum(IDEcoefs)*X)
  # Create interaction terms dynamically with all covariates in Cov
  interaction_terms <- paste("X +", colnames(Cov), collapse = " + ")
  target_cols <- c(mediator_indices, "Direct")
  formula <- as.formula(paste("Y ~ ", interaction_terms))
  if (type == "binary") {
    model <- glm(formula, data = data.frame(Y = Y, X = X, Cov), family = "binomial")
    result$totalLogOdds <- predict(model, newdata = data.frame(Y = Y, X = X, Cov), type = "link")
  } else {
    # Direct effect with continuous outcome
    # Fit the linear model using constructed formula and provided data
    if(is.Surv(Y)) {
      model <- coxph(formula, data = data.frame(Y = Y, X = X, Cov))
      result$totalHazardRatio <- predict(model)
      Y <- result$totalHazardRatio
    }
    else {
      model <- lm(formula, data = data.frame(Y = Y, X = X, Cov))
    }
  }
  individual_X_contribution <- predict(model, newdata = data.frame(Y = Y, X = X, Cov)) - (IDEs)
  result$Direct <- ifelse(result$Treatment != 0, individual_X_contribution, 0)
  result$Residuals <- Y - (rowSums(result[, mediator_indices]) + result$Direct)
  result$Adjusted <- rowSums(result[,c(mediator_indices, "Direct")])
  result$DominantEffect <- apply(result[, target_cols], 1, function(row) {
    target_cols[which.max(abs(row))]
  })
  return(result)
}
plot_examples <- function(data, plotname) {
  fullmodel <- lm(Outcome ~ Treatment + Sex + Age, data = data)
  mse <- mean((residuals(fullmodel))^2)
  mae <- mean(abs(residuals(fullmodel)))
  fullmodelpl <- ggplot(data.frame(Predicted = predict(fullmodel),
                                           data),
                                aes(x = Predicted, y = Outcome)) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = plotname, x = "Predicted", y = "Actual", color = "\n") +
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, label = paste("MSE:", round(mse, 2), "\nMAE:", round(mae, 2)),
             hjust = 1.1, vjust = 1.1, size = 4, color = "red")
  return(fullmodelpl)
}

plot_mediator_adjustments <- function(type = "continuous",
                          actual,
                          estimations,
                          mediator_indices,
                          plotname,
                          labels = NULL) {
  predicted <- rowSums(estimations[,mediator_indices]) + estimations$Direct
  regression_method = "lm"
  x_label = "Predicted"
  y_label = "Actual"
  pl <- NULL
  if(type == "binary") {
    x_label = "Treatment"
    y_label = "Probability"
    probs <- 1/(exp(-predicted)+1) 
    predoutcome <- ifelse(probs >= 0.5, 1, 0)
    logloss <- -mean(actual * log(probs) + (1 - actual) * log(1 - probs))
    cm <- table(predoutcome, actual)
    if (nrow(cm) == 1) {
      cm <- rbind(cm, c(0,0))
      rownames(cm)[2] <- 1-as.numeric(rownames(cm)[1])
    }
    print(cm)
    mcc <- ((cm["1","1"] 
               * cm["0","0"] 
               - cm["0","1"] 
               * cm["1","0"]) 
              / sqrt((cm["1","1"] 
                      + cm["1","0"]) 
                     * (cm["1","1"] 
                        + cm["0","1"]) 
                     * (cm["0","0"] 
                        + cm["1","0"]) 
                     * (cm["0","0"] 
                        + cm["0","1"])))
    print(data.frame(mcc))
    pl <- ggplot(data.frame(estimations,y_vals = probs),
                 aes(x = Treatment,
                     y = y_vals)) +
      geom_point(position = position_jitter(width=0.06, height = 0.045),aes(color = factor(actual)), size = 2) + # Use 'factor' to treat 'value' as categorical
      geom_hline(yintercept=0.5, linetype='dotted', linewidth = 1, col = "blue") +
      labs(title = plotname, x = x_label, y = y_label, color = labels[1]) +
      scale_color_manual(labels = labels[c(2,3)],values = c("0" = "red", "1" = "green")) + # Manually set colors  # Points colored by actual value
      theme_minimal() +
      annotate("text", x = 0.6, y = 1, label = paste("Log Loss:", round(logloss, 2), "\nMCC:", round(mcc, 2), "\nTP:", cm["1","1"],"\nTN:", cm["0","0"], "\nFP:", cm["1","0"], "\nFN:", cm["0","1"]),
               hjust = 1.1, vjust = 1.1, size = 4, color = "red")
  
  }
  else {
    if(is.Surv(actual)) {
      x_label = "Hazard"
      y_label = "Time"
      pl <- ggplot(data.frame(x_vals = predicted,
                              estimations),
                   aes(x = x_vals, y = Time)) +
        geom_point(color = "black", size = 2) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = plotname, x = "Hazard Ratio", y = "Time", color = "Disease\n") +
        theme_minimal()
    }
    else {
      mse <- mean((actual-predicted)^2)
      mae <- mean(abs(actual-predicted))
      pl <- ggplot(data.frame(Predicted = predicted,
                              Actual = actual),
                   aes(x = Predicted, y = Actual)) +
        geom_point(color = "black", size = 2) +
        geom_smooth(method = regression_method, se = FALSE, color = "blue") +
        labs(title = plotname, x = x_label, y = y_label, color = "\n") +
        theme_minimal() +
        annotate("text", x = Inf, y = Inf, label = paste("MSE:", round(mse, 2), "\nMAE:", round(mae, 2)),
                 hjust = 1.1, vjust = 1.1, size = 4, color = "red")
    }
  }
  return(pl)
}

do_tests <- function(data, y,  outcome_type = "continuous") {
	ks_result <- kstests(data)
	y_col <- data[,y]
	results <- list()
	is_significant <- function(p) {
		if( p <= 0.05) return(TRUE)
		else return(FALSE)
	}
	results[["Kolmogorov-Smirnov Test"]] <- ks_result
	is_normal <- all(sapply(ks_result[,"P_Value"], is_significant))
	if(is.Surv(y_col)) {
		y_col <- y_col[,"time"]
		hazard_CI <- boot_tests(data, "DominantEffect", "totalHazardRatio", mean)
		results[["95% CI Total Hazard Ratios ~ Dominant Effect"]] <- hazard_CI
		aov_hazard_ratio <- aov(totalHazardRatio ~ DominantEffect, data = data.frame(data))
		results[["ANOVA: Hazard Ratio ~ Dominant Effect"]] <- summary(aov_hazard)
	}
	if(is_normal & outcome_type != "binary") {
		aov_y <- aov(y_val ~ DominantEffect, data = data.frame(data, y_val = y_col))
		results[[paste("ANOVA:", y,"~ Dominant Effect")]] <- summary(aov_y)
		if(is_significant(summary(aov_y)[[1]][["Pr(>F)"]][1])) {
			tukey_y = TukeyHSD(aov_y, conf.level = 0.95)
			results[[paste("Tukey HSD:", y,"~ Dominant Effect")]] <- tukey_y
		}
		y_CI <- boot_tests(data, "DominantEffect", y, mean)
	}
	else {
		chisq_y <- chisq.test(y_col, data$DominantEffect)
		results[[paste("Chi-squared:",y,"~ Dominant Effect")]] <- list(chisq_y, chisq_y$observed)
		logodds_CI <- boot_tests(data, "DominantEffect", "totalLogOdds", mean)
		results[["95% CI Total Log Odds ~ Dominant Effect"]] <- logodds_CI
		aov_logodds <- aov(totalLogOdds ~ DominantEffect, data = data.frame(data))
		results[["ANOVA: Log Odds ~ Dominant Effect"]] <- summary(aov_logodds)

	}
	results[[paste("95% CI",y,"~ Dominant Effect")]] <- y_CI
	age_CI <- boot_tests(data, "DominantEffect", "Age", mean)
	results[[("95% CI Age~ Dominant Effect")]] <- age_CI
	aov_age <- aov(Age ~ DominantEffect, data = data.frame(data))
	results[["ANOVA: Age ~ Dominant Effect"]] <- summary(aov_age)
	if(is_significant(summary(aov_y)[[1]][["Pr(>F)"]][1])) {
		tukey_age = TukeyHSD(aov_age, conf.level = 0.95)
		results[["Tukey HSD: Age ~ Dominant Effect"]] <- tukey_age
	}
	chisq_sex <- chisq.test(data$Sex, data$DominantEffect)
	results[["Chi-squared: Sex ~ DominantEffect"]] <- list(chisq_sex, matrix(chisq_sex$observed))
	chisq_treatment <- chisq.test(data[data$DominantEffect != "Direct",]$Treatment, data[data$DominantEffect != "Direct",]$DominantEffect)
	results[["Chi-squared: Treatment ~ DominantEffect"]] <- list(chisq_treatment, chisq_treatment$observed)
	return(results)
}
