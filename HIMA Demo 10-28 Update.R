library(HIMA)
library(survival)
library(ggplot2)
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
  IDEs <- sum(IDEcoefs)*X
  # Create interaction terms dynamically with all covariates in Cov
  interaction_terms <- paste("X +", colnames(Cov), collapse = " + ")
  target_cols <- c(mediator_indices, "Direct")
  formula <- as.formula(paste("Y ~ ", interaction_terms))
  if (type == "binary") {
    model <- glm(formula, data = data.frame(Y = Y, X = X, Cov), family = "binomial")
    result$totalLogOdds <- predict(model, newdata = data.frame(Y = Y, X = X, Cov), type = "link") - coef(model)["(Intercept)"]
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
  result$Other <- Y - (rowSums(result[, mediator_indices]) + result$Direct)
  result$DominantEffect <- apply(result[, target_cols], 1, function(row) {
    target_cols[which.max(abs(row))]
  })
  return(result)
}
plot_and_save_full_examples <- function(data, plotname, filename) {
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
  ggsave(filename = filename, plot = fullmodelpl,  width = 8, height = 8, dpi = 300)
  return(fullmodelpl)
}

plot_and_save <- function(type = "continuous",
                          actual,
                          estimations,
                          mediator_indices,
                          filename,
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
    pl <- ggplot(data.frame(estimations,y_vals = probs),
                 aes(x = Treatment,
                     y = y_vals)) +
      geom_point(position = position_jitter(width=0.06, height = 0.045),aes(color = factor(actual)), size = 2) + # Use 'factor' to treat 'value' as categorical
      geom_hline(yintercept=0.5, linetype='dotted', linewidth = 1, col = "blue") +
      labs(title = plotname, x = x_label, y = y_label, color = labels[1]) +
      scale_color_manual(labels = labels[c(2,3)],values = c("0" = "red", "1" = "green")) + # Manually set colors  # Points colored by actual value
      theme_minimal() +
      annotate("text", x = 0.6, y = 1, label = paste("Log Loss:", round(logloss, 2), "\nMCC:", round(mcc, 2)),
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
  ggsave(filename = filename, plot = pl,  width = 8, height = 8, dpi = 300)
  return(pl)
}
## Not run:
# Note: In the following examples, M1, M2, and M3 are true mediators.
data(himaDat)
# When Y is continuous and normally distributed
# Example 1 (continuous outcome):
head(himaDat$Example1$PhenoData)
Example1fullmodelpl <- plot_and_save_full_examples(data = himaDat$Example1$PhenoData,
                                                   plotname = "Example One Full Model Outcome Predictions",
                                                   filename = "Example1plots/Example1fullmodel.png")
Example1fullmodelpl
hima.fit <- classicHIMA(
  X = himaDat$Example1$PhenoData$Treatment,
  Y = himaDat$Example1$PhenoData$Outcome,
  M = himaDat$Example1$Mediator,
  COV.XM = himaDat$Example1$PhenoData[, c("Sex", "Age")],
  Y.type = "continuous",
  scale = FALSE, # Disabled only for simulation data
  verbose = TRUE
)
hima.fit

E1classic <- estimatedEffects()
Example1classicmodelpl <- plot_and_save(actual = himaDat$Example1$PhenoData$Outcome,
                                        estimations = E1classic,
                                        mediator_indices = hima.fit[,"Index"],
                                        filename = "Example1plots/Example1classicmodel.png",
                                        plotname = "Example One Classic HIMA Outcome Predictions")
Example1classicmodelpl

# When Y is binary
# Example 2 (binary outcome):
head(himaDat$Example2$PhenoData)
E2pred <- predict(Example2fullmodel, type = "response")
E2predoutcome <- ifelse(E2pred >= 0.5, 1, 0)
E2logloss <- -mean(himaDat$Example2$PhenoData$Disease * log(E2pred) + (1 - himaDat$Example2$PhenoData$Disease) * log(1 - E2pred))
E2cm <- table(E2predoutcome, himaDat$Example2$PhenoData$Disease)
if (nrow(E2cm) == 1) {
  E2cm <- rbind(E2cm, c(0,0))
  rownames(E2cm)[2] <- 1-as.numeric(rownames(E2cm)[1])
}
E2mcc <- ((E2cm["1","1"] 
          * E2cm["0","0"] 
          - E2cm["0","1"] 
          * E2cm["1","0"]) 
          / sqrt((E2cm["1","1"] 
                  + E2cm["1","0"]) 
                 * (E2cm["1","1"] 
                    + E2cm["0","1"]) 
                 * (E2cm["0","0"] 
                    + E2cm["1","0"]) 
                 * (E2cm["0","0"] 
                    + E2cm["0","1"])))
Example2fullmodelpl <- ggplot(data.frame(himaDat$Example2$PhenoData, Response = predict(Example2fullmodel, type = "response")),
                              aes(x = Treatment,
                                  y = Response)) +
  geom_point(position = position_jitter(width=0.06, height = 0.045),aes(color = factor(Disease)), size = 2) + # Use 'factor' to treat 'value' as categorical
  geom_hline(yintercept=0.5, linetype='dotted', linewidth = 1, col = "blue") +
  labs(title = "Example Two Full Model Outcome Predictions", x = "Treatment", y = "Probability of Testing Positive", color = "Disease\n") +
  scale_color_manual(labels = c("Negative","Positive"),values = c("0" = "red", "1" = "green")) + # Manually set colors  # Points colored by actual value
  theme_minimal() +
  annotate("text", x = 0.6, y = 1, label = paste("Log Loss:", round(E2logloss, 2), "\nMCC:", round(E2mcc, 2)),
           hjust = 1.1, vjust = 1.1, size = 4, color = "red")
Example2fullmodelpl
ggsave(filename = "Example2plots/Example2fullmodel.png", plot = Example2fullmodelpl, width = 8, height = 8, dpi = 300)

hima.logistic.fit <- classicHIMA(
  X = himaDat$Example2$PhenoData$Treatment,
  Y = himaDat$Example2$PhenoData$Disease,
  M = himaDat$Example2$Mediator,
  COV.XM = himaDat$Example2$PhenoData[, c("Sex", "Age")],
  Y.type = "binary",
  scale = FALSE, # Disabled only for simulation data
  verbose = TRUE
)
hima.logistic.fit



E2classic <- estimatedEffects(fit = hima.logistic.fit,
                              type = "binary",
                              dat = himaDat$Example2,
                              X = himaDat$Example2$PhenoData$Treatment,
                              Y = himaDat$Example2$PhenoData$Disease,
                              Cov = himaDat$Example2$PhenoData[, c("Sex", "Age")])


E2classicmodelpl <- plot_and_save(type = "binary",
              actual = himaDat$Example2$PhenoData$Disease,
              estimations = E2classic,
              mediator_indices = hima.logistic.fit[,"Index"],
              filename = "Example2plots/Example2classicmodel.png",
              plotname = "Example Two Classic HIMA Outcome Predictions",
              labels = c("Disease", "Negative", "Positive"))
E2classicmodelpl
## Not run:
# Note: In the following examples, M1, M2, and M3 are true mediators.
# Y is continuous and normally distributed
# Example:
head(himaDat$Example1$PhenoData)
dblassohima.fit <- dblassoHIMA(
  X = himaDat$Example1$PhenoData$Treatment,
  Y = himaDat$Example1$PhenoData$Outcome,
  M = himaDat$Example1$Mediator,
  COV = himaDat$Example1$PhenoData[, c("Sex", "Age")],
  scale = FALSE, # Disabled only for simulation data
  FDRcut = 0.05,
  verbose = TRUE
)
dblassohima.fit
## End(Not run)
E1dblasso <- estimatedEffects(fit = dblassohima.fit)
E1dblassopl <- plot_and_save(actual = himaDat$Example1$PhenoData$Outcome,
                             estimations = E1dblasso,
                             mediator_indices = dblassohima.fit[,"Index"],
                             filename = "Example1plots/Example1dblassomodel.png",
                             plotname = "Example One Debiased Lasso Penalty HIMA Outcome Predictions")
E1dblassopl
## Not run:http://127.0.0.1:9731/graphics/c88e4117-1df4-49c8-834e-59864187f513.png
# Note: In the following example, M1, M2, and M3 are true mediators.
# Y is continuous and normally distributed
# Example (continuous outcome):
head(himaDat$Example1$PhenoData)
eHIMA.fit <- eHIMA(
  X = himaDat$Example1$PhenoData$Treatment,
  Y = himaDat$Example1$PhenoData$Outcome,
  M = himaDat$Example1$Mediator,
  COV = himaDat$Example1$PhenoData[, c("Sex", "Age")],
  scale = FALSE, # Disabled only for simulation data
  FDRcut = 0.05,
  verbose = TRUE
)
eHIMA.fit
## End(Not run)
E1efficient <- estimatedEffects(fit = eHIMA.fit)
E1efficientpl <- plot_and_save(actual = himaDat$Example1$PhenoData$Outcome,
                             estimations = E1efficient,
                             mediator_indices = eHIMA.fit[,"Index"],
                             filename = "Example1plots/Example1plotsExample1efficientmodel.png",
                             plotname = "Example One Efficient HIMA Outcome Predictions")
E1efficientpl

# Example 1 (continuous outcome - linear HIMA):
head(himaDat$Example1$PhenoData)
e1 <- himaFit(Outcome ~ Treatment + Sex + Age,
              data.pheno = himaDat$Example1$PhenoData,
              data.M = himaDat$Example1$Mediator,
              mediator.type = "gaussian",
              penalty = "MCP", # Can be "DBlasso" for dblassoHIMA
              scale = FALSE
) # Disabled only for simulation data
e1
E1mcp <- estimatedEffects(fit = e1, IDEcoefs = e1[,"alpha*beta"], betas = e1[,"beta"])
attributes(e1)$variable.labels
# Efficient HIMA (only applicable to mediators and outcomes that are
# both continuous and normally distributed.)
E1mcppl <- plot_and_save(actual = himaDat$Example1$PhenoData$Outcome,
                               estimations = E1mcp,
                               mediator_indices = e1[,"ID"],
                               filename = "Example1plots/Example1mcpmodel.png",
                               plotname = "Example One MCP HIMA Outcome Predictions")
E1mcppl
e1e <- himaFit(Outcome ~ Treatment + Sex + Age,
               data.pheno = himaDat$Example1$PhenoData,
               data.M = himaDat$Example1$Mediator,
               mediator.type = "gaussian",
               efficient = TRUE,
               penalty = "MCP", # Efficient HIMA does not support DBlasso
               scale = FALSE
) # Disabled only for simulation data
e1e
E1mcpEfficient<- estimatedEffects(fit = e1e, IDEcoefs = e1e[,"alpha*beta"], betas = e1e[,"beta"])
E1mcpefficientpl <- plot_and_save(actual = himaDat$Example1$PhenoData$Outcome,
                         estimations = E1mcpEfficient,
                         mediator_indices = e1e[,"ID"],
                         filename = "Example1plots/Example1mcpefficientmodel.png",
                         plotname = "Example One Efficient MCP HIMA Outcome Predictions")
E1mcpefficientpl
attributes(e1e)$variable.labels
# Example 2 (binary outcome - logistic HIMA):
head(himaDat$Example2$PhenoData)
e2 <- himaFit(Disease ~ Treatment + Sex + Age,
              data.pheno = himaDat$Example2$PhenoData,
              data.M = himaDat$Example2$Mediator,
              mediator.type = "gaussian",
              penalty = "MCP",
              scale = FALSE
) # Disabled only for simulation data
e2

E2mcp <- estimatedEffects(fit = e2,
                          betas = e2[,"beta"],
                          IDEcoefs = e2[,"alpha*beta"],
                          type = "binary",
                          dat = himaDat$Example2,
                          X = himaDat$Example2$PhenoData$Treatment,
                          Y = himaDat$Example2$PhenoData$Disease,
                          Cov = himaDat$Example2$PhenoData[, c("Sex", "Age")])

attributes(e2)$variable.labels

E2mcpmodelpl <- plot_and_save(type = "binary",
                                  actual = himaDat$Example2$PhenoData$Disease,
                                  estimations = E2mcp,
                                  mediator_indices = e2[,"ID"],
                                  filename = "Example2plots/Example2mcpmodel.png",
                                  plotname = "Example Two MCP HIMA Outcome Predictions",
                                  labels = c("Disease", "Negative", "Positive"))
E2mcpmodelpl

# Example 3 (time-to-event outcome - survival HIMA):
head(himaDat$Example3$PhenoData)
Example3fullmodel <- coxph(Surv(Time, Status) ~ Treatment + Sex + Age, data = himaDat$Example3$PhenoData)
Example3fullmodelpl <- ggplot(data.frame(Predicted = predict(Example3fullmodel),
                                         himaDat$Example3$PhenoData),
                              aes(x = Predicted, y = Time)) +
  geom_point(color = "black", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Example Three Full Model Outcome Predictions", x = "Hazard Ratio", y = "Time", color = "Disease\n") +
  theme_minimal()
Example3fullmodelpl
ggsave(filename = "Example3plots/Example3fullmodel.png", plot = Example3fullmodelpl,  width = 8, height = 8, dpi = 300)
e3 <- himaFit(Surv(Status, Time) ~ Treatment + Sex + Age,
              data.pheno = himaDat$Example3$PhenoData,
              data.M = himaDat$Example3$Mediator,
              mediator.type = "gaussian",
              penalty = "DBlasso",
              scale = FALSE
) # Disabled only for simulation data
e3
E3dblasso <- estimatedEffects(fit = e3,
                              betas = e3[,"beta"],
                              IDEcoefs = e3[,"alpha*beta"],
                              type = "continuous",
                              dat = himaDat$Example3,
                              X = himaDat$Example3$PhenoData$Treatment,
                              Y = Surv(himaDat$Example3$PhenoData$Time, himaDat$Example3$PhenoData$Status),
                              Cov = himaDat$Example3$PhenoData[, c("Sex", "Age")])

E3dblassomodelpl <- plot_and_save(actual = Surv(himaDat$Example3$PhenoData$Time, himaDat$Example3$PhenoData$Status),
                                  estimations = E3dblasso,
                                  mediator_indices = e3[,"ID"],
                                  filename = "Example3plots/Example3dblassomodel.png",
                                  plotname = "Example Three Debiased Lasso HIMA Survival Time Predictions")
E3dblassomodelpl
attributes(e3)$variable.labels

# Example 4 (compositional data as mediator, e.g., microbiome):
head(himaDat$Example4$PhenoData)
Example4fullmodelpl <- plot_and_save_full_examples(data = himaDat$Example4$PhenoData,
                                                   plotname = "Example Four Full Model Outcome Predictions",
                                                   filename = "Example4plots/Example4fullmodel.png")
Example4fullmodelpl
e4 <- himaFit(Outcome ~ Treatment + Sex + Age,
              data.pheno = himaDat$Example4$PhenoData,
              data.M = himaDat$Example4$Mediator,
              mediator.type = "compositional",
              penalty = "DBlasso"
) # Scaling is always enabled internally for microHIMA
e4
E4dblasso <- estimatedEffects(fit = e4,
                              betas = e4[,"beta"],
                              IDEcoefs = e4[,"alpha*beta"],
                              dat = himaDat$Example4,
                              X = himaDat$Example4$PhenoData$Treatment,
                              Y = himaDat$Example4$PhenoData$Outcome,
                              Cov = himaDat$Example4$PhenoData[, c("Sex", "Age")])
E4dblassopl <- plot_and_save(actual = himaDat$Example4$PhenoData$Outcome,
                             estimations = E4dblasso,
                             mediator_indices = e4[,"ID"],
                             filename = "Example4plots/Example4dblassomodel.png",
                             plotname = "Example Four Debiased Lasso Penalty HIMA Outcome Predictions")
E4dblassopl
attributes(e4)$variable.labels
# Example 5 (quantile mediation anlaysis - quantile HIMA):
head(himaDat$Example5$PhenoData)
# Note that the function will prompt input for quantile level.
Example5_30thfullmodelpl <- plot_and_save_full_examples(data = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.3),],
                                                   plotname = "Example Five 30th Percentile Full Model Outcome Predictions",
                                                   filename = "Example5plots/Example5_30thfullmodel.png")
Example5_30thfullmodelpl
Example5_50thfullmodelpl <- plot_and_save_full_examples(data = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.5),],
                                                        plotname = "Example Five 50th Percentile Full Model Outcome Predictions",
                                                        filename = "Example5plots/Example5_50thfullmodel.png")
Example5_50thfullmodelpl
Example5_70thfullmodelpl <- plot_and_save_full_examples(data = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.7),],
                                                        plotname = "Example Five 70th Percentile Full Model Outcome Predictions",
                                                        filename = "Example5plots/Example5_70thfullmodel.png")
Example5_70thfullmodelpl
e5 <- himaFit(Outcome ~ Treatment + Sex + Age,
              data.pheno = himaDat$Example5$PhenoData,
              data.M = himaDat$Example5$Mediator,
              mediator.type = "gaussian",
              quantile = TRUE,
              penalty = "MCP", # Quantile HIMA does not support DBlasso
              scale = FALSE, # Disabled only for simulation data
              tau = c(0.3, 0.5, 0.7)
) # Specify multiple quantile level
e5
max_rimp <- aggregate(`Relative Importance (%)` ~ tau + ID, e5, max)
e5 <- merge(e5, max_rimp)
e5 <- e5[order(e5$tau),c("ID", "alpha", "beta", "alpha*beta", "Relative Importance (%)", "p-value", "tau")]

attributes(e5)$variable.labels
## End(Not run)
E5MCP30th <- estimatedEffects(fit = e5[e5$tau == 0.3,],
                              IDEcoefs = e5[e5$tau==0.3,"alpha*beta"],
                              betas = e5[e5$tau == 0.3,"beta"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5MCP30th <- E5MCP30th[E5MCP30th$Outcome <= quantile(E5MCP30th$Outcome,0.3),]
E5mcp30thpl <- plot_and_save(actual = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.3),]$Outcome,
                                  estimations = E5MCP30th,
                                  mediator_indices = e5[e5$tau == 0.3,"ID"],
                                  filename = "Example5plots/Example5_30thmcpmodel.png",
                                  plotname = "Example Five 30th Percentile MCP HIMA Outcome Predictions")
E5mcp30thpl

E5MCP50th <- estimatedEffects(fit = e5[e5$tau == 0.5,],
                              IDEcoefs = e5[e5$tau == 0.5, "alpha*beta"],
                              betas = e5[e5$tau == 0.5,"beta"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5MCP50th <- E5MCP50th[E5MCP50th$Outcome <= quantile(E5MCP50th$Outcome,0.5),]
E5mcp50thpl <- plot_and_save(actual = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.5),]$Outcome,
                             estimations = E5MCP50th,
                             mediator_indices = e5[e5$tau == 0.5,"ID"],
                             filename = "Example5plots/Example5_50thmcpmodel.png",
                             plotname = "Example Five 50th Percentile MCP HIMA Outcome Predictions")
E5mcp50thpl
E5MCP70th <- estimatedEffects(fit = e5[e5$tau == 0.5,],
                              IDEcoefs = e5[e5$tau == 0.5, "alpha*beta"],
                              betas = e5[e5$tau == 0.5,"beta"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5MCP70th <- E5MCP70th[E5MCP70th$Outcome <= quantile(E5MCP70th$Outcome,0.7),]
E5mcp70thpl <- plot_and_save(actual = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.7),]$Outcome,
                             estimations = E5MCP70th,
                             mediator_indices = e5[e5$tau == 0.7,"ID"],
                             filename = "Example5plots/Example5_70thmcpmodel.png",
                             plotname = "Example Five 70th Percentile MCP HIMA Outcome Predictions")
E5mcp70thpl
## Not run:
# Note: In the following example, M1, M2, and M3 are true mediators.
data(himaDat)
head(himaDat$Example4$PhenoData)
microHIMA.fit <- microHIMA(
  X = himaDat$Example4$PhenoData$Treatment,
  Y = himaDat$Example4$PhenoData$Outcome,
  OTU = himaDat$Example4$Mediator,
  COV = himaDat$Example4$PhenoData[, c("Sex", "Age")],
  FDRcut = 0.05,
  verbose = TRUE
)
microHIMA.fit
## End(Not run)
E4micro <- estimatedEffects(fit = microHIMA.fit,
                            betas = microHIMA.fit[,"beta_hat"],
                            dat = himaDat$Example4,
                            X = himaDat$Example4$PhenoData$Treatment,
                            Y = himaDat$Example4$PhenoData$Outcome,
                            Cov = himaDat$Example4$PhenoData[, c("Sex", "Age")]
                            )
E4micropl <- plot_and_save(actual = himaDat$Example4$PhenoData$Outcome,
                             estimations = E4micro,
                             mediator_indices = microHIMA.fit[,"Index"],
                             filename = "Example4plots/Example4micro.png",
                             plotname = "Example Four Microbiome HIMA Outcome Predictions")
E4micropl
## Not run:
# Note: In the following example, M1, M2, and M3 are true mediators.

head(himaDat$Example5$PhenoData)

qHIMA.fit <- qHIMA(
  X = himaDat$Example5$PhenoData$Treatment,
  M = himaDat$Example5$Mediator,
  Y = himaDat$Example5$PhenoData$Outcome,
  COV = himaDat$Example5$PhenoData[, c("Sex", "Age")],
  tau = c(0.3, 0.5, 0.7),
  scale = FALSE, # Disabled only for simulation data
  Bonfcut = 0.05,
  verbose = TRUE
)
qHIMA.fit
max_rimp <- aggregate(rimp ~ tau + Index, qHIMA.fit, max)
qHIMA.fit <- merge(qHIMA.fit, max_rimp)
qHIMA.fit <- qHIMA.fit[order(qHIMA.fit$tau),c("Index", "alpha_hat", "alpha_se", "beta_hat", "beta_se", "IDE", "rimp", "pmax", "tau")]
## End(Not run)
E5qhima30th <- estimatedEffects(fit = qHIMA.fit[qHIMA.fit$tau == 0.3,],
                              betas = qHIMA.fit[qHIMA.fit$tau == 0.3,"beta_hat"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5qhima30th <- E5qhima30th[E5qhima30th$Outcome <= quantile(E5qhima30th$Outcome,0.3),]
E5q30thpl <- plot_and_save(actual = himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.3),]$Outcome,
                             estimations = E5qhima30th,
                             mediator_indices = qHIMA.fit[qHIMA.fit$tau == 0.3,"Index"],
                             filename = "Example5plots/Example5_30th_quantilemodel.png",
                             plotname = "Example Five 30th Percentile Quantile HIMA Outcome Predictions")
E5q30thpl
E5qhima50th <- estimatedEffects(fit = qHIMA.fit[qHIMA.fit$tau == 0.5,],
                              betas = qHIMA.fit[qHIMA.fit$tau == 0.5,"beta_hat"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5qhima50th <- E5qhima50th[E5qhima50th$Outcome <= quantile(E5qhima50th$Outcome,0.5),]
E5q50thpl <- plot_and_save(actual = (himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.5),]$Outcome),
                           estimations = E5qhima50th,
                           mediator_indices = qHIMA.fit[qHIMA.fit$tau == 0.5,"Index"],
                           filename = "Example5plots/Example5_50th_quantilemodel.png",
                           plotname = "Example Five 50th Percentile Quantile HIMA Outcome Predictions")
E5q50thpl
E5qhima70th <- estimatedEffects(fit = qHIMA.fit[qHIMA.fit$tau == 0.7,],
                              betas = qHIMA.fit[qHIMA.fit$tau == 0.7,"beta_hat"],
                              type = "continuous",
                              dat = himaDat$Example5,
                              X = himaDat$Example5$PhenoData$Treatment,
                              Y = himaDat$Example5$PhenoData$Outcome,
                              Cov = himaDat$Example5$PhenoData[, c("Sex", "Age")])
E5qhima70th <- E5qhima70th[E5qhima70th$Outcome <= quantile(E5qhima70th$Outcome,0.7),]
E5q70thpl <- plot_and_save(actual = (himaDat$Example5$PhenoData[himaDat$Example5$PhenoData$Outcome <= quantile(himaDat$Example5$PhenoData$Outcome, 0.7),]$Outcome),
                           estimations = E5qhima70th,
                           mediator_indices = qHIMA.fit[qHIMA.fit$tau == 0.7,"Index"],
                           filename = "Example5plots/Example5_70th_quantilemodel.png",
                           plotname = "Example Five 70th Percentile Quantile HIMA Outcome Predictions")
E5q70thpl
## Not run:
# Note: In the following example, M1, M2, and M3 are true mediators.
head(himaDat$Example3$PhenoData)
survHIMA.fit <- survHIMA(
  X = himaDat$Example3$PhenoData$Treatment,
  M = himaDat$Example3$Mediator,
  OT = himaDat$Example3$PhenoData$Time,
  status = himaDat$Example3$PhenoData$Status,
  COV = himaDat$Example3$PhenoData[, c("Sex", "Age")],
  scale = FALSE, # Disabled only for simulation data
  FDRcut = 0.05,
  verbose = TRUE
)
survHIMA.fit
E3survhima <- estimatedEffects(fit = survHIMA.fit,
                              betas = survHIMA.fit[,"beta_hat"],
                              type = "continuous",
                              dat = himaDat$Example3,
                              X = himaDat$Example3$PhenoData$Treatment,
                              Y = Surv(himaDat$Example3$PhenoData$Time, himaDat$Example3$PhenoData$Status),
                              Cov = himaDat$Example3$PhenoData[, c("Sex", "Age")])
E3survivalmodelpl <- plot_and_save(actual = Surv(himaDat$Example3$PhenoData$Time, himaDat$Example3$PhenoData$Status),
                                  estimations = E3survhima,
                                  mediator_indices = survHIMA.fit[,"Index"],
                                  filename = "Example3plots/Example3survivalmodel.png",
                                  plotname = "Example Three Survival HIMA Survival Time Predictions")
E3survivalmodelpl
## End(Not run)

write.csv(E1classic, "Example1/E1classic.csv")
write.csv(E1dblasso, "Example1/E1dblasso.csv")
write.csv(E1efficient, "Example1/E1efficient.csv")
write.csv(E1mcp, "Example1/E1mcp.csv")
write.csv(E1mcpEfficient, "Example1/E1mcpEfficient.csv")
write.csv(E2classic, "Example2/E2classic.csv")
write.csv(E2mcp, "Example2/E2mcp.csv")
write.csv(E3dblasso, "Example3/E3dblasso.csv")
write.csv(E3survhima, "Example3/E3survhima.csv")
write.csv(E4dblasso, "Example4/E4dblasso.csv")
write.csv(E4micro, "Example4/E4micro.csv")
write.csv(E5MCP30th, "Example5/E5MCP30th.csv")
write.csv(E5MCP50th, "Example5/E5MCP50th.csv")
write.csv(E5MCP70th, "Example5/E5MCP70th.csv")
write.csv(E5qhima30th, "Example5/E5qhima30th.csv")
write.csv(E5qhima50th, "Example5/E5qhima50th.csv")
write.csv(E5qhima70th, "Example5/E5qhima70th.csv")

