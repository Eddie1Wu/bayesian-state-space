

# Load libraries ---------------------------------------------------------------

rm(list = ls())
pacman::p_load(TanakaData, tidyr, dplyr, stargazer, leaps, glmnet,
               randomForest, naivebayes, nnet, MASS, splines2, ggplot2)
source("functions.R")


######## Options ########

seed <- 1744384695
set.seed(seed)
fold <- 10   # for cross validation

#########################


# Pre-process data -------------------------------------------------------------

df_bother <- reshape_SWET_to_wide(SWET)
df_sassad <- reshape_sassad_to_wide(extract_signs_SWET())
df <- full_join(df_sassad, df_bother)  # Merge Bother score with SASSAD
df <- drop_na(df)

# Drop irrelevant variables
df <- df[-c(3:8)]    # Drop area total score
df <- df[-c(10:36)]  # Drop Bother on days 1 - 27
df <- subset(df, (Week != 8) & (Week != 20) & (Week != 0))

# Drop missing observations
df <- df %>%
  drop_na(SASSAD) %>%
  drop_na(Bother.28)

summary(df)



# Run models ------------------------------------------------------------------- 

# Variable selection to choose the most important signs
regfit <- regsubsets(Bother.28 ~ Erythema_Total + Exudation_Total 
                     + Excoriation_Total + Dryness_Total + Cracking_Total 
                     + Lichenification_Total, data = df)
reg.summary = summary(regfit)
names(reg.summary)

par(mfrow = c(2,2))  # Plot the criteria
plot(reg.summary$rss, xlab = "Number of Variables", ylab = "RSS", type = "l", cex.lab=1.3)
plot(reg.summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l", cex.lab=1.3)
points(which.max(reg.summary$adjr2), reg.summary$adjr2[which.max(reg.summary$adjr2)],
       col = "red", cex = 2, pch = 20)
plot(reg.summary$cp, xlab = "Number of Variables", ylab = "Cp", type = "l", cex.lab=1.3)
points(which.min(reg.summary$cp), reg.summary$cp[which.min(reg.summary$cp)], 
       col = "red", cex = 2, pch = 20)
plot(reg.summary$bic, xlab = "Number of Variables", ylab = "BIC", type = "l", cex.lab=1.3)
points(which.min(reg.summary$bic), reg.summary$bic[which.min(reg.summary$bic)], 
       col = "red", cex = 2, pch = 20)



# Variable selection with LASSO
y = data.matrix(df[, colnames(df) == "Bother.28"])
x = data.matrix(dplyr::select(df, ends_with("Total")))
fit <- cv.glmnet(x, y, type.measure="mse", alpha=1, family = "gaussian")
lasso.fit <- glmnet(x, y, family = "gaussian", alpha = 1, lambda = fit$lambda.1se)
coef(lasso.fit)



# Cross validation using regression and classification models to predict Bother

result <- data.frame()
model_list <- c("naiveBayes", "naiveBayes(sum)", "randomForest", 
                "randomForest(sum)", "multi-logistic", "multi-logistic(sum)",
                "ordinal logistic", "ordinal logistic(sum)")

for (model in model_list) {
  
  folds <- make_folds(df, fold)
  out <- cross_validate_2(model, df, folds)
  result <- make_dataframe2(result, out, model)
  print(paste0(model, " finished"))
  
}

if (TRUE) {
  rmse_plot <- ggplot(result, aes(x = Model, y = RMSE)) + 
    geom_boxplot() + 
    coord_flip() + 
    theme_bw(base_size = 20) + 
    theme(panel.grid.minor.x = element_blank())
  
  acc_plot <- ggplot(result, aes(x = Model, y = Accuracy)) + 
    geom_boxplot() + 
    coord_flip() + 
    theme_bw(base_size = 20) + 
    theme(panel.grid.minor.x = element_blank())
  
  egg::ggarrange(rmse_plot, 
                 acc_plot + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank()),
                 nrow = 1)
}












