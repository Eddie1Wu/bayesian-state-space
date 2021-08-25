

# This file contains the other codes (just trying, not cleaned) that
# I have written during the preliminary study.
# This file is just keeping a record of the miscellaneous findings
# and records some of the other explorations of the SWET dataset that I did
# but outputs are not directly relevant hence excluded from final report


# Mixed effect models
mixed.lm1 <- lmer(SASSAD ~ WeekMean.1 + WeekMean.2 + WeekMean.3 + WeekMean.4
                  + (1|Patient), data = sb)
summary(mixed.lm1)
plot(mixed.lm1)
qqnorm(resid(mixed.lm1)) # Plotting check
qqline(resid(mixed.lm1))

mixed.lm2 <- lmer(SASSAD ~ OverallMean + (1|Patient), data = sb)
summary(mixed.lm2)
plot(mixed.lm2)
qqnorm(resid(mixed.lm2)) # Plotting check
qqline(resid(mixed.lm2))

mixed.lm2 <- lmer(SASSAD ~ OverallMean + (1 + OverallMean|Patient), data = sb)
summary(mixed.lm2)
plot(mixed.lm2)
qqnorm(resid(mixed.lm2)) # Plotting check
qqline(resid(mixed.lm2))





# function for k-fold cross validation
kfoldcv <- function(model_type, data, seed, fold) {
  
  set.seed(seed)
  patient_ID <- unique(data$Patient)
  folds <- createFolds(patient_ID, k = fold, list = TRUE, returnTrain = TRUE)
  
  RMSE <- rep(NA, fold)
  MAE <- rep(NA, fold)
  Rsquare <- rep(NA, fold)
  
  for (i in 1:fold) {
    train_id <- patient_ID[folds[[i]]]
    test_id <- patient_ID[!(patient_ID %in% train_id)]
    train <- subset(data, Patient%in%train_id)
    test <- subset(data, Patient%in%test_id)
    
    if (model_type == "intercept") {
      model <- lm(Median ~ 1, data = train)
    } else if (model_type == "linear sassad") {
      model <- lm(Median ~ SASSAD, data = train)
    } else if (model_type == "linear area") {
      model <- lm(Median ~ Head_Neck_Total + Hands_Total + Arms_Total 
                  + Trunk_Total + Feet_Total + Legs_Total, data = train)
    } else if (model_type == "linear sign") {
      model <- lm(Median ~ Erythema_Total + Exudation_Total + Excoriation_Total
                  + Dryness_Total + Cracking_Total + Lichenification_Total, data = train)
    } else if (model_type == "spline sassad") {
      model <- lm(Median ~ iSpline(SASSAD, knots = c(20), degree = 2), data = sb.dropna)
    } else if (model_type == "spline area") {
      model <- lm(Median ~ iSpline(Head_Neck_Total, knots = c(5), degree = 2) +
                    iSpline(Hands_Total, knots = c(5), degree = 2) +
                    iSpline(Arms_Total, knots = c(5), degree = 2) +
                    iSpline(Trunk_Total, knots = c(5), degree = 2) +
                    iSpline(Feet_Total, knots = c(5), degree = 2) +
                    iSpline(Legs_Total, knots = c(5), degree = 2), data = sb.dropna)
    } else if (model_type == "spline sign") {
      model <- lm(Median ~ iSpline(Erythema_Total, knots = c(6), degree = 2) +
                    iSpline(Exudation_Total, knots = c(2), degree = 2) +
                    iSpline(Excoriation_Total, knots = c(5), degree = 2) +
                    iSpline(Dryness_Total, knots = c(7), degree = 2) +
                    iSpline(Cracking_Total, knots = c(2), degree = 2) +
                    iSpline(Lichenification_Total, knots = c(7), degree = 2), data = sb.dropna)
    } else if (model_type == "lmer sassad") {
      model <- lmer(Median ~ SASSAD + (1 | Patient), data = train)
    } else if (model_type == "lmer area") {
      model <- lmer(Median ~ Head_Neck_Total + Hands_Total + Arms_Total 
                    + Trunk_Total + Feet_Total + Legs_Total + (1 | Patient), data = train)
    } else if (model_type == "lmer sign") {
      model <- lmer(Median ~ Erythema_Total + Exudation_Total + Excoriation_Total
                    + Dryness_Total + Cracking_Total + Lichenification_Total 
                    + (1 | Patient), data = train)
    }
    
    prediction <- predict(model, newdata = test, allow.new.levels = TRUE)
    error <- RMSE(prediction, test$Median)
    abs_error <- MAE(prediction, test$Median)
    r2 <- R2(prediction, test$Median)
    
    RMSE[i] <- error
    MAE[i] <- abs_error
    Rsquare[i] <- r2
  }
  
  return(c(RMSE, MAE, Rsquare))
}

# function for creating dataframe for plotting
make_dataframe <- function(out, name) {
  df2 <- data.frame(out[1:10])
  names(df2) <- c("RMSE")
  df2$MAE <- out[11:20]
  df2$R2 <- out[21:30]
  df2$Model <- name
  
  df <- bind_rows(df, df2)
  
  return(df)
}

# 10-fold cv of linear models
fold = 10
seed = 5
df <- data.frame()

out <- kfoldcv(model_type = "intercept", sb.dropna, seed, fold)
df <- make_dataframe(out, "Intercept Only")

out <- kfoldcv(model = "linear sassad", sb.dropna, seed, fold)
df <- make_dataframe(out, "lm(SASSAD)")

out <- kfoldcv(model = "linear area", sb.dropna, seed, fold)
df <- make_dataframe(out, "lm(6 areas)")

out <- kfoldcv(model = "linear sign", sb.dropna, seed, fold)
df <- make_dataframe(out, "lm(6 signs)")

# 10-fold cv of i-splines
out <- kfoldcv(model = "spline sassad", sb.dropna, seed, fold)
df <- make_dataframe(out, "iSpline(SASSAD)")

out <- kfoldcv(model = "spline area", sb.dropna, seed, fold)
df <- make_dataframe(out, "iSpline(6 areas)")

out <- kfoldcv(model = "spline sign", sb.dropna, seed, fold)
df <- make_dataframe(out, "iSpline(6 signs)")

# 10-fold cv of mixed effects models
out <- kfoldcv(model = "lmer sassad", sb.dropna, seed, fold)
df <- make_dataframe(out, "lmer(SASSAD)")

out <- kfoldcv(model = "lmer area", sb.dropna, seed, fold)
df <- make_dataframe(out, "lmer(6 areas)")

out <- kfoldcv(model = "lmer sign", sb.dropna, seed, fold)
df <- make_dataframe(out, "lmer(6 signs)")





# Mixed effect models with mixed levels (patients) in both training and test set
set.seed(seed)
sb.dropna2 <- sb.dropna %>%
  mutate(ID = 1 + 1*(0:(dim(sb.dropna)[1]-1)))

folds <- createFolds(sb.dropna2$ID, k=fold, list = TRUE, returnTrain = TRUE)

RMSE <- rep(NA, fold)
MAE <- rep(NA, fold)
Rsquare <- rep(NA, fold)
for (i in 1:fold) {
  train_id <- folds[[i]]
  test_id <- sb.dropna2$ID[!(sb.dropna2$ID %in% folds[[i]])]
  train <- subset(sb.dropna2, ID%in%train_id)
  test <- subset(sb.dropna2, ID%in%test_id)
  
  lmm <- lmer(Median ~ SASSAD + (1|Patient), data = train)
  prediction <- predict(lmm, newdata = test, allow.new.levels = TRUE)
  error <- RMSE(prediction, test$Median)
  abs_error <- MAE(prediction, test$Median)
  r2 <- R2(prediction, test$Median)
  
  RMSE[i] <- error
  MAE[i] <- abs_error
  Rsquare[i] <- r2
}

out <- c(RMSE, MAE, Rsquare)
df <- make_dataframe(out, "lmer(SASSAD)**")

lmm <- lm(Median ~ SASSAD, data = sb.dropna)
summary(lmm)

# plot the results
ggplot(df, aes(x = Model, y = RMSE)) +
  geom_boxplot()

ggplot(df, aes(x = Model, y = MAE)) +
  geom_boxplot()

ggplot(df, aes(x = Model, y = R2)) +
  geom_boxplot()





# compare ordinal logit regression models with mixed effects
model <- clmm(Median_factor ~ SASSAD,   # 1 regressor only
              data = sb.dropna, Hess = TRUE)
summary(model)

# 4 signs
model <- clmm(Median_factor ~ Erythema_Total + Exudation_Total 
              + Dryness_Total+ Lichenification_Total + (1|Patient), 
              data = sb.dropna, Hess = TRUE)
summary(model)

# 6 signs
model <- clmm(Median_factor ~ Erythema_Total + Exudation_Total + Excoriation_Total 
              + Dryness_Total + Cracking_Total + Lichenification_Total + (1|Patient), 
              data = sb.dropna, Hess = TRUE)
warnings()
summary(model)





# Round the .5 median scores to the next largest integer
sb.dropna$Median <- ceiling(sb.dropna$Median)
sb.dropna$Patient_factor <- as.factor(sb.dropna$Patient)
sb.dropna$Median_factor <- as.factor(sb.dropna$Median)

# function for creating dataframe for plotting
make_dataframe <- function(out, name) {
  df2 <- data.frame(out[1:10])
  names(df2) <- c("RMSE")
  df2$MAE <- out[11:20]
  df2$Model <- name
  
  df <- bind_rows(df, df2)
  
  return(df)
}

# calculate RMSE from a table
calc_table_rmse <- function(table, table_length) {
  total <- 0
  for (i in 1:table_length) {
    for (j in 1:table_length) {
      error <- table[i,j]*((i-1)-(j-1))^2
      total <- total + error
    }
  }
  
  mse <- total/(table_length^2)
  rmse <- sqrt(mse)
  
  return(rmse)
}

# calculate MAE from a table
calc_table_mae <- function(table, table_length) {
  total <- 0
  for (i in 1:table_length) {
    for (j in 1:table_length) {
      error <- table[i,j]*abs((i-1)-(j-1))
      total <- total + error
    }
  }
  
  mae <- total/(table_length^2)
  
  return(mae)
}

# function for k-fold cross validation
kfoldcv <- function(model_type, type, data, seed, fold) {
  
  set.seed(seed)
  patient_ID <- unique(data$Patient)
  folds <- createFolds(patient_ID, k = fold, list = TRUE, returnTrain = TRUE)
  
  RMSE <- rep(NA, fold)
  MAE <- rep(NA, fold)
  
  for (i in 1:fold) {
    train_id <- patient_ID[folds[[i]]]
    test_id <- patient_ID[!(patient_ID %in% train_id)]
    train <- subset(data, Patient%in%train_id)
    test <- subset(data, Patient%in%test_id)
    
    if (model_type == "intercept") {
      model <- lm(Median ~ 1, data = train)
    } else if (model_type == "linear") {
      model <- lm(Median ~ Erythema_Total + Exudation_Total + Excoriation_Total
                  + Dryness_Total + Cracking_Total + Lichenification_Total, data = train)
    } else if (model_type == "ordinal") {
      model <- polr(Median_factor ~ Erythema_Total + Exudation_Total + Excoriation_Total
                    + Dryness_Total + Cracking_Total + Lichenification_Total, data = train, Hess = TRUE)
    } else if (model_type == "linear lasso") {
      x = data.matrix(train[c(46:51)])
      y = data.matrix(train[,59])
      new_x = data.matrix(test[c(46:51)])
      new_y = data.matrix(test[,59])
      fit <- cv.glmnet(x, y, type.measure="mse", alpha=1, family = "gaussian")
      prediction <- predict(fit, s=fit$lambda.1se, newx = new_x)
      prediction <- round(prediction)   # Round to the nearest whole number for comparison with cls
      RMSE[i] <- RMSE(prediction, new_y)
      MAE[i] <- MAE(prediction, new_y)
    } else if (model_type == "ordinal lasso") {
      model <- polr(Median_factor ~ Erythema_Total + Exudation_Total
                    + Dryness_Total + Lichenification_Total, data = train, Hess = TRUE)
    }
    
    if (type == "reg") {
      prediction <- predict(model, newdata = test, allow.new.levels = TRUE)
      prediction <- round(prediction)   # Round to the nearest whole number for comparison with cls
      RMSE[i] <- RMSE(prediction, test$Median)
      MAE[i] <- MAE(prediction, test$Median)
    } else if (type == "cls") {
      prediction <- predict(model, test)
      tab <- table(test$Median_factor, prediction)
      RMSE[i] <- calc_table_rmse(tab, 11)
      MAE[i] <- calc_table_mae(tab, 11)
    }
  }
  
  return(c(RMSE, MAE))
}

# Run the models
fold = 10
seed = 222
df <- data.frame()

out <- kfoldcv(model_type = "intercept", "reg", sb.dropna, seed, fold)
df <- make_dataframe(out, "Intercept")

out <- kfoldcv(model_type = "linear", "reg", sb.dropna, seed, fold)
df <- make_dataframe(out, "Linear")

out <- kfoldcv(model_type = "ordinal", "cls", sb.dropna, seed, fold)
df <- make_dataframe(out, "Ordinal")

out <- kfoldcv(model_type = "linear lasso", "-", sb.dropna, seed, fold)
df <- make_dataframe(out, "Linear lasso")

out <- kfoldcv(model_type = "ordinal lasso", "cls", sb.dropna, seed, fold)
df <- make_dataframe(out, "Ordinal lasso")

# Plot the results
ggplot(df, aes(x = Model, y = RMSE)) +
  geom_boxplot()

ggplot(df, aes(x = Model, y = MAE)) +
  geom_boxplot()

# Variable selection with LASSO
x = data.matrix(sb.dropna[c(46:51)])
y = data.matrix(sb.dropna[,59])
fit <- cv.glmnet(x, y, type.measure="mse", alpha=1, family = "gaussian")
fit2 <- glmnet(x, y, family = "gaussian", alpha = 1, lambda = fit$lambda.1se)
coef(fit2)








