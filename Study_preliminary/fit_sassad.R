

# Load libraries ---------------------------------------------------------------

rm(list = ls())
pacman::p_load(TanakaData, tidyr, dplyr, stargazer, glmnet,
               ggplot2, lme4, splines, splines2)
source("functions.R")


######## Options ########

seed <- 1744834695
set.seed(seed)
fold <- 10   # for cross validation

#########################


# Pre-process data -------------------------------------------------------------

df_bother <- reshape_SWET_to_wide(SWET)
df_sassad <- extract_SASSAD_SWET() 
df_sassad$Patient <- as.factor(df_sassad$Patient)  # convert numeric into factor
df <- full_join(df_sassad, df_bother)  # Merge Bother score with SASSAD

if (TRUE) {
  df$WeekMean.1 <- rowMeans(df[c(4:10)], na.rm = TRUE)  # Generate weekly mean
  df$WeekMean.2 <- rowMeans(df[c(11:17)], na.rm = TRUE)
  df$WeekMean.3 <- rowMeans(df[c(18:24)], na.rm = TRUE)
  df$WeekMean.4 <- rowMeans(df[c(25:31)], na.rm = TRUE)
  
  df$WeekMax.1 <- apply(df[,4:10], 1, max, na.rm = TRUE)  # Generate weekly max
  df$WeekMax.2 <- apply(df[,11:17], 1, max, na.rm = TRUE)
  df$WeekMax.3 <- apply(df[,18:24], 1, max, na.rm = TRUE)
  df$WeekMax.4 <- apply(df[,25:31], 1, max, na.rm = TRUE)
  
  df$WeekMax.1[df$WeekMax.1 == -Inf] <- NaN 
  df$WeekMax.2[df$WeekMax.2 == -Inf] <- NaN 
  df$WeekMax.3[df$WeekMax.3 == -Inf] <- NaN 
  df$WeekMax.4[df$WeekMax.4 == -Inf] <- NaN 
}

df <- drop_na(df)
df$ID <- 1:length(df$Patient)



# Run models ------------------------------------------------------------------- 

# Linear model of SASSAD on 28 days of Bother
tmp <- select(df, SASSAD | starts_with("Bother"))
lm.fit = lm(SASSAD ~ ., data = tmp)

out <- coef(summary(lm.fit)) # Extract coefficients and standard errors
estimate <- out[,1]
stderror <- out[,2]
df <- data.frame(estimate, stderror) # Create data frame for plotting
df$indep_var <- rownames(df)
df <- df[2:length(df$estimate),]
df$max <- df$estimate + 2*df$stderror
df$min <- df$estimate - 2*df$stderror
df$day <- 1:length(df$estimate)
df_highlight <- data.frame(subset(df, df$indep_var == "Bother.28"))

ggplot(data = df, aes(x = day, y = estimate, ymin = min, ymax = max)) +
  geom_pointrange(alpha = 0.5) +
  geom_pointrange(data = df_highlight, 
                  aes(x = day, y = estimate,
                      ymin = min, ymax = max), color = "red") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme_bw(base_size = 20) +
  # theme(panel.grid.minor.x = element_blank(),
  #       axis.text.y = element_blank()) +
  ylab("Coefficient size with 2 std dev.") + 
  xlab("Day in ascending order") +
  scale_x_continuous(breaks = seq(1, 28, by = 1))
# stargazer(lm.fit, type = "html", align = TRUE, out = "28days.htm")



# Best variable selection by LASSO
tmp.y <- data.matrix(tmp[,colnames(tmp) == "SASSAD"])
tmp.x <- data.matrix(tmp[, !(colnames(tmp) == "SASSAD")])
lasso.fit <- cv.glmnet(tmp.x, tmp.y, type.measure = "mse", alpha = 1, family = "gaussian")
out <- coef(lasso.fit)
Day <- 1:(length(out)-1)
df <- data.frame(Day)
for (i in 1:length(df$Day)) {
  df$estimate[i] <- out[i+1]
}


ggplot(data = df, aes(x=Day, y = estimate)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw(base_size = 20) +
  # theme(panel.grid.minor.x = element_blank(),
  #       axis.text.y = element_blank()) +
  ylab("Coefficient size") + 
  xlab("Day in ascending order") +
  scale_x_continuous(breaks = seq(1, 28, by = 1))
  



# Linear model of SASSAD on weekly mean and max of Bother
mean.fit = lm(SASSAD ~ WeekMean.1 + WeekMean.2 + WeekMean.3 + WeekMean.4, data = df)
summary(mean.fit)

max.fit = lm(SASSAD ~ WeekMax.1 + WeekMax.2 + WeekMax.3 + WeekMax.4, data = df)
summary(max.fit)


# Linear mixed effects model of SASSAD
lmer1.fit <- lmer(SASSAD ~ WeekMean.1 + WeekMean.2 + WeekMean.3 + WeekMean.4
                   + (1|Patient), data = df)

lmer2.fit <- lmer(SASSAD ~ WeekMax.1 + WeekMax.2 + WeekMax.3 + WeekMax.4
                  + (1|Patient), data = df)

stargazer(mean.fit, lmer1.fit, max.fit, lmer2.fit, type = "html", 
          align = TRUE, out = "weekly.htm")



# Cross validation to check non-linearity
result <- data.frame()
model_list <- c("mixed intercept", "mixed slope", "I-splines", "natural splines", 
                "polynomial", "linear", "intercept only")

for (model in model_list) {
  if (startsWith(model, "mixed")) {
    by_patient = FALSE
  } else {
    by_patient = TRUE
  }
  
  folds <- make_folds(df, fold, by_patient = by_patient)
  out <- cross_validate_1(model, df, folds, by_patient = by_patient)
  result <- make_dataframe(result, out, model)
}

if (TRUE) {
  rmse_plot <- ggplot(result, aes(x = Model, y = RMSE)) + 
    geom_boxplot() + 
    coord_flip() + 
    theme_bw(base_size = 20) + 
    theme(panel.grid.minor.x = element_blank())
  
  mae_plot <- ggplot(result, aes(x = Model, y = MAE)) + 
    geom_boxplot() + 
    coord_flip() + 
    theme_bw(base_size = 20) + 
    theme(panel.grid.minor.x = element_blank())
  
  r2_plot <- ggplot(result, aes(x = Model, y = R2)) + 
    geom_boxplot() + 
    coord_flip() + 
    theme_bw(base_size = 20) + 
    theme(panel.grid.minor.x = element_blank())
  
  egg::ggarrange(rmse_plot, 
                 mae_plot + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank()),
                 r2_plot + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank()),
                 nrow = 1)
}



# Visualise the relationship between sassad and Bother
ggplot(df, aes(x = Bother.28, y = SASSAD)) +
  geom_point() +
  geom_smooth(method = "gam", aes(colour = "Non-linear")) +
  geom_smooth(method = "lm", aes(colour = "Linear")) + 
  scale_colour_manual(name="legend", values=c("blue", "red")) +
  labs(x = "Bother") +
  theme_bw(base_size = 20)
  
  
  
  



