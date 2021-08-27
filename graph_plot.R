
# Plot graphs for the report

# Load packages ---------------------------------------------------------------

rm(list = ls())

seed <- 1744834695 # seed used throughout the project
set.seed(seed)
pacman::p_load("ggplot2", "cowplot", "reshape2")
color_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Validation metrics for Bother -----------------------------------------------

######## Options ########

metric <- "RPS" # MSE or RPS
horizon <- "1day" # 1day or 1wk

#########################

if (FALSE) {
  par_file <- file.path("Result_eval", paste0("summary_", metric, "_", horizon, ".rds"))
  par <- readRDS(par_file)
  par$iteration <- par$iteration + 1

  # Reshape into long data for easier plotting
  par_long <- melt(par, id.vars = c("iteration"), 
                   variable.name = "model", 
                   value.name = "value")
  
  ggplot(data = par_long, aes(x = iteration, y = value, color = model)) +
    # geom_line(size = 1.5) + 
    geom_smooth(size = 1.5, method = "loess", span = 0.8, se = FALSE) +
    scale_color_manual(values = color_palette) + 
    labs(y = expression("Ranked Probability Score"), x = expression("Week")) +
    scale_x_continuous(name="Week", breaks = 1:10) +
    guides(color = guide_legend(title = "Model")) + 
    theme_bw(base_size = 22)
}



# Validation MSE for SASSAD -----------------------------------------------

if (FALSE) {
  par_file <- file.path("Result_eval-sassad", paste0("summary_MSE", ".rds"))
  par_min_file <- file.path("Result_eval-sassad", paste0("summary_MSE_min", ".rds"))
  par_max_file <- file.path("Result_eval-sassad", paste0("summary_MSE_max", ".rds"))
  
  par <- readRDS(par_file)
  par$iteration <- par$iteration + 1
  par <- par[, !(colnames(par) %in% c("Ordinal_OrdinalMC"))]
  # Reshape into long data for easier plotting
  par_long <- melt(par, id.vars = c("iteration"), 
                   variable.name = "model", 
                   value.name = "value")
  
  par_min <- readRDS(par_min_file)
  par_min$iteration <- par_min$iteration + 1
  par_min <- par_min[, !(colnames(par_min) %in% c("Ordinal_OrdinalMC"))]
  # Reshape into long data for easier plotting
  par_min_long <- melt(par_min, id.vars = c("iteration"),
                   variable.name = "model",
                   value.name = "value")

  par_max <- readRDS(par_max_file)
  par_max$iteration <- par_max$iteration + 1
  par_max <- par_max[, !(colnames(par_max) %in% c("Ordinal_OrdinalMC"))]
  # Reshape into long data for easier plotting
  par_max_long <- melt(par_max, id.vars = c("iteration"),
                   variable.name = "model",
                   value.name = "value")

  df <- merge(par_long, par_min_long, by = c("iteration", "model"))
  df <- merge(df, par_max_long, by = c("iteration", "model"))
  
  # Sort the data frame to get min and max values
  iteration <- df$iteration
  model <- df$model
  df_plot <- data.frame(iteration, model)
  out <- matrix(NA, 8, 3)
  for (i in 1:length(df$value)) {
    x <- sort(df[i, c(3,4,5)])
    out[i,1] <- x[[1]]
    out[i,2] <- x[[2]]
    out[i,3] <- x[[3]]
  }
  df_plot <- cbind(df_plot, out)
  colnames(df_plot) <- c("iteration", "model", "val_min", "val", "val_max")
  df_plot$val_min <- sqrt(df_plot$val_min)
  df_plot$val_max <- sqrt(df_plot$val_max)
  df_plot$val <- sqrt(df_plot$val)
  
  ggplot(data = df_plot, aes(x = iteration, group = model, color = model)) +
    geom_pointrange(aes(y = val, ymin = val_min, ymax = val_max), alpha = 0.5) +
    labs(y = expression("Root Mean Square Error"), x = expression("Iteration")) +
    theme_grey(base_size = 22)
}



# Compare fitted models --------------------------------------------------------

mdl1 <- "Lme_OrdinalAR"
mdl2 <- "Lme_OrdinalMC"
par <- "mu_wY"

file1 <- file.path("Result", paste0("fit_", mdl1, ".rds"))
file2 <- file.path("Result", paste0("fit_", mdl2, ".rds"))
fit1 <- readRDS(file1)
fit2 <- readRDS(file2)

if (FALSE) {

  par1 <- rstan::extract(fit1, pars = par)[[1]]
  par2 <- rstan::extract(fit2, pars = par)[[1]]
  
  df <- data.frame(par1)
  df$Model <- rep("Lme_OrdinalAR",length(df$par1))
  colnames(df) <- c("Value", "Model")
  df2 <- data.frame(par2)
  df2$Model <- rep("Lme_OrdinalMC",length(df2$par2))
  colnames(df2) <- c("Value", "Model")
  df <- rbind(df, df2)
  
  ggplot(df, aes(x=Value, color=Model)) +
    geom_histogram(fill="white", alpha=0.5, position="identity") +
    theme_grey(base_size = 22)
  
}



# Plot latent states
pacman::p_load(TanakaData)
mdl <- "Lme_OrdinalMC"
pt <- 1003
par_file <- file.path("Result", paste0("par_", mdl, ".rds"))

if (FALSE) {
  par <- readRDS(par_file)
  latent_mean <- par$Mean[(par$Variable == "y_lat") & (par$Patient == pt)]
  latent_sd <- par$sd[(par$Variable == "y_lat") & (par$Patient == pt)]
  latent_plus <- latent_mean + latent_sd
  latent_minus <- latent_mean - latent_sd
  
  bother <- SWET$Bother[SWET$Patient == pt]
  factor <- mean(bother, na.rm = TRUE) / mean(latent_mean)
  latent_mean <- latent_mean * factor
  latent_plus <- latent_plus * factor
  latent_minus <- latent_minus * factor
  
  df <- data.frame(latent_mean)
  df$plus <- latent_plus
  df$minus <- latent_minus
  df$Observed <- bother
  df$Day <- 1:length(df$plus)
  colnames(df) <- c("Latent", "Plus", "Minus", "Observed", "Day")
  
  colors <- c("Observed" = "black", "Latent" = "red")
  ggplot(df, aes(x=Day, ymin = Minus, ymax = Plus)) + 
    geom_line(aes(y = Observed, color="Observed")) +
    geom_line(aes(y = Latent, color = "Latent")) + 
    geom_ribbon(alpha=0.2, fill = "red") +
    labs(x = "Day", y = "Magnitude", color = "Legend") +
    scale_color_manual(values = colors) +
    theme_grey(base_size = 22)
}




# Param plot Lme_OrdinalMC -----------------------------------------------------

## Mixed effects model params plot
mdl <- "Lme_OrdinalMC"
par <- c("wS", "mu_bs", "sigma_bs")

file <- file.path("Result", paste0("fit_", mdl, ".rds"))
fit <- readRDS(file)

if (FALSE) {
  
  par <- rstan::extract(fit, pars = par)
  df1 <- data.frame(par[[1]])
  df1$mu_bs <- par[[2]]
  df1$sigma_bs <- par[[3]]
  colnames(df1) <- c("wS", "mu_bs", "sigma_bs")
  
  ggplot(df1, aes(x = wS)) +
    geom_histogram(fill="white", color = "orange", position="identity", bins=35) +
    theme_grey(base_size = 22)
  
  ggplot(df1, aes(x = mu_bs)) +
    geom_histogram(fill="white", color = "orange", position="identity", bins=35) +
    theme_grey(base_size = 22)
  
  ggplot(df1, aes(x = sigma_bs)) +
    geom_histogram(fill="white", color = "orange", position="identity", bins=35) +
    theme_grey(base_size = 22)
}

## Bother ordinal logistic cut points 
# I forgot to save the cut points for the fully estimated Lme_OrdinalMC, as such, 
# to plot the density for Lme_OrdinalMC, one could take Lme_OrdinalAR fit to plot
# those and the output will be very similar. Alternatively, One could take the
# 10th iteration par file from forward chaining for Lme_OrdinalMC 
# and plot that. The output will be the same.

# mdl <- "Lme_OrdinalAR"
# par <- c("ct")
# file <- file.path("Result", paste0("fit_", mdl, ".rds"))
# fit <- readRDS(file)
# 
# if (FALSE) {
#   
#   par <- rstan::extract(fit, pars = par)
#   df <- data.frame(par[[1]])
#   colnames(df) <- c("point1", "point2", "point3", "point4", "point5", "point6",
#                     "point7", "point8", "point9", "point10")
#   df$id <- 1:length(df$point1)
#   
#   df <- melt(df, id.vars = c("id"), 
#              variable.name = "model", 
#              value.name = "value")
#   
#   ggplot(df, aes(x=value, color=model)) +
#     geom_histogram(fill="white", alpha=0.5, position="identity", binwidth = 0.05) +
#     theme_grey(base_size = 22)
# }








