
# Fit model and obtain parameter estimates

# Load packages ---------------------------------------------------------------

rm(list=ls()) # Clear workspace

seed <- 1744834695 # Seed used throughout this project
set.seed(seed)

library(TanakaData)
library(HuraultMisc)
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled models
options(mc.cores = parallel::detectCores()) # Parallelise computation
# library(shinystan)
source("function_data.R")
source("function_eval.R")


########### Options ###########
mdl_name = "Lme_OrdinalMC"
run <- FALSE
n_chains <- 4
n_iter <- 2000
###############################


# Define model and parameters -------------------------------------------------

mdl_name = match.arg(mdl_name, c("OrdinalRW", "OrdinalAR", "OrdinalMC",
                                 "Lme_OrdinalAR", "Lme_OrdinalMC",
                                 "Ordinal_OrdinalMC"))
stan_code <- file.path("Model", paste0(mdl_name, ".stan"))
suff <- paste0(mdl_name, ".rds")
res_file <- file.path("Result", paste0("fit_", suff))
par_file <- file.path("Result", paste0("par_", suff))


if (mdl_name == "OrdinalRW") {
  param_pop <- c("sigma_ylat", "ct")
  param_ind <- c()
  param_latent <- c("y_lat")
} else if (mdl_name == "OrdinalAR") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                 "mu_wT", "sigma_wT", "ct")
  param_ind <- c("wY", "wT")
  param_latent <- c("y_lat")
} else if (mdl_name == "OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                 "mu_wT", "sigma_wT", "sigma_lam", "ct")
  param_ind <- c("wY", "wT", "lam")
  param_latent <- c("y_lat")
} else if (mdl_name == "Lme_OrdinalAR") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", 
                 "sigma_wT", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
  param_ind <- c("wY", "wT", "b_s")
  param_latent <- c("y_lat")
} else if (mdl_name == "Lme_OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", "sigma_wT", 
                 "sigma_lam", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
  param_ind <- c("wY", "wT", "lam", "b_s")
  param_latent <- c("y_lat")
} else if(mdl_name == "Ordinal_OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT",
                 "sigma_wT", "sigma_lam", "y_ct", "z_ct")
  param_ind <- c("wY", "wT", "lam")
  param_latent <- c("y_lat")
}

param <- c(param_pop, param_ind, param_latent)


# Functions for plotting ------------------------------------------------------

plot_latent_state <- function(pt) {
  # Plot the fitted latent states for a patient
  # Black line is observed bother, red line is latent state
  # 
  # Arg:
  # pt: patient ID
  
  dev.off()
  latent <- par$Mean[(par$Variable == "y_lat") & (par$Patient == pt)]
  bother <- SWET$Bother[SWET$Patient == pt]
  factor <- max(bother, na.rm = TRUE) / max(latent)
  latent <- latent * factor
  plot(latent, type = "l", col = "red")
  lines(bother, type = "l")
  
}


plot_patient_coef <- function(fit, par_names, pt, CI = c(.05, .95)) {
  # Plot patient coefficient estimates from fitted model
  # Patient are ordered according to the first parameter in par_names
  #
  # Arg:
  # fit: stanfit object
  # par_names: vector of names of the patient-dependent parameters to plot
  # pt: vector of patient ID (same order as the patients in the model)
  # CI: (optional) vector of length two indicating the credible interval lower and upper bounds
  #
  # Returns:
  # List of ggplot of patient coefficient estimates
  
  library(ggplot2)
  
  # Extract and summarise posterior
  tmp <- rstan::extract(fit, pars = par_names)
  d <- do.call(rbind,
               lapply(1:length(tmp),
                      function(i) {
                        data.frame(Patient = factor(pt, levels = rev(pt)),
                                   Mean = apply(tmp[[i]], 2, mean),
                                   Lower = apply(tmp[[i]], 2, function(x) {quantile(x, probs = min(CI))}),
                                   Upper = apply(tmp[[i]], 2, function(x) {quantile(x, probs = max(CI))}),
                                   Variable = names(tmp)[i])
                      }))
  
  # Order patients by the mean value of par_names[1]
  par1 <- subset(d, Variable == par_names[1])
  d$Patient <- factor(d$Patient, levels = par1$Patient[order(par1$Mean, decreasing = TRUE)])
  
  # Plot
  lapply(1:length(par_names),
         function(i) {
           ggplot(data = subset(d, Variable == par_names[i]),
                  aes(x = Patient, y = Mean, ymin = Lower, ymax = Upper)) +
             geom_pointrange() +
             coord_flip() +
             theme_bw(base_size = 20) +
             theme(panel.grid.minor.x = element_blank(),
                   axis.text.y = element_blank())
         })
}


# Pre-processing --------------------------------------------------------------

if (mdl_name == "OrdinalRW" | mdl_name == "OrdinalAR" | mdl_name == "OrdinalMC") {
  df <- process1_SWET(SWET)
} else if (mdl_name == "Lme_OrdinalAR" | mdl_name == "Lme_OrdinalMC") {
  df <- process2_SWET(SWET, extract_SASSAD_SWET())
} else if (mdl_name == "Ordinal_OrdinalMC") {
  df <- process3_SWET(SWET, extract_SASSAD_SWET())
}

if (mdl_name == "OrdinalRW") {
  data_stan <- format_stan_data_RW(df)
} else if (mdl_name == "OrdinalAR") {
  data_stan <- format_stan_data_AR(df)
} else if (mdl_name == "OrdinalMC") {
  data_stan <- format_stan_data_MC(df)
} else if (mdl_name == "Lme_OrdinalAR") {
  data_stan <- format_stan_data_LmeAR(df)
} else if (mdl_name == "Lme_OrdinalMC") {
  data_stan <- format_stan_data_LmeMC(df)
} else if (mdl_name == "Ordinal_OrdinalMC") {
  data_stan <- format_stan_data_OrdMC(df)
}

pt <- unique(df[["Patient"]])


# Fit model -------------------------------------------------------------------

if (run) {
  
  # Fit model
  fit <- stan(file = stan_code,
              data = data_stan,
              iter = n_iter,
              chains = n_chains,
              pars = param,
              seed = seed,
              control = list(adapt_delta = 0.9, max_treedepth = 10)) # Change hyperparams here
  saveRDS(fit, file = res_file) # Save model fit
  
  # Extract parameters
  par <- extract_parameters(fit,
                            param = param,
                            param_ind = param_ind,
                            param_latent = param_latent,
                            pt = pt,
                            data_stan = data_stan)
  saveRDS(par, file = par_file) # Save parameters
} else {
  fit <- readRDS(res_file) # Load fitted model
  par <- readRDS(par_file) # Load model param
}


# Results ---------------------------------------------------------------------

if (FALSE) {
  
  #Diagnostics
  check_hmc_diagnostics(fit)
  get_num_divergent(fit)

  param_plot <- c("sigma_ylat", "b_ylat", "mu_wY", "mu_wT", "sigma_lam") # Choose which param to plot
  pairs(fit, pars = param_plot) # Pairs plot
  plot(fit, pars = param_plot, plotfun = "trace") # Trace plot
  
  plot_latent_state(pt = 1001)
  
  pl <- plot_patient_coef(fit, c("wY", "wT", "lam"), pt)
  cowplot::plot_grid(pl[[1]] 
                     + labs(y = expression(paste("Persistence (", w[Y]^(k), ")", sep = ""))) 
                     + coord_flip(ylim = c(0, 1)),
                     pl[[2]] 
                     + labs(y = expression(paste("Treatment response (", w[T]^(k), ")", sep = ""))) 
                     + coord_flip(ylim = c(-2.5, 1.5))
                     + geom_hline(yintercept = 0, color = "red"),
                     pl[[3]] 
                     + labs(y = expression(paste("Flare trigger (", Lam^(k), ")", sep = ""))) 
                     + coord_flip(ylim = c(1, 8)),
                     nrow = 1, labels = "auto")
}


