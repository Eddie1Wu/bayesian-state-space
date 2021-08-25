
# Run prior predictive check and fake data check

# Load packages ---------------------------------------------------------------

rm(list=ls()) # Clear workspace

seed <- 1744834965 # Seed used throughout this project
set.seed(seed)

library(HuraultMisc)
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled models
options(mc.cores = parallel::detectCores()) # Parallelise computation
library(ggplot2)
library(cowplot)
source("function_eval.R")
source("function_data.R")


########### Options ###########
mdl_name = "Lme_OrdinalMC"
run_prior <- FALSE
run_fake <- FALSE
n_pt <- 2
n_days <- 7
n_chains <- 4
n_iter <- 2000
###############################


# Define model and parameters -------------------------------------------------

mdl_name = match.arg(mdl_name, c("OrdinalRW", "OrdinalAR", "OrdinalMC",
                                 "Lme_OrdinalAR", "Lme_OrdinalMC",
                                 "Ordinal_OrdinalMC"))
stan_code <- file.path("Model_check", paste0(mdl_name, "_check", ".stan"))
prior_file <- file.path("Result", paste0("prior_", mdl_name, ".rds"))
par0_file <- file.path("Result", paste0("par0_", mdl_name, ".rds"))
fake_file <- file.path("Result", paste0("fake_", mdl_name, ".rds"))

if (any(c(run_prior, run_fake))) {
  compiled_model <- stan_model(stan_code)
}


if (mdl_name == "OrdinalRW") {
  param_pop <- c("sigma_ylat", "ct")
  param_ind <- c()
  param_latent <- c("y_lat")
  param_pred <- c("y_pred")
} else if (mdl_name == "OrdinalAR") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                 "mu_wT", "sigma_wT", "ct")
  param_ind <- c("wY", "wT")
  param_latent <- c("y_lat")
  param_pred <- c("y_pred")
} else if (mdl_name == "OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                 "mu_wT", "sigma_wT", "sigma_lam", "ct")
  param_ind <- c("wY", "wT", "lam")
  param_latent <- c("y_lat")
  param_pred <- c("y_pred")
} else if (mdl_name == "Lme_OrdinalAR") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", 
                 "sigma_wT", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
  param_ind <- c("wY", "wT", "b_s")
  param_latent <- c("y_lat")
  param_pred <- c("y_pred", "sassad_pred")
} else if (mdl_name == "Lme_OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", "sigma_wT", 
                 "sigma_lam", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
  param_ind <- c("wY", "wT", "lam", "b_s")
  param_latent <- c("y_lat")
  param_pred <- c("y_pred", "sassad_pred")
} else if(mdl_name == "Ordinal_OrdinalMC") {
  param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT",
                 "sigma_wT", "sigma_lam", "y_ct", "z_ct")
  param_ind <- c("wY", "wT", "lam")
  param_latent <- c("y_lat")
  param_pred <- c("y_pred", "sassad_pred")
}

param <- c(param_pop, param_ind, param_latent, param_pred)



# Pre-processing --------------------------------------------------------------

df <- expand.grid(Patient = 1:n_pt, Day = 1:n_days) # Generate data
df <- df[order(df$Patient, df$Day),]
df$Bother <- NA
df$Bother[df$Day == 1] <- sample(0:10, n_pt, replace = TRUE)
df$Treatment <- do.call(c,  # Generate treatment by assuming some persistence over time
                    lapply(1:n_pt,
                           function(x) {
                             generate_treatment(c(rbeta(1, 2, 3), rbeta(1, 3, 2)), n_days)
                           }))

# Generate SASSAD for models requiring SASSAD
if (mdl_name %in% c("Lme_OrdinalMC", "Lme_OrdinalAR")) {
  df$SASSAD <- NA
  df$SASSAD[df$Day == 1] <- do.call(c, 
                                    lapply(1:n_pt,
                                           function(i) {
                                             max(round(rnorm(1, mean = df$Bother[df$Patient==i]*2, sd = 0.5) +
                                               rnorm(1, mean = 8, sd = 2)), 3)
                                           }))
} else if (mdl_name == "Ordinal_OrdinalMC") {
  df$SASSAD <- NA
  df$SASSAD[df$Day == 1] <- do.call(c,
                                    lapply(1:n_pt,
                                           function(i) {
                                             min(round(rnorm(1, mean = df$Bother[df$Patient==i]*2, sd = 0.5) +
                                               rnorm(1, mean = 8, sd = 2)), 27)
                                           }))
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

pt <- unique(df$Patient)


# Prior predictive check ------------------------------------------------------

if (run_prior) {
  fit_prior <- sampling(compiled_model,
                        data = data_stan,
                        iter = n_iter,
                        chains = n_chains,
                        pars = param,
                        seed = seed,
                        control = list(adapt_delta = 0.9, max_treedepth = 10)) # Tune hyperparams accordingly here
  saveRDS(fit_prior, file = prior_file)
  par0 <- extract_parameters_check(fit_prior,
                                   param = param,
                                   param_ind = param_ind,
                                   param_latent = param_latent,
                                   param_pred = param_pred,
                                   pt = pt,
                                   data_stan = data_stan)
  saveRDS(par0, file = par0_file)
} else {
  fit_prior <- readRDS(prior_file)
  par0 <- readRDS(par0_file)
}


# Analyse results
if (FALSE) {
  check_hmc_diagnostics(fit_prior)
  param_plot <- param_pop[!(param_pop %in% c("ct"))]  # Indicate which params to plot
  pairs(fit_prior, pars = param_plot)

  # Distribution of parameters
  plot(fit_prior, pars = param_plot, cex.lab = 1.5)
  plot(fit_prior, pars = c(param_plot, paste0(param_ind, "[1]")), plotfun = "hist")
  
  # Posterior predictive distribution
  lapply(pt[1:length(pt)],
         function(i) {
           ggplot(data = subset(par0, Patient == i & Variable == "y_lat"),
                  aes(x = Day, y = Mean, ymin = `5%`, ymax = `95%`)) +
             geom_line() +
             geom_ribbon(alpha = .5) +
             scale_y_continuous(breaks = -2.5:15.5, limits = c(-2.5, 15.5)) +
             theme_bw(base_size = 20) +
             theme(panel.grid.minor.y = element_blank())
         })
}



# Fake data check -------------------------------------------------------------

df_fake <- df # Initialise a data frame for fake data
data_pred <- rstan::extract(fit_prior, pars = param_pred)
bother_pred <- data_pred[[1]]
if (mdl_name %in% c("Ordinal_OrdinalMC", "Lme_OrdinalMC", "Lme_OrdinalAR")) {
  sassad_pred <- data_pred[[2]]
}



# Create fake data by drawing from predicted values
draw <- sample(1:nrow(bother_pred), 1)
df_fake$Bother <- bother_pred[draw, ]
if (mdl_name %in% c("Ordinal_OrdinalMC", "Lme_OrdinalMC", "Lme_OrdinalAR")) {
  draw <- sample(1:nrow(sassad_pred), 1)
  for (i in pt) {
    df_fake$SASSAD[(df$Patient==i) & (df$Day==1)] <- round(sassad_pred[i])
  }
}
df_fake[as.logical(rbinom(nrow(df), 1, 0.1)), "Bother"] <- NA # Generate missing values


# Extract true parameters values
true_param <- par0[, c("Mean", "Variable", "Index")]
colnames(true_param) <- c("Value", "Parameter", "Index")
# Assigne patient ID to patient-dependent params
true_param[["Patient"]] <- NA
id <- (true_param[["Parameter"]] %in% param_ind)
true_param[id, "Patient"] <- pt[true_param[id, "Index"]]

dev.off()
# Look at the data
lapply(pt[1:length(pt)],
       function(patientID) {
         ggplot(data = subset(df_fake, Patient == patientID),
                aes(x = Day, y = Bother)) +
           geom_path() +
           scale_y_continuous(limits = c(0, 10), breaks = 0:10) +
           labs(y = "Bother (fake)") +
           theme_bw(base_size = 15) +
           theme(panel.grid.minor.y = element_blank())
       })


# Fit model with fake data
if (mdl_name == "OrdinalRW") {
  data_fake <- format_stan_data_RW(df_fake)
} else if (mdl_name == "OrdinalAR") {
  data_fake <- format_stan_data_AR(df_fake)
} else if (mdl_name == "OrdinalMC") {
  data_fake <- format_stan_data_MC(df_fake)
} else if (mdl_name == "Lme_OrdinalAR") {
  data_fake <- format_stan_data_LmeAR(df_fake)
} else if (mdl_name == "Lme_OrdinalMC") {
  data_fake <- format_stan_data_LmeMC(df_fake)
} else if (mdl_name == "Ordinal_OrdinalMC") {
  data_fake <- format_stan_data_OrdMC(df_fake)
}


if (run_fake) {
  fit_fake <- sampling(compiled_model,
                       data = data_fake,
                       iter = n_iter,
                       chains = n_chains,
                       pars = param,
                       seed = seed,
                       control = list(adapt_delta = 0.9, max_treedepth = 10))
  saveRDS(fit_fake, file = fake_file)
} else {
  fit_fake <- readRDS(fake_file)
}

# Analyse results
if (FALSE) {
  
  check_hmc_diagnostics(fit_fake)
  param_plot <- param_pop[param_pop != "ct"]
  pairs(fit_fake, pars = param_plot)
  plot(fit_prior, pars = param_plot)
  
  par_fake <- extract_parameters_check(fit_fake,
                                 param = param,
                                 param_ind = param_ind,
                                 param_latent = param_latent,
                                 param_pred = param_pred,
                                 pt = pt,
                                 data_stan = data_fake)
  
  # Recover the population parameters
  tmp <- merge(subset(par_fake, Variable %in% c(param_pop, param_ind)),
               change_colnames(true_param, c("Parameter", "Value"), c("Variable", "True")),
               by = c("Variable", "Patient"))
  tmp <- subset(tmp, !(Variable %in% "ct"))
  tmp$Patient <- factor(tmp$Patient, levels = pt)
  # Plot the population parameters
  ggplot(data = subset(tmp, Variable %in% param_pop),
         aes(x = Variable)) +
    geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
    geom_point(aes(y = True), colour = "#E69F00", size = 2) +
    coord_flip() +
    labs(x = "", y = "Estimate") +
    theme_bw(base_size = 15)
  
  # Recover the patient-dependent parameters
  lapply(intersect(c("wY", "wT", "lam", "b_s"), param_ind), # Change the variables in the vector accordingly
         function(var_name) {
           # Coefficient plot
           tmp <- subset(tmp, Variable == var_name)
           tmp$Patient <- factor(tmp$Patient, levels = tmp[order(tmp$Mean), "Patient"])
           p1 <- ggplot(data = tmp,
                        aes(x = Patient)) +
             geom_pointrange(aes(y = Mean, ymin = `5%`, ymax = `95%`)) +
             geom_point(aes(y = True), colour = "#E69F00") +
             coord_flip() +
             labs(x = "", y = "Estimate") +
             theme_bw(base_size = 15)
           
           # Coverage plot
           p2 <- HuraultMisc::plot_coverage(rstan::extract(fit_fake, pars = var_name)[[1]],
                                            true_param[true_param[["Parameter"]] == var_name, "Value"])
           
           plot_grid(p1, p2, ncol = 2)
         })
  
  # Posterior predictive checks
  ppc <- prepare_ppc(fit_fake, df_fake, par_fake, get_index(pt, data_fake$t_max, 1))
  ppc$Patient <- as.character(ppc$Patient)
  lapply(sample(pt, 1),
         function(pid) {
           plot_ppc(ppc, patientID = pid)
         })
  
}






