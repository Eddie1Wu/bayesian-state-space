
# Validate model with forward chaining to predict SASSAD

# Load packages ---------------------------------------------------------------

rm(list=ls()) # Clear workspace

seed <- 1744834695 # Seed used throughout this project
set.seed(seed)

library(TanakaData)
library(HuraultMisc)
library(rstan)
rstan_options(auto_write = TRUE) # Save compiled models
options(mc.cores = parallel::detectCores()) # Parallelise computation
library(foreach)
library(doParallel)
source("function_data.R")
source("function_eval.R")


########### Options ###########
mdl_name = "Ordinal_OrdinalMC"
run <- TRUE
n_chains <- 4
n_iter <- 2000
n_cluster <- 2  # floor(parallel::detectCores() / n_chains)
t_horizon <- 1
max_training <- 3
###############################



# Define model and parameters -------------------------------------------------

mdl_name = match.arg(mdl_name, c("Lme_OrdinalAR", "Lme_OrdinalMC",
                                 "Ordinal_OrdinalMC"))
is_stan_model <- TRUE
suff <- paste0(mdl_name, "_eval.rds")


if (is_stan_model) {
  
  stan_code <- file.path("Model_eval-sassad", paste0(mdl_name, "_eval.stan"))
  
  if (mdl_name == "Lme_OrdinalAR") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", 
                   "sigma_wT", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
    param_ind <- c("wY", "wT", "b_s")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "sassad_pred")
  } else if (mdl_name == "Lme_OrdinalMC") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", "sigma_wT", 
                   "sigma_lam", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
    param_ind <- c("wY", "wT", "lam", "b_s")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "sassad_pred")
  } else if(mdl_name == "Ordinal_OrdinalMC") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT",
                   "sigma_wT", "sigma_lam", "y_ct", "z_ct")
    param_ind <- c("wY", "wT", "lam")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "sassad_pred")
  }
  
  param <- c(param_pop, param_ind, param_latent, param_pred)
  
  compiled_model <- stan_model(stan_code)
}



# Functions for forward-chaining ----------------------------------------------

if (is_stan_model) {
  
  format_stan_test <- function(df) {
    with(df, 
         list(N_test = length(Bother),
              test_t_max = aggregate(Day ~ Patient, FUN = length)$Day,
              y_test = Bother))
  }
  
  
  get_stan_data <- function(mdl, df_train, df_test) {
    if (mdl == "Lme_OrdinalAR") {
      train <- format_stan_data_LmeAR(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    } else if (mdl == "Lme_OrdinalMC") {
      train <- format_stan_data_LmeMC(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    } else if (mdl == "Ordinal_OrdinalMC") {
      train <- format_stan_data_OrdMC(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    }
    
    return(data_stan)
  }
}



# Pre-processing --------------------------------------------------------------

preprocess_data <- function(mdl_name) {
  
  if (mdl_name == "Lme_OrdinalAR" | mdl_name == "Lme_OrdinalMC") {
    df <- process2_SWET(SWET, extract_SASSAD_SWET())
  } else if (mdl_name == "Ordinal_OrdinalMC") {
    df <- process3_SWET(SWET, extract_SASSAD_SWET())
  }
  
  df <- fill_a_patient(df)
  df <- drop_missing_SASSADs(df)
  
  return(df)
}

df <- preprocess_data(mdl_name)



# Run forward chaining --------------------------------------------------------

if (run) {
  
  duration <- Sys.time()
  cl <- makeCluster(n_cluster)
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  
  out <- foreach(iter = max_training:0) %dopar% {
    
    library(rstan) # Reload libraries
    rstan_options(auto_write = TRUE) # Save compiled model
    options(mc.cores = parallel::detectCores()) # Parallel computing
    source("function_data.R")
    source("function_eval.R")
    
    sink("log.txt", append = TRUE)
    cat(paste("Starting model", iter, "\n"))
    
    # Split dataset
    if (iter == 0) { # Generate prior data for 1 day so that no real data is used here
      n_pt <- length(unique(df$Patient))
      df_train <- expand.grid(Patient = 1:n_pt, Day = 1:2) # Generate data
      df_train <- df_train[order(df_train$Patient, df_train$Day),]
      df_train$Bother <- NA
      df_train$Bother[df_train$Day == 1] <- sample(0:10, n_pt, replace = TRUE)
      df_train$Treatment <- do.call(c,  # Generate treatment by assuming some persistence over time
                                    lapply(1:n_pt,
                                           function(x) {
                                             generate_treatment(c(rbeta(1, 2, 3), rbeta(1, 3, 2)), 2)
                                           }))
      if (mdl_name %in% c("Lme_OrdinalMC", "Lme_OrdinalAR")) {
        df_train$SASSAD <- NA
        df_train$SASSAD[df_train$Day == 1] <- do.call(c, 
                                                      lapply(1:n_pt,
                                                             function(i) {
                                                               max(ceiling(rnorm(1, mean = df_train$Bother[df_train$Patient==i]*3.5, sd = 1.5) +
                                                                           rnorm(1, mean = 10, sd = 3)), 0)
                                                             }))
      } else if (mdl_name == "Ordinal_OrdinalMC") {
        df_train$SASSAD <- NA
        df_train$SASSAD[df_train$Day == 1] <- do.call(c,
                                                      lapply(1:n_pt,
                                                             function(i) {
                                                               min(ceiling(rnorm(1, mean = df_train$Bother[df_train$Patient==i]*3.5, sd = 1.5) +
                                                                           rnorm(1, mean = 10, sd = 3)), 27)
                                                             }))
      }
      test_days <- 1
      
    } else if (iter == 1) {
      iter <- 4
      train_days <- 1:(iter*7 - 1)
      df_train <- df[df$Day %in% train_days, ]
      test_days <- iter*7
    } else if (iter == 2) {
      iter <- 12
      train_days <- 1:(iter*7 - 1)
      df_train <- df[df$Day %in% train_days, ]
      test_days <- iter*7
    } else if (iter == 3) {
      iter <- 16
      train_days <- 1:(iter*7 - 1)
      df_train <- df[df$Day %in% train_days, ]
      test_days <- iter*7
    }
    
    df_test <- df[df$Day %in% test_days,]
    df_test$Bother[is.na(df_test$Bother)] <- 99
    
    pt <- unique(df_train[["Patient"]])
    par_file <- file.path("Result_eval-sassad", paste0("par", iter, "_", suff))

    if (is_stan_model) {
      
      data_stan <- get_stan_data(mdl_name, df_train, df_test)
      
      fit <- sampling(compiled_model,
                      data = data_stan,
                      iter = n_iter,
                      chains = n_chains,
                      pars = param,
                      control = list(adapt_delta = 0.9, max_treedepth = 10))
      
      par <- extract_parameters_eval(fit,
                                     param = param,
                                     param_ind = param_ind,
                                     param_latent = param_latent,
                                     param_pred = param_pred,
                                     pt = pt,
                                     data_stan = data_stan)
      saveRDS(par, file = par_file) # Save parameters
      
    }
    
    cat(paste("Ending model", iter, "\n"))
    NULL # Return
  }
  
  stopCluster(cl)
  (duration = Sys.time() - duration)
}



# Analyse results -------------------------------------------------------------

########### Options ###########

metric = "MSE"
iterations = 3
mdl_list = c("Lme_OrdinalAR", "Lme_OrdinalMC", "Ordinal_OrdinalMC")

###############################

if (metric == "MSE") {
  
  # Initialise output data frame
  output <- data.frame(matrix(vector(), (iterations+1), length(mdl_list), 
                              dimnames=list(c(), mdl_list)),
                       stringsAsFactors=F)
  output$iteration <- 0:iterations
  output_file <- file.path("Result_eval-sassad", paste0("summary_MSE_max.rds"))
  
  for (model in mdl_list) {
    # Load the data for this model
    df <- preprocess_data(model)
    i = 0 # Keep count
    for (t in 0:iterations) {
      i = i + 1
      # Load test data for this iteration
      if (t == 0) {
        test_days <- 1
      } else if (t == 1) {
        t <- 4
        test_days <- t * 7
      } else if (t == 2) {
        t <- 12
        test_days <- t * 7
      } else if (t == 3) {
        t <- 16
        test_days <- t * 7
      }
      df_test <- df[df$Day %in% test_days,]
      pt <- unique(df_test$Patient)
      
      suff <- paste0(model, "_eval.rds")   # Get file path of the parameters file
      par_file <- file.path("Result_eval-sassad", paste0("par", t, "_", suff))
      par <- readRDS(par_file) # Load the parameters file
      
      # Calculate MSE
      pred <- par$`75%`[(par$Variable == "sassad_pred")]
      actual <- df_test$SASSAD
      both <- data.frame(actual, pred)
      MSE <- mean((both$actual - both$pred)^2)
      print(MSE)
      
      # Append to the output data frame
      output[(i), model] = MSE
      print(paste0("Finished ", t, " out of ", iterations, " for ", model))
    }
  }
  saveRDS(output, output_file)
}






