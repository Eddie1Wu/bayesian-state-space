
# Validate model with forward chaining and analyse results at the end

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
mdl_name = "OrdinalMC"
run <- FALSE
n_chains <- 4
n_iter <- 2000
n_cluster <- 2  # floor(parallel::detectCores() / n_chains)
t_horizon <- 7
max_week <- 10
###############################



# Define model and parameters -------------------------------------------------

mdl_name = match.arg(mdl_name, c("OrdinalRW", "OrdinalAR", "OrdinalMC",
                                 "Lme_OrdinalAR", "Lme_OrdinalMC",
                                 "Ordinal_OrdinalMC", "Historical"))

is_stan_model <- !(mdl_name %in% c("Historical"))
suff <- paste0(mdl_name, "_eval.rds")


if (is_stan_model) {
  
  stan_code <- file.path("Model_eval", paste0(mdl_name, "_eval.stan"))
  
  if (mdl_name == "OrdinalRW") {
    param_pop <- c("sigma_ylat", "ct")
    param_ind <- c()
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
  } else if (mdl_name == "OrdinalAR") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                   "mu_wT", "sigma_wT", "ct")
    param_ind <- c("wY", "wT")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
  } else if (mdl_name == "OrdinalMC") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", 
                   "mu_wT", "sigma_wT", "sigma_lam", "ct")
    param_ind <- c("wY", "wT", "lam")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
  } else if (mdl_name == "Lme_OrdinalAR") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", 
                   "sigma_wT", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
    param_ind <- c("wY", "wT", "b_s")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
  } else if (mdl_name == "Lme_OrdinalMC") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT", "sigma_wT", 
                   "sigma_lam", "wS", "mu_bs", "sigma_bs", "sigma_S", "ct")
    param_ind <- c("wY", "wT", "lam", "b_s")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
  } else if(mdl_name == "Ordinal_OrdinalMC") {
    param_pop <- c("sigma_ylat", "b_ylat", "mu_wY", "sigma_wY", "mu_wT",
                   "sigma_wT", "sigma_lam", "y_ct", "z_ct")
    param_ind <- c("wY", "wT", "lam")
    param_latent <- c("y_lat")
    param_pred <- c("y_lat_pred", "lpd", "y_pred", "cdf")
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
    if (mdl == "OrdinalRW") {
      train <- format_stan_data_RW(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    } else if (mdl == "OrdinalAR") {
      train <- format_stan_data_AR(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    } else if (mdl == "OrdinalMC") {
      train <- format_stan_data_MC(df_train)
      test <- format_stan_test(df_test)
      data_stan <- c(train, test)
    } else if (mdl == "Lme_OrdinalAR") {
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
  
  if (mdl_name == "OrdinalRW" | mdl_name == "OrdinalAR" | 
      mdl_name == "OrdinalMC" | mdl_name == "Historical") {
    df <- process1_SWET(SWET)
    df <- fill_a_patient(df)
  } else if (mdl_name == "Lme_OrdinalAR" | mdl_name == "Lme_OrdinalMC") {
    df <- process2_SWET(SWET, extract_SASSAD_SWET())
    df <- fill_a_patient(df)
    df <- drop_missing_SASSAD(df)
  } else if (mdl_name == "Ordinal_OrdinalMC") {
    df <- process3_SWET(SWET, extract_SASSAD_SWET())
    df <- fill_a_patient(df)
  }
  
  return(df)
}

df <- preprocess_data(mdl_name)



# Run forward chaining --------------------------------------------------------

if (run) {
  
  duration <- Sys.time()
  cl <- makeCluster(n_cluster)
  registerDoParallel(cl)
  
  writeLines(c(""), "log.txt")
  
  out <- foreach(iter = max_week:0) %dopar% {
    
    library(rstan) # Reload libraries
    rstan_options(auto_write = TRUE) # Save compiled model
    options(mc.cores = parallel::detectCores()) # Parallel computing
    source("function_data.R")
    source("function_eval.R")
    
    sink("log.txt", append = TRUE)
    cat(paste("Starting model", iter, "\n"))

    # Split dataset
    if (iter == 0) {
      train_days <- 1:2
      test_days <- 3:9
    } else {
      train_days <- 1:(iter*t_horizon)
      test_days <- (iter*t_horizon + 1):((iter+1)*t_horizon)
    }
    df_train <- df[df$Day %in% train_days, ]
    df_test <- df[df$Day %in% test_days,]
    df_test$Bother[is.na(df_test$Bother)] <- 99
    
    pt <- unique(df_train[["Patient"]])
    par_file <- file.path("Result_eval", paste0("par", iter, "_", suff))
    
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
      
    } else if (mdl_name == "Historical") {
      out <- calculate_lpd_and_ypred(df_train, df_test, pt)
      cdf <- calculate_cdf(df_train, df_test, pt)
      out <- dplyr::bind_rows(out, cdf)
      saveRDS(out, file = par_file) # Save parameters
    }
    
    cat(paste("Ending model", iter, "\n"))
    NULL # Return
  }
  
  stopCluster(cl)
  (duration = Sys.time() - duration)
}



# Analyse results -------------------------------------------------------------

########### Options ###########

metric = "RPS"  # Choose which metric: RPS, MSE, LPD
t_pred = 7  # Number of days to predict in [1,7]
iterations = 9 # Number of weeks trained
mdl_list = c("Historical", "OrdinalRW", "OrdinalAR", "OrdinalMC",
             "Lme_OrdinalAR", "Lme_OrdinalMC", "Ordinal_OrdinalMC")

###############################

if (metric == "RPS") {
  
  # Initialise output data frame
  output <- data.frame(matrix(vector(), (iterations+1), length(mdl_list), 
                           dimnames=list(c(), mdl_list)),
                    stringsAsFactors=F)
  output$iteration <- 0:iterations
  output_file <- file.path("Result_eval", paste0("summary_RPS.rds"))
  
  for (model in mdl_list) {
    # Load the data for this model
    df <- preprocess_data(model)
    
    for (t in 0:iterations) {
      # Load test data for this iteration
      if (t == 0) {
        test_days <- 3:(3+(t_pred-1))
      } else {
        test_days <- (t*t_horizon + 1):(t*t_horizon + t_pred)
      }
      df_test <- df[df$Day %in% test_days,]
      pt <- unique(df_test$Patient)
      
      suff <- paste0(model, "_eval.rds")  # Get file path of the parameters file
      par_file <- file.path("Result_eval", paste0("par", t, "_", suff))
      par <- readRDS(par_file) # Load the parameters file
      
      # Calculate actual CDF from test set
      test_t_max <- aggregate(df_test$Day ~ df_test$Patient, FUN = length)$`df_test$Day`
      test_start <- rep(NA, length(pt))
      test_end <- rep(NA, length(pt))
      for (k in 1:length(pt)) {
        if (k == 1) {
          test_start[k] = 1
        } else {
          test_start[k] = test_end[k-1] + 1
        }
        test_end[k] = test_start[k] - 1 + test_t_max[k]
      }
      out <- get_index(pt, test_t_max, 11) # Create an empty output data frame
      out$Mean <- rep(NA, length(out$Index)) 
      
      for (k in 1:length(pt)) {
        for (i in 1:test_t_max[k]) {
          if (is.na(df_test$Bother[df_test$Patient == pt[k]][i])) {
            
            out$Mean[(out$Day == i) & (out$Patient == pt[k])] = rep(99, 11)
            
          } else {
          
            out$Mean[(out$Day == i) & (out$Patient == pt[k])] = rep(0, 11)
            observed <- df_test$Bother[df_test$Patient == pt[k]][i]
            
            for (j in (observed+1):11) {
              out$Mean[(out$Day == i) & (out$Patient == pt[k])][j] = 1
            }
          }
        }
      }
      
      # Calculate RPS
      pred <- par$Mean[startsWith(par$Variable, "cdf") & (par$Day %in% 1:t_pred)]
      both <- data.frame(out$Mean, pred)
      colnames(both) <- c("actual", "pred")
      both <- both[!(both$actual > 98 | both$pred > 98),]
      RPS <- sum((both$actual - both$pred)^2)
      RPS <- RPS / ((11-1) + length(both$pred))
      print(RPS)
      
      # Append to the output data frame
      output[(t+1), model] = RPS
      print(paste0("Finished ", t, " out of ", iterations, " for ", model))
    }
  }
  saveRDS(output, output_file)

  
} else if (metric == "LPD") {
  
  # Initialise output data frame
  output <- data.frame(matrix(vector(), (iterations+1), length(mdl_list), 
                              dimnames=list(c(), mdl_list)),
                       stringsAsFactors=F)
  output$iteration <- 0:iterations
  output_file <- file.path("Result_eval", paste0("summary_LPD.rds"))
  
  
  for (model in mdl_list) {
    # Loop over model
    for (t in 0:iterations) {
      # Loop over iterations
      suff <- paste0(model, "_eval.rds")  # Get file path of the parameters file
      par_file <- file.path("Result_eval", paste0("par", t, "_", suff))
      par <- readRDS(par_file) # Load the parameters file

      # Calculate log predictive density
      x <- par$Mean[(par$Variable == "lpd") & (par$Day %in% (1:t_pred))]
      x <- x[(x<99)]
      lpd <- mean(x)
      print(lpd)
      
      # Append to the output data frame
      output[(t+1), model] <- lpd
      print(paste0("Finished ", t, " out of ", iterations, " for ", model))
    }
  }
  saveRDS(output, output_file)
  
  
} else if (metric == "MSE") {
  
  # Initialise output data frame
  output <- data.frame(matrix(vector(), (iterations+1), length(mdl_list), 
                              dimnames=list(c(), mdl_list)),
                       stringsAsFactors=F)
  output$iteration <- 0:iterations
  output_file <- file.path("Result_eval", paste0("summary_MSE.rds"))
  
  for (model in mdl_list) {
    # Load the data for this model
    df <- preprocess_data(model)
    
    for (t in 0:iterations) {
      # Load test data for this iteration
      if (t == 0) {
        test_days <- 3:(3+(t_pred-1))
      } else {
        test_days <- (t*t_horizon + 1):(t*t_horizon + t_pred)
      }
      df_test <- df[df$Day %in% test_days,]
      df_test$Bother[is.na(df_test$Bother)] <- 99
      pt <- unique(df_test$Patient)
      
      suff <- paste0(model, "_eval.rds")   # Get file path of the parameters file
      par_file <- file.path("Result_eval", paste0("par", t, "_", suff))
      par <- readRDS(par_file) # Load the parameters file
      
      # Calculate MSE
      pred <- par$Mean[(par$Variable == "y_pred") & (par$Day %in% 1:t_pred)]
      actual <- df_test$Bother
      both <- data.frame(actual, pred)
      both <- both[!(both$actual > 98 | both$pred > 98),]
      MSE <- mean((both$actual - both$pred)^2)
      print(MSE)
      
      # Append to the output data frame
      output[(t+1), model] = MSE
      print(paste0("Finished ", t, " out of ", iterations, " for ", model))
    }
  }
  saveRDS(output, output_file)
}





