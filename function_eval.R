
# Extract parameters ----------------------------------------------------------

get_index <- function(pt, t_max, multiple) {
  # Translate unique indices (used in Stan model) into (patient, day) pair
  #
  # Args:
  # pt: Vector of patient ID (same order as the patient parameters in stanfit)
  # t_max: Vector of time-series length (including missing values) for each patient
  # multiple: Repeat each row of the data frame how many times
  #
  # Returns:
  # Data frame with columns: Patient, Day, Index
  
  stopifnot(length(pt) == length(t_max),
            is.numeric(t_max),
            all(t_max == round(t_max)))
  
  out <- data.frame(Patient = rep(pt, t_max),
                    Day = do.call(c, lapply(t_max, function(x) {1:x})))
  out[["Patient"]] <- as.character(out[["Patient"]])
  out <- out[rep(seq_len(nrow(out)), each = multiple), ]
  out[["Index"]] <- 1:nrow(out)
  return(out)
}


extract_parameters <- function(fit, param, param_ind, param_latent, pt, data_stan) {
  # Extract parameters' summary 
  #
  # Args:
  # fit: stanfit object
  # param: all params to extract
  # param_ind: individual parameters
  # param_latent: latent parameters
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the Stan function
  #
  # Returns: data frame containing posterior summary statistics of parameters 

  par <- HuraultMisc::summary_statistics(fit, param)
  par$Patient <- NA
  par$Day <- NA
  
  pt <- as.character(pt)
  
  ## Patient-dependent parameter
  for (i in intersect(param_ind, param)) {
    idx <- which(par$Variable == i)
    par$Patient[idx] <- pt[par$Index[idx]]
  }
  
  ## Patient and time-dependent parameter (latent states)
  dict <- get_index(pt, data_stan$t_max, 1)
  for (i in intersect(param_latent, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Day")] <- dict[, c("Patient", "Day")]
  }
  
  return(par)
}


extract_parameters_eval <- function(fit, param, param_ind, param_latent,
                                    param_pred, pt, data_stan) {
  # Extract parameters' summary for forward chaining
  #
  # Args:
  # fit: stanfit object
  # param: all params to extract
  # param_ind: individual parameters
  # param_latent: latent parameters
  # param_pred: the generated quantities
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the Stan function
  #
  # Returns: data frame containing posterior summary statistics of parameters 
  
  par <- HuraultMisc::summary_statistics(fit, param)
  par$Patient <- NA
  par$Day <- NA
  
  # Assign index to the Variable "cdf"
  x <- length(is.na(par$Index[startsWith(par$Variable, "cdf")]))
  par$Index[is.na(par$Index) & startsWith(par$Variable, "cdf")] <- 1:x
  pt <- as.character(pt)
  
  ## Patient-dependent parameter
  for (i in intersect(param_ind, param)) {
    idx <- which(par$Variable == i)
    par$Patient[idx] <- pt[par$Index[idx]]
  }
  
  ## Patient and time-dependent parameter (latent states)
  dict <- get_index(pt, data_stan$t_max, 1)
  for (i in intersect(param_latent, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Day")] <- dict[, c("Patient", "Day")]
  }
  
  ## Predictions
  dict1 <- get_index(pt, data_stan$test_t_max, 1)
  dict2 <- get_index(pt, data_stan$test_t_max, 11)
  for (i in intersect(param_pred, param)) {
    if (startsWith(i, "cdf")) {
      idx <- sort(which(startsWith(par$Variable, "cdf")))
      par[idx, c("Patient", "Day")] <- dict2[, c("Patient", "Day")]
    } else {
      idx <- sort(which(par$Variable == i))
      par[idx, c("Patient", "Day")] <- dict1[, c("Patient", "Day")]
    }
  }
  
  return(par)
}


extract_parameters_check <- function(fit, param, param_ind, param_latent,
                                    param_pred, pt, data_stan) {
  # Extract parameters' summary for forward chaining
  #
  # Args:
  # fit: stanfit object
  # param: all params to extract
  # param_ind: individual parameters
  # param_latent: latent parameters
  # param_pred: the generated quantities
  # pt: vector of patients ID (same order as the patient parameters in stanfit)
  # data_stan: data input to the Stan function
  #
  # Returns: data frame containing posterior summary statistics of parameters 
  
  par <- HuraultMisc::summary_statistics(fit, param)
  par$Patient <- NA
  par$Day <- NA

  pt <- as.character(pt)
  
  ## Patient-dependent parameter
  for (i in intersect(param_ind, param)) {
    idx <- which(par$Variable == i)
    par$Patient[idx] <- pt[par$Index[idx]]
  }
  
  ## Patient and time-dependent parameter (latent states)
  dict <- get_index(pt, data_stan$t_max, 1)
  for (i in intersect(param_latent, param)) {
    idx <- sort(which(par$Variable == i))
    par[idx, c("Patient", "Day")] <- dict[, c("Patient", "Day")]
  }
  
  ## Predictions
  dict1 <- get_index(pt, data_stan$t_max, 1)
  for (i in intersect(param_pred, param)) {
    if (i == "sassad_pred") {
      dict1 <- get_index(pt, c(1,1), 1)
      idx <- sort(which(par$Variable == i))
      par[idx, c("Patient", "Day")] <- dict1[, c("Patient", "Day")]
    } else {
      idx <- sort(which(par$Variable == i))
      par[idx, c("Patient", "Day")] <- dict1[, c("Patient", "Day")]
    }
  }
  
  return(par)
}



# Functions for historical model forward-chaining -----------------------------

calculate_lpd_and_ypred <- function(df_train, df_test, pt) {
  # Calculate log predictive density based on historical occurrences
  #
  # Args:
  # df_train: training set
  # df_test: test set
  # pt: vector of patient IDs
  #
  # Return: data frame in the same format as Stan param output
  
  ## Initialisation
  # Find maximum number of days for each patient in training and test set
  t_max <- aggregate(df_train$Day ~ df_train$Patient, FUN = length)$`df_train$Day`
  test_t_max <- aggregate(df_test$Day ~ df_test$Patient, FUN = length)$`df_test$Day`
  
  # Get indices of start and end day of each patient
  start <- rep(NA, length(pt))
  end <- rep(NA, length(pt))
  test_start <- rep(NA, length(pt))
  test_end <- rep(NA, length(pt))
  
  for (k in 1:length(pt)) {
    if (k == 1) {
      start[k] = 1
      test_start[k] = 1
    } else {
      start[k] = end[k-1] + 1
      test_start[k] = test_end[k-1] + 1
    }
    end[k] = start[k] - 1 + t_max[k]
    test_end[k] = test_start[k] - 1 + test_t_max[k]
  }
  
  out <- get_index(pt, test_t_max, 1) # Create an empty output data frame
  out$Variable <- rep("lpd", length(out$Index))
  out$Mean <- rep(NA, length(out$Index))
  
  out_y_pred <- get_index(pt, test_t_max, 1) # Create an empty output data frame
  out_y_pred$Variable <- rep("y_pred", length(out_y_pred$Index))
  out_y_pred$Mean <- rep(NA, length(out_y_pred$Index))
  
  for (k in 1:length(pt)) {
    # Looping over patient
    x <- df_train$Bother[start[k]:end[k]]
    x <- x[!(is.na(x))]
    data_vector <- c(x, 0:10)
    pmf <- rep(NA, 11)
    y <- length(data_vector)
    for (i in 1:11) {
      pmf[i] <- sum(data_vector == (i-1)) / y
    }

    for (t in test_start[k]:test_end[k]) {
      # Looping over test time horizon
      if (df_test$Bother[t] == 99) {
        out$Mean[t] <- 99
        out_y_pred$Mean[t] <- 99
      } else {
        p <- pmf[(df_test$Bother[t] + 1)]
        out$Mean[t] <- log(p) # Find log(p) based on historical frequency
        
        out_y_pred$Mean[t] <- sum(pmf * (0:10)) # Find expected value of y_pred
      }
    }
  }
  
  out <- dplyr::bind_rows(out, out_y_pred)
  return(out)
}


calculate_cdf <- function(df_train, df_test, pt) {
  # Calculate cdf based on historical occurrences
  #
  # Args:
  # df_train: training set
  # df_test: test set
  # pt: vector of patient IDs
  #
  # Return: data frame in the same format as Stan param output
  
  ## Initialisation
  # Find maximum number of days for each patient in training and test set
  t_max <- aggregate(df_train$Day ~ df_train$Patient, FUN = length)$`df_train$Day`
  test_t_max <- aggregate(df_test$Day ~ df_test$Patient, FUN = length)$`df_test$Day`
  
  # Get indices of start and end day of each patient
  start <- rep(NA, length(pt))
  end <- rep(NA, length(pt))
  test_start <- rep(NA, length(pt))
  test_end <- rep(NA, length(pt))
  
  for (k in 1:length(pt)) {
    if (k == 1) {
      start[k] = 1
      test_start[k] = 1
    } else {
      start[k] = end[k-1] + 1
      test_start[k] = test_end[k-1] + 1
    }
    end[k] = start[k] - 1 + t_max[k]
    test_end[k] = test_start[k] - 1 + test_t_max[k]
  }
  
  out <- get_index(pt, test_t_max, 11) # Create an empty output data frame
  out$Variable <- rep("cdf", length(out$Index))
  out$Mean <- rep(NA, length(out$Index))
  
  for (k in 1:length(pt)) {
    # Looping over patient
    x <- df_train$Bother[start[k]:end[k]]
    x <- x[!(is.na(x))]
    data_vector <- c(x, 0:10)
    pmf <- rep(NA, 11)
    y <- length(data_vector)
    for (i in 1:11) {
      pmf[i] <- sum(data_vector == (i-1)) / y
    }
    
    for (t in 1:test_t_max[k]) {

      out$Mean[(out$Day == t) & (out$Patient == pt[k])][1] = pmf[1]
      
      for (i in 2:11) {
        
        out$Mean[(out$Day == t) & (out$Patient == pt[k])][i] = 
          out$Mean[(out$Day == t) & (out$Patient == pt[k])][i-1] + pmf[i]
      }
    }
  }
  return(out)
}


# Functions for prior predictive check and fake data checl --------------------

generate_treatment <- function(p, tmax) {
  # Generate treatment from a markov chain
  #
  # Args:
  # p: vector of length 2 containing the probabilities of using treatment when treatment was not used (p01) or used (p11) the day before respectively
  # To simulate "sticky" (realistic behaviour), p[1] << p[2]
  # tmax: length of time series
  #
  # Returns:
  # Time series of treatment usage
  
  stopifnot(is.numeric(p),
            length(p) == 2,
            min(p) >= 0,
            max(p) <= 1,
            is.numeric(tmax),
            length(tmax) == 1,
            tmax > 0,
            tmax == round(tmax))
  
  Treat <- rep(0, tmax)
  for (i in 2:tmax) {
    if (Treat[i - 1] == 0){
      Treat[i] <- rbinom(1, 1, p[1])
    } else {
      Treat[i] <- rbinom(1, 1, p[2])
    }
  }
  return(Treat)
}


compute_pmf <- function(ps, pred = FALSE) {
  # Compute probability mass function from matrix of posterior samples:
  # - Truncate: discard value outside 0-10
  # - (Truncate trajectories: discard future when past is outside 0-10)
  # - Discretise
  # - Correct rounding at the bounds
  # - Add artificial sample (to avoid probability to be exactly 0)
  #
  # Args:
  # ps: matrix of posterior samples
  # pred: logical indicating whether we are in "prediction mode"
  #
  # Returns:
  # Matrix: observations in rows, pmf in columns
  
  stopifnot(is.matrix(ps),
            is.logical(pred))
  
  # Truncate values outside 0-10
  ps[ps < 0 | ps > 10] <- NA # Truncate value outside 0-10
  
  # If prediction mode, only accept samples when past samples are in 0-10 as well
  if (pred & ncol(ps) > 1) {
    for (i in 1:(ncol(ps) - 1)) {
      ps[is.na(ps[, i]), (i + 1):(ncol(ps))] <- NA
    }
  }
  
  # Discretise/Round
  ps <- round(ps)
  
  # Compute probability table
  pmf <- do.call(rbind,
                 lapply(1:ncol(ps),
                        function(i) {
                          x <- na.omit(ps[, i])
                          n_samp <- length(x)
                          x <- c(x, x[x %in% c(0, 10)]) #  Double the number of samples at the bounds (correct for half support)
                          x <- c(x, 0:10) # Add artificial samples to avoid probability of 0
                          tbl <- table(x) / length(x) # Probability table
                          # If original x doesn't contain enough samples return NA
                          if (n_samp <= 100) {
                            tbl <- tbl * NA
                          }
                          return(tbl)
                        }))
  return(pmf)
}


prepare_ppc <- function(fit, df, par, idx_pred) {
  # Prepare a dataset for posterior predictive checks plot
  #
  # Args:
  # fit: Stanfit object
  # df: Dataframe of the data
  # par: Dataframe of parameters' summary
  # idx_pred: Dataframe indexing predictions (cf. predictions_dictionary)
  #
  # Returns:
  # Dataframe with columns: Patient, Day, Severity (data), Fit (Posterior fit), 0 to 10 (Predictions probability)
  
  stopifnot(class(fit) == "stanfit",
            is.data.frame(df),
            all(c("Patient", "Day", "Bother") %in% colnames(df)),
            is.data.frame(par),
            all(c("Patient", "Day", "Mean", "Variable") %in% colnames(par)),
            nrow(subset(par, Variable == "y_pred")) > 0,
            is.data.frame(idx_pred),
            all(c("Patient", "Day", "Index") %in% colnames(idx_pred)))
  
  stopifnot(nrow(par) > 0)
  
  # Compute predictions
  yrep <- rstan::extract(fit, pars = "y_pred")[[1]]
  prob <- compute_pmf(yrep, pred = FALSE)
  pred <- cbind(idx_pred, prob)
  
  
  
  # Merge df and fit (in par)
  out <- merge(df[, c("Patient", "Day", "Bother")],
               par[par[["Variable"]] == "y_pred", c("Patient", "Day", "Mean")],
               by = c("Patient", "Day"), all = TRUE)
  # Merge results with predictions
  out <- merge(out, pred, by = c("Patient", "Day"), all.x = TRUE, all.y = FALSE)
  # Minor processing
  out <- out[order(out[["Patient"]], out[["Day"]]), ]
  out <- change_colnames(out, "Mean", "Fit")
  
  return(out)
}



plot_ppc <- function(ppc, patientID) {
  # Plot Posterior predictive distribution (density/mass plot)
  #
  # Args:
  # ppc: Dataframe, output from prepare_ppc
  # patientID: Patient ID
  #
  # Returns:
  # Ggplot
  
  library(ggplot2)
  palette <- c("#FFFFFF", RColorBrewer::brewer.pal(n = 6, "Blues")) # white-blue
  gamma <- 1.5
  col_correction <- scales::rescale(seq(0, 1, .1)^gamma)
  
  stopifnot(is.data.frame(ppc),
            all(c("Patient", "Day", "Bother", "Fit", 0:10) %in% colnames(ppc)),
            patientID %in% unique(ppc[["Patient"]]))

  tmp <- subset(ppc, Patient == as.character(patientID))
  
  ## Process trajectories (observed and fit)
  traj <- tmp[, c("Patient", "Day", "Bother", "Fit")]
  
  # If we have a pattern of points Observed Missing Observed (O-M-O), the O-M line will be marked as Observed and the M-O line will be marked as Missing
  # instead, we want both line to be marked as Missing
  traj[["Observed"]] <- !is.na(traj[["Bother"]])
  traj$ColourLine <- as.logical(traj$Observed)
  traj$ColourLine[which(!as.logical(traj$ColourLine)) - 1] <- FALSE
  traj$ColourLine <- factor(traj$ColourLine, levels = c(TRUE, FALSE))
  
  # Only show fit when missing
  traj$Fit[traj$Observed] <- round(traj$Fit[traj$Observed])
  # Offset by 1 to align with factors
  traj$Fit <- traj$Fit + 1
  
  ## Process predictions
  pred <- tmp[, c("Patient", "Day", 0:10)]
  pred<- reshape2::melt(pred,
                        id.vars = c("Patient", "Day"),
                        variable.name = "S",
                        value.name = "Probability")
  
  ## Plot
  # Heatmap
  p <- ggplot() +
    geom_tile(data = pred, aes(x = Day, y = S, fill = Probability)) +
    scale_fill_gradientn(colours = palette, limits = c(0, 1), breaks = c(0, .5, 1), values = col_correction) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
  # Overlay trajectory
  p <- p +
    geom_line(data = traj, aes(x = Day, y = Fit, colour = ColourLine, group = 1), lwd = 1.5) + # Trajectory (observed and missing)
    geom_point(data = traj[stats::filter(as.numeric(as.logical(traj$Observed)), c(-1, 1, -1), method = "convolution", sides = 2) == 1, ],
               aes(x = Day, y = Fit), colour = "black") + # show the observed points in the case of a pattern M-O-M
    scale_color_manual("", labels = c("Observed", "Missing"), values = c("black","grey"))
  # Formatting
  p <-  p +
    labs(y = "Bother score") +
    theme_classic(base_size = 20) +
    theme(legend.position = "top")
  
  return(p)
}




