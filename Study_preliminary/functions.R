
pacman::p_load(dplyr, tidyr)

# Functions for processing raw data --------------------------------------------

extract_signs_SWET <- function() {
  # Define the extract_signs_SWET function to get SASSAD
  file <- get_datapath("SWET/sassad_anonymised.csv")
  df <- readr::read_csv(file)
  
  df$visdes <- factor(df$visdes,
                      levels = c("eligibility criteria and baseline",
                                 "four week crf", "12 week crf", "16 week crf"),
                      labels = c(0, 4, 12, 16))
  df <- rename(df, Patient = refno, Week = visdes)
  
  delete <- c("trialref", "siteid", "patno", "visno", "pageno", "repsite", "photo", "verify", "verify2", "extstamp")
  for (i in delete) {df[, i] <- NULL}
  
  df <- reshape2::melt(df, id = c("Patient", "Week"), value.name = "Score")
  
  tmp <- as.data.frame(do.call(rbind, strsplit(as.character(df$variable), "_")))
  colnames(tmp) <- c("Area", "Sign")
  
  df <- cbind(df, tmp)
  df$variable <- NULL
  
  df$Area <- factor(df$Area, levels = c("hn", "ha", "ar", "tr", "fe", "le", "bd", "rs"),
                    labels = c("Head_Neck", "Hands", "Arms", "Trunk", "Feet", "Legs", "Body", "Representative_Site"))
  df$Sign <- factor(df$Sign, levels = c("ery", "exu", "exc", "dry", "cra", "lic", "oed", "tot"),
                    labels = c("Erythema", "Exudation", "Excoriation", "Dryness", "Cracking", "Lichenification", "Oedema", "Total"))
  return(df)
}



extract_SASSAD_SWET <- function() {
  # Define the extract_SASSAD_SWET function
  extract_signs_SWET() %>%
    filter(Sign %in% c("Erythema", "Exudation", "Excoriation", "Dryness", "Cracking", "Lichenification"),
           Area %in% c("Head_Neck", "Hands", "Arms", "Trunk", "Feet", "Legs")) %>%
    group_by(Patient, Week) %>%
    summarise(SASSAD = sum(Score)) %>%
    ungroup()
}



# Functions for pre-processing SWET data ---------------------------------------

reshape_SWET_to_wide <- function(df) {
  # Reshape SWET into wide format for merging with SASSAD
  # Arg: SWET dataset from TanakaData
  # 
  # Return: 
  # Wide data frame of SWET
  
  out <- df %>% # Generate variables for reshaping to wide
    mutate(Week = case_when(
      df$Day == 1 ~ 1,
      df$Day == 29 ~ 2,
      df$Day == 57 ~ 3,
      df$Day == 85 ~ 4,
      df$Day == 113 ~ 5
    )) %>%
    fill(Week) %>%  
    group_by(Patient, Week) %>%  
    mutate(Day_Count = 1 + 1*(0:(length(Week)-1))) %>%   # Create the day counts for every four weeks
    ungroup()
  
  out <- out[, !names(out) %in% c("Day")] # Drop "Day" variable
  
  out <- out %>%
    group_by(Week) %>%
    pivot_wider(
      id_cols = c(Patient, Week),
      names_from = Day_Count,
      values_from = Bother,
      names_prefix = "Bother."
    ) %>%
    ungroup()
  
  out$Week <- factor(out$Week, # Convert numeric into factor to match the SASSAD dataset
                     levels = c(1, 2, 3, 4, 5),
                     labels = c(4, 8, 12, 16, 20))
  
  out$Patient <- as.factor(out$Patient)  # Convert numeric into factor
  
  return(out)
}



reshape_sassad_to_wide <- function(df) {
  # Reshape sassad into wide format and generate relevant variables
  #
  # Return:
  # data frame of wide sassad
  df <- extract_signs_SWET()
  df <- subset(df, Area != "Representative_Site" & Area != "Body")
  df$Score[is.na(df$Score)] <- 0  # Fill NA with 0
  
  # Concatenate the Area and Sign variables
  df$AreaSign <- paste(df$Area, df$Sign) 
  df <- df[-c(4:5)]  # Drop Area and Sign
  
  # Reshape into wide panel
  out <- df %>%
    group_by(Week) %>%
    pivot_wider(
      id_cols = c(Patient, Week),
      names_from = AreaSign,
      values_from = Score
    ) %>%
    ungroup()
  
  # Calculate the sign-wise sums
  out$Erythema_Total <- rowSums(dplyr::select(out, ends_with("Erythema")))
  out$Exudation_Total <- rowSums(dplyr::select(out, ends_with("Exudation")))
  out$Excoriation_Total <- rowSums(dplyr::select(out, ends_with("Excoriation")))
  out$Dryness_Total <- rowSums(dplyr::select(out, ends_with("Dryness")))
  out$Cracking_Total <- rowSums(dplyr::select(out, ends_with("Cracking")))
  out$Lichenification_Total <- rowSums(dplyr::select(out, ends_with("Lichenification")))
  
  # Only keep the relevant variables
  out <- dplyr::select(out, starts_with("Patient") | starts_with("Week") | ends_with("Total"))
  
  out$SASSAD <- rowSums(dplyr::select(out, ends_with("Total"))) / 2  # Calculate total
  out$Patient <- as.factor(out$Patient)  # Convert into factors for merging

  return(out)
}



# Functions for cross validation -----------------------------------------------

make_folds <- function(df, fold, by_patient = TRUE) {
  # Create folds for cross validation
  # 
  # Return:
  # A list containing the IDs for each fold
  
  if (by_patient) {
    ID <- unique(df$Patient)
  } else {
    ID <- df$ID
  }
  
  folds <- caret::createFolds(ID, k = fold, list = TRUE, returnTrain = TRUE)
  return(folds)
}



make_dataframe <- function(result, out, model) {
  # Create data frame for plotting
  #
  # Return:
  # Data frame with Model, RMSE, MAE, R2
  
  df <- data.frame(out[1:10])
  names(df) <- c("RMSE")
  df$MAE <- out[11:20]
  df$R2 <- out[21:30]
  df$Model <- model
  
  out <- bind_rows(result, df)
  
  return(out)
}



square_table <- function(x,y) {
  # Always create a square table
  #
  # Return:
  # A square table
  x <- factor(x)
  y <- factor(y)
  
  commonLevels <- sort(unique(c(levels(x), levels(y))))

  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  
  return(table(x,y))
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



make_dataframe2 <- function(result, out, model) {
  # Create data frame for plotting
  #
  # Return:
  # Data frame with Model, RMSE and Accuracy
  
  df <- data.frame(out[1:10])
  names(df) <- c("RMSE")
  df$Accuracy <- out[11:20]
  df$Model <- model
  
  out <- bind_rows(result, df)
  
  return(out)
}



cross_validate_1 <- function(model, df, folds, by_patient = TRUE) {
  # Cross validation for checking non-linearity in SASSAD against Bother
  #
  # Return:
  # A list containing RMSE, MAE and R2
  
  RMSE <- rep(NA, length(folds))
  MAE <- rep(NA, length(folds))
  Rsquare <- rep(NA, length(folds))
  
  for (i in 1:length(folds)) {
    
    # Split training and test set
    if (by_patient) {
      patient_ID <- unique(df$Patient)
      train_id <- patient_ID[folds[[i]]]
      test_id <- patient_ID[!(patient_ID %in% train_id)]
      train <- subset(df, Patient%in%train_id)
      test <- subset(df, Patient%in%test_id)
    } else {
      train_id <- folds[[i]]
      test_id <- df$ID[!(df$ID %in% folds[[i]])]
      train <- subset(df, ID%in%train_id)
      test <- subset(df, ID%in%test_id)
    }
    
    # Train model 
    if (model == "intercept only") {
      fit <- lm(SASSAD ~ 1, data = train)
    } else if (model == "linear") {
      fit <- lm(SASSAD ~ Bother.28, data = train)
    } else if (model == "polynomial") {
      fit <- lm(SASSAD ~ poly(Bother.28, 4), data = train)
    } else if (model == "natural splines") {
      fit <- lm(SASSAD ~ ns(Bother.28, 3), data = train)
    } else if (model == "I-splines") {
      fit <- lm(SASSAD ~ iSpline(Bother.28, knots = c(5,7), degree = 3), data = train)
    } else if (model == "mixed slope") {
      fit <- lmer(SASSAD ~ Bother.28 + (Bother.28 | Patient), data = train)
    } else if (model == "mixed intercept") {
      fit <- lmer(SASSAD ~ Bother.28 + (1 | Patient), data = train)
    }
    
    # Make prediction and calculate metrics
    prediction <- predict(fit, newdata = test, allow.new.levels = TRUE)
    error <- caret::RMSE(prediction, test$SASSAD)
    abs_error <- caret::MAE(prediction, test$SASSAD)
    r2 <- caret::R2(prediction, test$SASSAD)
    
    RMSE[i] <- error
    MAE[i] <- abs_error
    Rsquare[i] <- r2
  }
  
  return(c(RMSE, MAE, Rsquare))
}



cross_validate_2 <- function(model, df, folds, factor = TRUE, by_patient = TRUE) {
  # Cross validation for classification to predict Bother
  #
  # Return:
  # A list containing RMSE and Accuracy

  if (factor) {
    df$Bother <- as.factor(df$Bother.28)
  } else {
    df$Bother <- df$Bother.28
  }
  
  RMSE <- rep(NA, length(folds))
  Accuracy <- rep(NA, length(folds))

  for (i in 1:length(folds)) {
    
    # Split training and test set
    if (by_patient) {
      patient_ID <- unique(df$Patient)
      train_id <- patient_ID[folds[[i]]]
      test_id <- patient_ID[!(patient_ID %in% train_id)]
      train <- subset(df, Patient%in%train_id)
      test <- subset(df, Patient%in%test_id)
    } else {
      train_id <- folds[[i]]
      test_id <- df$ID[!(df$ID %in% folds[[i]])]
      train <- subset(df, ID%in%train_id)
      test <- subset(df, ID%in%test_id)
    }
    
    # Train model 
    if (model == "naiveBayes") {
      fit <- naive_bayes(Bother ~ Erythema_Total + Exudation_Total 
                         + Excoriation_Total + Dryness_Total + Cracking_Total 
                         + Lichenification_Total, data = train, laplace = 3)
    } else if (model == "naiveBayes(sum)") {
      fit <- naive_bayes(Bother ~ SASSAD, data = train, laplace = 3)
    } else if (model == "randomForest") {
      fit <- randomForest(Bother ~ Erythema_Total + Exudation_Total 
                          + Excoriation_Total + Dryness_Total + Cracking_Total 
                          + Lichenification_Total, data = train)
    } else if (model == "randomForest(sum)") {
      fit <- randomForest(Bother ~ SASSAD, data = train)
    } else if (model == "multi-logistic") {
      fit <- multinom(Bother ~ Erythema_Total + Exudation_Total 
                      + Excoriation_Total + Dryness_Total + Cracking_Total 
                      + Lichenification_Total, data = train)
    } else if (model == "multi-logistic(sum)") {
      fit <- multinom(Bother ~ SASSAD, data = train)
    } else if (model == "ordinal logistic") {
      fit <- polr(Bother ~ Erythema_Total + Exudation_Total + Excoriation_Total
                  + Dryness_Total + Cracking_Total + Lichenification_Total, 
                  data = train, Hess = TRUE)
    } else if (model == "ordinal logistic(sum)") {
      fit <- polr(Bother ~ SASSAD, data = train, Hess = TRUE)
    }

    # Make prediction and calculate metrics
    if (model == "multi-logistic") {
      prediction <- predict(fit, newdata = test, "class")
    } else {
      prediction <- predict(fit, test)
    }
    
    tab <- square_table(test$Bother, prediction) # Building classification table
    RMSE[i] <- calc_table_rmse(tab, dim(tab)[1]) # Calculating RMSE
    Accuracy[i] <- round((sum(diag(tab))/sum(tab))*100, 2) # Calculating accuracy
  }
  
  return(c(RMSE, Accuracy))
}






