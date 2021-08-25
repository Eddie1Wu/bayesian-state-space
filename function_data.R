
# Load packages
pacman::p_load(dplyr, tidyr)


extract_signs_SWET <- function() {
  # Extract SASSAD area-signs in long format from csv data file
  #
  # Returns:
  # Data frame
  
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
  # Extract SASSAD total score in long format
  
  extract_signs_SWET() %>%
    filter(Sign %in% c("Erythema", "Exudation", "Excoriation", "Dryness", "Cracking", "Lichenification"),
           Area %in% c("Head_Neck", "Hands", "Arms", "Trunk", "Feet", "Legs")) %>%
    group_by(Patient, Week) %>%
    summarise(SASSAD = sum(Score)) %>%
    ungroup()
}


drop_patient <- function(df) {
  # Drop patients with more than 50% obs missing
  #
  # Returns:
  # Data frame
  
  t_max <- aggregate(Day ~ Patient, data = df, FUN = length)$Day
  drop_idx <- which(t_max<(16*7/2))
  patient_no <-unique(df$Patient)[drop_idx] # Get patient IDs
  df <- df[!(df$Patient %in% patient_no), ]
  
  # Remove irrelevant variables
  df <- df[, !(colnames(df) %in% c("Date", "CI", "SU", "Home"))]
  df <- rename(df, Treatment = CS) # Rename CS to Treatment
  df$Treatment[is.na(df$Treatment)] <- 0 # Impute missing treatment by 0
  df <- HuraultMisc::factor_to_numeric(df, "Treatment")
  return(df)
}


fill_missing_day1 <- function(df) {
  # Fill in missing Bother score for patient 1041, 3048 and 3085
  #
  # Return:
  # Data frame
  
  df$Bother[(df$Patient==1041) & (df$Day==1)] <- 7
  df$Bother[(df$Patient==3048) & (df$Day==1)] <- 1
  df$Bother[(df$Patient==3085) & (df$Day==1)] <- 2
  
  return(df)
}


merge_sassad <- function(df1, df2) {
  # Merge Bother with SASSAD 
  # By Patient, Day
  # 
  # Returns:
  # Data frame
  
  df_right <- df2 %>%
    mutate(Day = case_when(
      df2$Week == 0 ~ 1,
      df2$Week == 4 ~ 28,
      df2$Week == 12 ~ 84,
      df2$Week == 16 ~ 112
    ))
  
  # Remove irrelevant variables
  df_right <- df_right[, !(colnames(df_right) %in% c("Week"))]
  df_right$Patient <- as.factor(df_right$Patient)
  # Merge horizontally
  joint <- merge(df1, df_right, by = c("Patient", "Day"), all.x = TRUE)
  return(joint)
}


process1_SWET <- function(df) {
  # Processing SWET data for Bother score with ordinal logistic and no SASSAD
  #
  # df: SWET main dataset
  #
  # Returns:
  # Data frame
  
  df <- drop_patient(df)  # Remove patients with more than 50% obs missing
  df <- fill_missing_day1(df)
  return(df)
}


process2_SWET <- function(df1, df2) {
  # Processing SWET data for Bother with ordinal logistic and 
  # SASSAD with linear mixed effects
  #
  # df1: SWET main dataset
  # df2: SASSAD total score dataset
  #
  # Returns:
  # Data frame
  
  df_left <- drop_patient(df1) # Remove patients with more than 50% obs missing
  df_left <- fill_missing_day1(df_left)
  
  out <- merge_sassad(df_left, df2) # Merge horizontally
  
  # Remove Patient 3022 and 3023 who have no SASSAD measurement
  patient_no = c(3022, 3023)
  out <- out[!(out$Patient %in% patient_no), ]
  return(out)
}


process3_SWET <- function(df1, df2) {
  # Processing SWET data for Bother with ordinal logistic and
  # SASSAD with ordinal logistic in bin size = 3, 28 groups in total from 0-27
  #
  # df1: SWET main dataset
  # df2: SASSAD total score dataset
  #
  # Returns:
  # Data frame
  
  df_left <- drop_patient(df1) # Remove patients with more than 50% obs missing
  df_left <- fill_missing_day1(df_left)
  
  df2$SASSAD[df2$SASSAD>81] = 81  # Put patients with SASSAD > 81 in one group
  df2$SASSAD <- ceiling(df2$SASSAD/3)  # Put patients into bins of 3
  
  out <- merge_sassad(df_left, df2)  # Merge horizontally
  
  # Remove Patient 3022 and 3023 who have no SASSAD measurement
  patient_no = c(3022, 3023)
  out <- out[!(out$Patient %in% patient_no), ]
  return(out)
}


format_stan_data_RW <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            M = 10))
}


format_stan_data_AR <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            y0 = Bother[Day == 1],
            Treat = Treatment,
            M = 10))
}


format_stan_data_MC <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            Treat = Treatment,
            M = 10))
}


format_stan_data_LmeAR <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            z_N_obs = sum(!is.na(SASSAD)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            z_t_max = aggregate(SASSAD ~ Patient, FUN = length)$SASSAD,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            z_idx_obs = which(!is.na(SASSAD)),
            z_obs = na.omit(SASSAD),
            y0 = Bother[Day == 1],
            Treat = Treatment,
            M = 10))
}


format_stan_data_LmeMC <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            z_N_obs = sum(!is.na(SASSAD)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            z_t_max = aggregate(SASSAD ~ Patient, FUN = length)$SASSAD,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            z_idx_obs = which(!is.na(SASSAD)),
            z_obs = na.omit(SASSAD),
            Treat = Treatment,
            M = 10))
}


format_stan_data_OrdMC <- function(df) {
  with(df,
       list(N = length(Bother),
            N_obs = sum(!is.na(Bother)),
            z_N_obs = sum(!is.na(SASSAD)),
            N_pt = length(unique(Patient)),
            t_max = aggregate(Day ~ Patient, FUN = length)$Day,
            idx_obs = which(!is.na(Bother)),
            y_obs = na.omit(Bother),
            z_idx_obs = which(!is.na(SASSAD)),
            z_obs = na.omit(SASSAD),
            Treat = Treatment,
            M = 10,
            D = 27))
}


fill_a_patient <- function(df) {
  # Fill in the missing days for a patient with NA
  
  Patient <- factor(rep(5049, 4))
  Day <- 68:71
  Treatment <- rep(0, 4)
  Bother <- rep(NA, 4)
  x <- data.frame(Patient, Day, Treatment, Bother)
  out <- dplyr::bind_rows(df, x)
  
  return(out)
}


drop_missing_SASSAD <- function(df) {
  # Drop patients whose first SASSAD measure is missing
  
  for (i in unique(df$Patient)) {
    if (is.na((df$SASSAD[(df$Patient == i) & (df$Day == 1)]))) {
      df <- df[!(df$Patient == i),]
    }
  }
  
  return(df)
}


drop_missing_SASSADs <- function(df) {
  # Drop patients whose SASSAD measures are incomplete in the 16 weeks
  
  days <- c(1,28,84,112)
  for (i in unique(df$Patient)) {
    if (any(is.na(df$SASSAD[(df$Patient == i) & (df$Day %in% days)]))) {
      df <- df[!(df$Patient == i),]
    }
    if (length(df$SASSAD[(df$Patient == i) & (df$Day %in% days)]) < 4) {
      df <- df[!(df$Patient == i),]
    }
  }
  
  return(df)
}






