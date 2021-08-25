data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing Bother obs
  int<lower = 0> z_N_obs; // Number of non-missing Sassad obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing Bother obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing Bother obs
  int<lower = 1, upper = N> z_idx_obs[z_N_obs]; // Vector of indices of non-missing Sassad obs
  int<lower = 0, upper = 27> z_obs[z_N_obs]; // Non-missing Sassad obs
  real<lower = 0, upper = 1> Treat[N]; // Daily treatment usage
  
  int<lower = 2, upper = 10> M; // Upper bound of Bother
  int<lower = 2, upper = 27> D; // Upper bound of Sassad
}



transformed data {
  int yc_obs[N_obs]; // Categorical y_obs
  int zc_obs[z_N_obs]; // Categorial z_obs
  
  int start[N_pt]; // Index of first observation for each patient
  int end[N_pt]; // Index of last observation for each patient
  
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
  }
  
  for (i in 1:z_N_obs) {
    zc_obs[i] = z_obs[i] + 1;
  }
  
  for (k in 1:N_pt) {
    if (k == 1) {
      start[k] = 1;
    } else {
      start[k] = end[k-1] + 1;
    }
    end[k] = start[k] - 1 + t_max[k];
  }
}



parameters {
  real<lower = 0> sigma_ylat; // Standard deviation of y_lat
  real<lower = 0, upper = 10> y_lat[N]; // Latent score
  
  real<lower = 2, upper = 8> mu_y0; // Population mean of y_lat at t0
  real<lower = 0> sigma_y0; // Population standard deviation of y_lat at t0
  
  real mu_wY; // Population autocorrelation logit mean
  real<lower = 0> sigma_wY; // Population autocorrelation logit standard dev
  real eta_wY[N_pt]; // Non-centred parameterisation for autocorrelation
  
  real mu_wT; // Population treatment response mean
  real<lower = 0> sigma_wT; // Population treatment response standard dev
  real eta_wT[N_pt]; // Non-centred parameterisation for treatment response
  
  real<lower = 0> sigma_lam; // Population standard deviation of lambda
  real<lower = 0> eta_lam[N_pt]; // Non-centred parameterisation for lambda
  real<lower = 0> eta_R[N]; // Non-centred parameterisation for R
  
  real<lower = 0> y_c0; // First cutpoint for Bother
  real<lower = 0> y_delta[M-1]; // Difference between cutpoints for Bother
  real<lower = 0> z_c0; // First cutpoint for Sassad
  real<lower = 0> z_delta[D-1]; // Difference between cutpoints for Sassad
}



transformed parameters {
  vector[M] y_ct; // Cutpoints for Bother
  vector[D] z_ct; // Cutpoints for Sassad

  real wY[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient treatment response
  real lam[N_pt]; // Patient flare exponential distribution rate parameter
  real R[N]; // Flare intensity
  
  y_ct[1] = y_c0;
  for (i in 2:M) {
    y_ct[i] = y_ct[i-1] + y_delta[i-1];
  }
  
  z_ct[1] = z_c0;
  for (i in 2:D) {
    z_ct[i] = z_ct[i-1] + z_delta[i-1];
  }
  
  for (k in 1:N_pt) {
    wY[k] = inv_logit(mu_wY + sigma_wY*eta_wY[k]);
    wT[k] = mu_wT + sigma_wT*eta_wT[k];
    lam[k] = sigma_lam*eta_lam[k];
    
    for (t in (start[k]+1):end[k]) {
      R[t] = lam[k] * eta_R[t];
    }
  }
}



model {
  sigma_ylat ~ normal(0, 2);
  mu_y0 ~ normal(5, 2);
  sigma_y0 ~ normal(0, 2);
  
  mu_wY ~ normal(0, 1);
  sigma_wY ~ normal(0, 1.5);
  
  mu_wT ~ normal(0, 1);
  sigma_wT ~ normal(0, 0.5);
  
  sigma_lam ~ normal(0, 1);
  
  y_c0 ~ normal(0.5, 0.1);
  y_delta ~ normal(1, 1);
  z_c0 ~ normal(0.17, 0.04);
  z_delta ~ normal(0.35, 0.35);
  
  eta_wY ~ std_normal();
  eta_wT ~ std_normal();
  eta_lam ~ std_normal();
  eta_R ~ exponential(1);
  
  for (k in 1:N_pt) {
    y_lat[start[k]] ~ normal(mu_y0, sigma_y0);
    
    for (t in (start[k]+1):end[k]) {
      y_lat[t] ~ normal(wY[k] * y_lat[t-1] + wT[k]*Treat[t-1] + R[t], sigma_ylat);
    }
  }
  for (t in 1:N_obs) {
    yc_obs[t] ~ ordered_logistic(y_lat[idx_obs[t]], y_ct);
  }
  for (t in 1:z_N_obs) {
    zc_obs[t] ~ ordered_logistic(y_lat[z_idx_obs[t]], z_ct);
  }
}

