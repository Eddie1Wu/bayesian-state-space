data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing Bother obs
  int<lower = 0> z_N_obs; // Number of non-missing Sassad obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 0> z_t_max[N_pt]; // Time series length of SASSAD for each patient
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing Bother obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing Bother obs
  int<lower = 1, upper = N> z_idx_obs[z_N_obs]; // Vector of indices of non-missing Sassad obs
  int<lower = 0> z_obs[z_N_obs]; // Non-missing Sassad obs
  real<lower = 0, upper = 1> Treat[N]; // Daily treatment usage
  
  int<lower = 2, upper = 10> M; // Upper bound of Bother
  
  int<lower = 0> N_test; // Total number of test obs
  int<lower = 0> test_t_max[N_pt];
  int<lower = 0> y_test[N_test];
}



transformed data {
  int yc_obs[N_obs]; // Categorical y_obs
  int start[N_pt]; // Index of first observation for each patient
  int end[N_pt]; // Index of last observation for each patient
  int z_start[N_pt]; // Index of first SASSAD ob for each patient
  int z_end[N_pt]; // Index of last SASSAD ob for each patient
  
  int test_start[N_pt]; 
  int test_end[N_pt];
  
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
  }
  
  for (k in 1:N_pt) {
    if (k == 1) {
      start[k] = 1;
      z_start[k] = 1;
      test_start[k] = 1;
    } else {
      start[k] = end[k-1] + 1;
      z_start[k]= z_end[k-1] + 1;
      test_start[k] = test_end[k-1] + 1;
    }
    end[k] = start[k] - 1 + t_max[k];
    z_end[k] = z_start[k] - 1 + z_t_max[k];
    test_end[k] = test_start[k] - 1 + test_t_max[k];
  }
}



parameters {
  real<lower = 0> sigma_ylat; // Standard deviation of y_lat
  real b_ylat; // Intercept of y_lat
  real eta[N]; // Non-centred parameterisation for the gaussian process
  
  real<lower = 3, upper = 7> mu_y0; // Population mean of y_lat at t0
  real<lower = 0> sigma_y0; // Population standard deviation of y_lat at t0
  real<lower = 0, upper = 10> y0[N_pt]; // Initial latent score for each patient
  
  real mu_wY; // Population autocorrelation logit mean
  real<lower = 0> sigma_wY; // Population autocorrelation logit standard dev
  real eta_wY[N_pt]; // Non-centred parameterisation for autocorrelation
  
  real mu_wT; // Population treatment response mean
  real<lower = 0> sigma_wT; // Population treatment response standard dev
  real eta_wT[N_pt]; // Non-centred parameterisation for treatment response
  
  real<lower = 0> sigma_lam; // Population standard deviation of lambda
  real<lower = 0> eta_lam[N_pt]; // Non-centred parameterisation for lambda
  real<lower = 0> eta_R[N]; // Non-centred parameterisation for R
  
  real<lower = 0> c0; // First cutpoint for Bother
  real<lower = 0> delta[M-1]; // Difference between cutpoints for Bother
  
  real wS; // Coefficient of Sassad measurement process
  real mu_bs; // Population mean of Sassad measurement intercept
  real<lower = 0> sigma_bs; // Population standard deviation of Sassad intercept
  real eta_bs[N_pt]; // Non-centred parameterisation for Sassad intercept
  real<lower = 0> sigma_S; // Standard deviation of Sassad measurement
}



transformed parameters {
  vector[M] ct; // Cutpoints for Bother
  vector[N] y_lat; // Latent score

  real wY[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient treatment response
  real lam[N_pt]; // Patient flare exponential distribution rate parameter
  real R[N]; // Flare intensity
  real b_s[N_pt]; // Sassad measurement intercept
  
  ct[1] = c0;
  for (i in 2:M) {
    ct[i] = ct[i-1] + delta[i-1];
  }
  
  for (k in 1:N_pt) {
    wY[k] = inv_logit(mu_wY + sigma_wY*eta_wY[k]);
    wT[k] = mu_wT + sigma_wT*eta_wT[k];
    lam[k] = 1 + sigma_lam*eta_lam[k];
    b_s[k] = mu_bs + sigma_bs*eta_bs[k];
    y_lat[start[k]] = y0[k];
    
    for (t in (start[k]+1):end[k]) {
      // The exponentially perturbed Gaussian autoregressive process
      R[t] = lam[k] * eta_R[t];
      y_lat[t] = wY[k] * y_lat[t-1] + wT[k]*Treat[t-1] + b_ylat + R[t] + sigma_ylat*eta[t];
    }
  }
}



model {
  sigma_ylat ~ normal(0, 1.5);
  b_ylat ~ normal(0, 1);
  mu_y0 ~ normal(5, 1.5);
  sigma_y0 ~ normal(0, 2);
  y0 ~ normal(mu_y0, sigma_y0);
  
  mu_wY ~ normal(0, 1);
  sigma_wY ~ normal(0, 1.5);
  
  mu_wT ~ normal(0, 1);
  sigma_wT ~ normal(0, 0.5);
  
  sigma_lam ~ normal(0, 1);
  
  c0 ~ normal(0.5, 0.1);
  delta ~ normal(1, 0.5);
  
  wS ~ normal(2, 0.5);
  mu_bs ~ normal(8, 2);
  sigma_bs ~ normal(0, 2);
  sigma_S ~ normal(0, 2.5);
  
  eta ~ std_normal();
  eta_wY ~ std_normal();
  eta_wT ~ std_normal();
  eta_lam ~ std_normal();
  eta_R ~ exponential(1);
  eta_bs ~ std_normal();
  
  for (k in 1:N_pt) {
    (b_ylat + wT[k]) ~ normal(0, 2); // Prior on the constant terms
    
    for (t in z_start[k]:z_end[k]) {
      // SASSAD measurement process
      z_obs[t] ~ normal(y_lat[z_idx_obs[t]]*wS + b_s[k], sigma_S);
    }
  }

  yc_obs ~ ordered_logistic(y_lat[idx_obs], ct); // Bother measurement process
}



generated quantities {
  real y_lat_pred[N_pt];
  real sassad_pred[N_pt];
  
  {
    for (k in 1:N_pt) {
      y_lat_pred[k] = exp_mod_normal_rng(wY[k] * y_lat[end[k]] 
                                            + wT[k] * Treat[end[k]]
                                            + b_ylat, sigma_ylat, 1/lam[k]);
      sassad_pred[k] = normal_rng(y_lat_pred[k]*wS + b_s[k], sigma_S);
    }
  }
}



