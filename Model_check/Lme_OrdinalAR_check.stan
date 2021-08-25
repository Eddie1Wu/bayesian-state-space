data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing obs
  int<lower = 0> z_N_obs; // Number of non-missing SASSAD obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 0> z_t_max[N_pt]; // Time-series length of SASSAD
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing obs
  int<lower = 0, upper = N> z_idx_obs[z_N_obs];
  int<lower = 0> z_obs[z_N_obs];
  int<lower = 0, upper = 10> y0[N_pt];
  real<lower = 0, upper = 1> Treat[N]; // Daily treatment usage
  
  int<lower = 2, upper = 10> M; // Upper bound of Bother
}



transformed data {
  int yc_obs[N_obs]; // Categorical y_obs
  
  int start[N_pt]; // Index of first observation for each patient
  int end[N_pt]; // Index of last observation for each patient
  int z_start[N_pt]; // Index of first SASSAD obs for each patient
  int z_end[N_pt]; // Index of last SASSAD obs for each patient
  
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
  }
  
  for (k in 1:N_pt) {
    if (k == 1) {
      start[k] = 1;
      z_start[k] = 1;
    } else {
      start[k] = end[k-1] + 1;
      z_start[k] = z_end[k-1] + 1;
    }
    end[k] = start[k] - 1 + t_max[k];
    z_end[k] = z_start[k] -1 + z_t_max[k];
  }
}



parameters {
  real<lower = 0> sigma_ylat; // Standard deviation of y_lat
  real b_ylat; // Intercept of y_lat
  real eta[N]; // Non-centred parameterisation for the gaussian process
  
  real mu_wY; // Population autocorrelation logit mean
  real<lower = 0> sigma_wY; // Population autocorrelation logit standard dev
  real eta_wY[N_pt]; // Non-centred parameterisation for autocorrelation
  
  real mu_wT; // Population treatment response mean
  real<lower = 0> sigma_wT; // Population treatment response standard dev
  real eta_wT[N_pt]; // Non-centred parameterisation for treatment response
  
  real<lower = 0> c0; // First cutpoint
  real<lower = 0> delta[M-1]; // Difference between cutpoints
  
  real wS;
  real mu_bs;
  real<lower = 0> sigma_bs;
  real eta_bs[N_pt];
  real<lower = 0> sigma_S; 
}



transformed parameters {
  vector[M] ct; // Cutpoints
  vector[N] y_lat; // Latent score

  real wY[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient treatment response
  real b_s[N_pt];
  
  ct[1] = c0;
  for (i in 2:M) {
    ct[i] = ct[i-1] + delta[i-1];
  }
  
  
  for (k in 1:N_pt) {
    wY[k] = inv_logit(mu_wY + sigma_wY*eta_wY[k]);
    wT[k] = mu_wT + sigma_wT*eta_wT[k];
    b_s[k] = mu_bs + sigma_bs*eta_bs[k];
  
    y_lat[start[k]] = y0[k];
    
    for (t in (start[k]+1):end[k]) {
      // The autoregressive process
      y_lat[t] = wY[k] * y_lat[t-1] + wT[k]*Treat[t-1] + b_ylat + sigma_ylat*eta[t];
    }
  }
}



model {
  b_ylat ~ normal(1, 1);
  sigma_ylat ~ normal(0, 1.5);
  
  mu_wY ~ normal(0, 1);
  sigma_wY ~ normal(0, 1.5);
  
  mu_wT ~ normal(0, 1);
  sigma_wT ~ normal(0, 0.5);
  
  c0 ~ normal(0.5, 0.5);
  delta ~ normal(1, 1);
  
  wS ~ normal(2, 1);
  mu_bs ~ normal(10, 3);
  sigma_bs ~ normal(0, 5);
  sigma_S ~ normal(0, 3.5);

  eta ~ std_normal();
  eta_wY ~ std_normal();
  eta_wT ~ std_normal();
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
  int y_pred[N];
  real sassad_pred[N_pt];
  
  for (n in 1:N) {
    y_pred[n] = ordered_logistic_rng(y_lat[n], ct) - 1;
  }
  
  for (k in 1:N_pt) {
    sassad_pred[k] = normal_rng(y_lat[start[k]]*wS + b_s[k], sigma_S);
  }
}


