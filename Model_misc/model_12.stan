data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing obs
  int<lower = 0> z_N_obs;
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 0> z_t_max[N_pt];
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
  int z_start[N_pt];
  int z_end[N_pt];
  
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
  real eta[N]; // Non-centred parameterisation for the gaussian process
  
  real mu_wY; // Population autocorrelation logit mean
  real<lower = 0> sigma_wY; // Population autocorrelation logit standard dev
  real eta_wY[N_pt]; // Non-centred parameterisation for autocorrelation
  
  real mu_wT; // Population treatment response mean
  real<lower = 0> sigma_wT; // Population treatment response standard dev
  real eta_wT[N_pt]; // Non-centred parameterisation for treatment response
  
  real<lower = 0> c0; // First cutpoint
  real<lower = 0> delta[M-1]; // Difference between cutpoints
  
  real<lower = 0> sigma_lam;
  real<lower = 0> eta_lam[N_pt];
  real<lower = 0> eta_R[N];

  real mu_wS;
  real<lower = 0> sigma_wS;
  real eta_wS[N_pt];
  real<lower = 0> b_S;
  real<lower = 0> sigma_S; 
}



transformed parameters {
  vector[M] ct; // Cutpoints
  vector[N] y_lat; // Latent score

  real wY[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient treatment response
  real lam[N_pt];
  real R[N];
  real wS[N_pt];
  
  ct[1] = c0;
  for (i in 2:M) {
    ct[i] = ct[i-1] + delta[i-1];
  }
  
  
  for (k in 1:N_pt) {
    wY[k] = inv_logit(mu_wY + sigma_wY*eta_wY[k]);
    wT[k] = mu_wT + sigma_wT*eta_wT[k];
    wS[k] = mu_wS + sigma_wS*eta_wS[k];
    lam[k] = sigma_lam*eta_lam[k];
  
    y_lat[start[k]] = y0[k];
    
    for (t in (start[k]+1):end[k]) {
      R[t] = lam[k] * eta_R[t];
      y_lat[t] = wY[k] * y_lat[t-1] + wT[k]*Treat[t-1] + R[t] + sigma_ylat*eta[t];
    }
  }
}


model {
  sigma_ylat ~ normal(0, 1.5);
  
  mu_wY ~ normal(0, 1);
  sigma_wY ~ normal(0, 1.5);
  
  mu_wT ~ normal(0, 1);
  sigma_wT ~ normal(0, 0.5);
  
  sigma_lam ~ normal(0, 1.5);
  
  c0 ~ normal(0.5, 0.5);
  delta ~ normal(1, 1);
  
  mu_wS ~ normal(2, 1);
  sigma_wS ~ normal(0, 1.5);
  b_S ~ normal(0, 10);
  sigma_S ~ normal(0, 5);

  eta ~ std_normal();
  eta_wY ~ std_normal();
  eta_wT ~ std_normal();
  eta_lam ~ std_normal();
  eta_R ~ exponential(1);
  eta_wS ~ std_normal();
  
  for (k in 1:N_pt) {

    for (t in z_start[k]:z_end[k]) {
      z_obs[t] ~ normal(y_lat[z_idx_obs[t]]*wS[k] + b_S, sigma_S);
    }
  }

  yc_obs ~ ordered_logistic(y_lat[idx_obs], ct);
}



