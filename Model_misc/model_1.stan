data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing obs

  int<lower = 2, upper = 10> M; // Upper bound of Bother
}



transformed data {
  int yc_obs[N_obs]; // Categorical y_obs
  
  int start[N_pt]; // Index of first observation for each patient
  int end[N_pt]; // Index of last observation for each patient
  
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
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
  real eta[N]; // Non-centred parameterisation for the gaussian process
  
  real<lower = 0, upper = 10> mu_y0; // Population mean of y_lat at t0
  real<lower = 0> sigma_y0; // Population standard deviation of y_lat at t0

  real<lower = 0> c0; // First cutpoint
  real<lower = 0> delta[M-1]; // Difference between cutpoints
}



transformed parameters {
  vector[M] ct; // Cutpoints
  vector[N] y_lat; // Latent score

  ct[1] = c0;
  for (i in 2:M) {
    ct[i] = ct[i-1] + delta[i-1];
  }
  
  for (k in 1:N_pt) {
    y_lat[start[k]] = mu_y0 + sigma_y0*eta[start[k]];
    
    for (t in (start[k]+1):end[k]) {
      y_lat[t] = y_lat[t-1] + sigma_ylat*eta[t];
    }
  }
}



model {
  sigma_ylat ~ normal(0, 2);
  mu_y0 ~ normal(5, 2);
  sigma_y0 ~ normal(0, 2.5);
  
  c0 ~ normal(0.5, 0.1);
  delta ~ normal(1, 0.1);
  
  eta ~ std_normal();

  yc_obs ~ ordered_logistic(y_lat[idx_obs], ct);
}




