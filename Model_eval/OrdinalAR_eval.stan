data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing obs
  int<lower = 0, upper = 10> y0[N_pt];
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
  
  int test_start[N_pt]; 
  int test_end[N_pt];
  
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
  }
  
  for (k in 1:N_pt) {
    if (k == 1) {
      start[k] = 1;
      test_start[k] = 1;
    } else {
      start[k] = end[k-1] + 1;
      test_start[k] = test_end[k-1] + 1;
    }
    end[k] = start[k] - 1 + t_max[k];
    test_end[k] = test_start[k] - 1 + test_t_max[k];
  }
}



parameters {
  real<lower = 0> delta[M-1]; // Difference between cutpoints

  real<lower = 0> sigma_ylat; // Standard deviation of y_lat
  real b_ylat; // Intercept of y_lat
  
  real mu_wY; // Population autocorrelation logit mean
  real<lower = 0> sigma_wY; // Population autocorrelation logit standard dev
  real eta_wY[N_pt]; // Non-centred parameterisation for autocorrelation
  
  real mu_wT; // Population treatment response mean
  real<lower = 0> sigma_wT; // Population treatment response standard dev
  real eta_wT[N_pt]; // Non-centred parameterisation for treatment response
  
  real eta[N]; // Non-centred parameterisation for the gaussian process
  real<lower = 0> c0; // First cutpoint
}



transformed parameters {
  vector[M] ct; // Cutpoints
  vector[N] y_lat; // Latent score
  
  real wY[N_pt]; // Patient autocorrelation
  real wT[N_pt]; // Patient treatment response
  
  ct[1] = c0;
  for (i in 2:M) {
    ct[i] = ct[i-1] + delta[i-1];
  }
  
  
  for (k in 1:N_pt) {
    wY[k] = inv_logit(mu_wY + sigma_wY*eta_wY[k]);
    wT[k] = mu_wT + sigma_wT*eta_wT[k];
    y_lat[start[k]] = y0[k];
    
    for (t in (start[k]+1):end[k]) { 
      // The autoregressive process
      y_lat[t] = wY[k] * y_lat[t-1] + wT[k]*Treat[t-1] + b_ylat + sigma_ylat*eta[t];
    }
  }
}


model {
  c0 ~ normal(0.5, 0.5);
  eta ~ std_normal();
  delta ~ normal(1, 1);
  
  eta_wY ~ std_normal();
  eta_wT ~ std_normal();
  
  b_ylat ~ normal(1, 1);
  sigma_ylat ~ normal(0, 1.5);
  mu_wY ~ normal(0, 1);
  sigma_wY ~ normal(0, 1.5);
  mu_wT ~ normal(0, 1);
  sigma_wT ~ normal(0, 0.5);
  
  for (k in 1:N_pt) {
    (b_ylat + wT[k]) ~ normal(0, 2); // Prior on the constant terms
  }

  yc_obs ~ ordered_logistic(y_lat[idx_obs], ct); // Bother measurement process
}



generated quantities {
  real y_lat_pred[N_test];
  real lpd[N_test];
  int y_pred[N_test];
  real cdf[N_test, M+1];
  
  {
    for (k in 1:N_pt) {
      y_lat_pred[test_start[k]] = normal_rng(wY[k] * y_lat[end[k]] 
                                            + wT[k] * Treat[end[k]]
                                            + b_ylat, sigma_ylat);
      
      for (t in (test_start[k]+1):test_end[k]) {
        y_lat_pred[t] = normal_rng(wY[k] * y_lat_pred[t-1] + b_ylat, sigma_ylat);
      }
    }
  }
  
  for (n in 1:N_test) {
    if (y_test[n] == 99) {
      lpd[n] = 99;
      for (i in 1:M+1) {
        cdf[n, i] = 99;
      }
      y_pred[n] = 99;
      
    } else {
      lpd[n] = ordered_logistic_lpmf(y_test[n] + 1 | y_lat_pred[n], ct);
      y_pred[n] = ordered_logistic_rng(y_lat_pred[n], ct) - 1;
      
      cdf[n, 1] = exp(ordered_logistic_lpmf( 1 | y_lat_pred[n], ct));
      for (i in 1:M) {
        cdf[n, i+1] = cdf[n, i] + exp(ordered_logistic_lpmf( i+1 | y_lat_pred[n], ct));
      }
    }
  }
}


