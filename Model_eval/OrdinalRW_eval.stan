data {
  int<lower = 0> N; // Total number of obs
  int<lower = 0> N_obs; // Number of non-missing obs
  int<lower = 0> N_pt; // Number of patients
  
  int<lower = 0> t_max[N_pt]; // Time-series length for each patient
  int<lower = 1, upper = N> idx_obs[N_obs]; // Vector of indices of non-missing obs
  int<lower = 0, upper = 10> y_obs[N_obs]; // Non-missing obs
  
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
  real<lower = 0> sigma_ylat; // Standard deviation of y_lat
  real eta[N]; // Non-centred parameterisation for the gaussian process
  
  real<lower = 2, upper = 8> mu_y0; // Population mean of y_lat at t0
  real<lower = 0> sigma_y0; // Population standard deviation of y_lat at t0
  real<lower = 0, upper = 10> y0[N_pt]; // Initial latent score for each patient

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
    y_lat[start[k]] = y0[k];
    
    for (t in (start[k]+1):end[k]) {
      y_lat[t] = y_lat[t-1] + sigma_ylat*eta[t];
    }
  }
}



model {
  sigma_ylat ~ normal(0, 1.5);
  mu_y0 ~ normal(5, 1.5);
  sigma_y0 ~ normal(0, 1);
  y0 ~ normal(mu_y0, sigma_y0);
  
  c0 ~ normal(0.5, 0.1);
  delta ~ normal(1, 1);
  
  eta ~ std_normal();

  yc_obs ~ ordered_logistic(y_lat[idx_obs], ct);
}



generated quantities {
  real y_lat_pred[N_test];
  real lpd[N_test];
  int y_pred[N_test];
  real cdf[N_test, M+1];

  {
    for (k in 1:N_pt) {
      y_lat_pred[test_start[k]] = normal_rng(y_lat[end[k]], sigma_ylat);
      
      for (t in (test_start[k]+1):test_end[k]) {
        y_lat_pred[t] = normal_rng(y_lat_pred[t-1], sigma_ylat);
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


