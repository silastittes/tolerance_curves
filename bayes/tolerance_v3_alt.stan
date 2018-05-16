functions{
  
  real kumara(real x, real a, real b, real c){
    return log(c) + log(a) + log(b) + (a-1)*log(x) + (b-1)*log1m(pow(x,a));
  }
  
  real exp_trans(real x, real lim, int upper){
    if(upper == 1){
      return lim + exp(-x);
    } else {
      return lim - exp(-x);
    }
  }
  
}

data {
  
	int N;
	int numSpp;
	real minx[numSpp];
	real maxx[numSpp];
	int <lower = 1, upper = numSpp> sppint[N];
	real y[N];
	real x[N];
	int zig;
	real a_pr_mu;
	real a_pr_sig;
	real b_pr_mu;
	real b_pr_sig;
	real c_pr_mu;
	real c_pr_sig;
	real d_pr_mu;
	real d_pr_sig;
	real e_pr_mu;
	real e_pr_sig;
}

transformed data {
  
	real min_shp = 1;
	int is_zero[N];

  for (i in 1:N) {
    if (y[i] == 0) {
      is_zero[i] = 1;
    } else {
      is_zero[i] = 0;
    }
  }
}

parameters {

  //species level
  real a_t[numSpp];
  real b_t[numSpp];
  real c_t[numSpp];
  real d_t[numSpp];
  real e_t[numSpp];
  real beta_0[numSpp];
  real beta_1[numSpp];
  real <lower = 0> nu[numSpp];
  
  real mean_a;
  real mean_b;
  real mean_c;
  real mean_d;
  real mean_e;
  real mean_beta_0;
  real mean_beta_1;
  real <lower = 0> mean_nu;
  
  /*
  real <lower = 0> sd_a;
  real <lower = 0> sd_b;
  real <lower = 0> sd_c;
  real <lower = 0> sd_d;
  real <lower = 0> sd_e;
  real <lower = 0> sd_beta_0;
  real <lower = 0> sd_beta_1;
  real <lower = 0> sd_nu;
  */
  
  //gamma shape parameter
  //real nu_t;
  //real nu_t[numSpp];
 

 

}

transformed parameters {
  
  vector<lower=0>[N] mu;
  real <lower=0, upper=1> xs[N];
  vector[N] logit_p_zero;
  vector<lower = min_shp> [numSpp] a;
  vector<lower = min_shp> [numSpp] b;
  vector<lower = 0> [numSpp] c;
  vector[numSpp] d;
  vector[numSpp] e;
  
  //real <lower = 0> nu = exp_trans(nu_t, 0, 1);
  //real <lower = 0> nu[numSpp];
  //nu[i] = exp_trans(nu_t[sppint[i]], 0, 1);

  for(i in 1:numSpp){
    a[i] = exp_trans(a_pr_mu + a_pr_sig*a_t[i], min_shp, 1);
    b[i] = exp_trans(b_pr_mu + b_pr_sig*b_t[i], min_shp, 1);
    c[i] = exp_trans(c_pr_mu + c_pr_sig*c_t[i], 0, 1);
    d[i] = exp_trans(d_pr_mu + d_pr_sig*d_t[i], minx[i], 0);
    e[i] = exp_trans(e_pr_mu + e_pr_sig*e_t[i], maxx[i], 1);
  }
  
  
  for(i in 1:N){
    xs[i]  = (x[i] - d[sppint[i]]) / (e[sppint[i]] - d[sppint[i]]);
    mu[i] = exp(kumara(xs[i], a[sppint[i]], b[sppint[i]], c[sppint[i]]));  
    logit_p_zero[i] = beta_0[sppint[i]] + beta_1[sppint[i]] * mu[i];
  }
}



model {
  
  mean_a ~ normal(0, 1);
  //sd_a ~ normal(0, 1);
  a_t ~ normal(mean_a, 1);

  mean_b ~ normal(0, 1);
  //sd_b ~ normal(0, 1);
  b_t ~ normal(mean_b, 1);

  mean_c ~ normal(0, 1);
  //sd_c ~ normal(0, 1);
  c_t ~ normal(mean_c, 1);

  //mean_d ~ normal(d_pr_mu, 1);
  mean_d ~ normal(0, 1);
  //sd_d ~ normal(0, 1);
  d_t ~ normal(mean_d, 1);

  //mean_e ~ normal(e_pr_mu, 1);
  mean_e ~ normal(0, 1);
  //sd_e ~ normal(0, 1);
  e_t ~ normal(mean_e, 1);

  mean_beta_0 ~ normal(0, 1);
  //sd_beta_0 ~ normal(0, 1);
  beta_0 ~ normal(mean_beta_0, 1);

  mean_beta_1 ~ normal(0, 1);
  //sd_beta_1 ~ normal(0, 1);
  beta_1 ~ normal(mean_beta_1, 1);
  
  mean_nu ~ normal(2, 1);
  //sd_nu ~ normal(0, 1);
  nu ~ gamma(mean_nu, 1);
  //nu_t ~ normal(0, 1);
  

  target += bernoulli_logit_lpmf(is_zero | logit_p_zero);
  
  for (n in 1:N) {
    if (y[n] > 0) {
      if(zig == 0){
        //target += gamma_lpdf( y[n] | nu, (nu * (1 - inv_logit(logit_p_zero[n]))) / mu[n]);  
        target += gamma_lpdf( y[n] | nu[sppint[n]], (nu[sppint[n]] * (1 - inv_logit(logit_p_zero[n]))) / mu[n]);  
      } else {
        //target += gamma_lpdf( y[n] | nu, nu  / mu[n]);
        target += gamma_lpdf( y[n] | nu[sppint[n]], nu[sppint[n]]  / mu[n]);
      }
    }
  }
}


//log like for psis-loo
generated quantities {
  
  vector[N] log_lik;
  for (n in 1:N){
    if(y[n] > 0){
      if(zig == 0){
        //log_lik[n] = gamma_lpdf(y[n] | nu, (nu * (1 - inv_logit(logit_p_zero[n]))) / mu[n]);
        log_lik[n] = gamma_lpdf(y[n] | nu[sppint[n]], (nu[sppint[n]] * (1 - inv_logit(logit_p_zero[n]))) / mu[n]);
      } else {
        //log_lik[n] = gamma_lpdf(y[n] | nu, nu/mu[n]);
        log_lik[n] = gamma_lpdf(y[n] | nu[sppint[n]], nu[sppint[n]]/mu[n]);
      }
    } else {
      log_lik[n] = bernoulli_logit_lpmf(is_zero[n] | logit_p_zero[n]);
    }
  }
}

