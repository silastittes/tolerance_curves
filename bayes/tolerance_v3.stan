functions{
  real kimura(real x, real a, real b, real c){
    //return c*a*b*pow(x,a-1)*pow(1-pow(x,a),b-1);
    return log(c) + log(a) + log(b) + (a-1)*log(x) + (b-1)*log1m(pow(x,a));
  }
}

data {
	int N;
	int numSpp;
	int<lower = 1, upper = numSpp> sppint[N];
	real y[N];
	real x[N];
}

transformed data {
  real minx;
	real maxx;
	int is_zero[N];

  minx = min(x);
  maxx = max(x);
  for (i in 1:N) {
    if (y[i] == 0) {
      is_zero[i] = 1;
    } else {
      is_zero[i] = 0;
    }
  }
}

parameters {
   
  #species level
  real <lower = 2> a[numSpp];
  real <lower = 2> b[numSpp];
  real <lower = 0> c[numSpp];
  real <lower = 0> mean_c;
  real <lower = 0> var_c;

  real <upper = minx> d[numSpp];
  real <lower = maxx> e[numSpp];

  # gamma shape parameter
  real <lower = 0> nu;
  
  # zero-inflation coefficients
  real beta_0[numSpp];
  real <upper = 0> beta_1[numSpp];
}

transformed parameters {
  vector<lower=0>[N] mu;
  real <lower=0, upper=1> xs[N];
  vector[N] logit_p_zero;
  real e1[numSpp];

  for (i in 1:numSpp) {
    e1[i] = e[i] - d[i];
  }

  for(i in 1:N){
    xs[i]  = (x[i] - d[sppint[i]]) / e1[sppint[i]];
    mu[i] = exp(kimura(xs[i], a[sppint[i]], b[sppint[i]], c[sppint[i]]));
    logit_p_zero[i] = beta_0[sppint[i]] + beta_1[sppint[i]] * mu[i];
  }

}

model {
  
  a ~ normal(4, 1);

  b ~ normal(3, 1);
  
  c ~ normal(mean_c, var_c);
  mean_c ~ normal(0, 10);
  var_c ~ normal(0, 10);
  
  d ~ normal(minx, 2);
  
  e ~ normal(maxx, 2);

  nu ~ gamma(20, 0.2);
  beta_0 ~ normal(0, 2);
  beta_1 ~ normal(0, 2);

  // likelihood
  is_zero ~ bernoulli_logit(logit_p_zero);
  for (i in 1:N) {
    if (y[i] > 0) {
      // y[i] ~ gamma(nu, nu * / mu[i]);
      y[i] ~ gamma(nu, (nu * (1 - inv_logit(logit_p_zero[i]))) / mu[i]);
      }
    }
}
