functions{
  real kimura(real x, real a, real b, real c){
    //return c*a*b*pow(x,a-1)*pow(1-pow(x,a),b-1);
    return log(c) + log(a) + log(b) + (a-1) * log(x) + (b-1) * log1m(pow(x,a));
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
  minx = min(x);
  maxx = max(x);
}

parameters {
   
  #species level
  real <lower = 2> a[numSpp];
  real <lower = 2> b[numSpp];
  real <lower = 1> c[numSpp];

  real <upper = minx> d[numSpp];
  real <lower = maxx> e[numSpp];

  real <lower = 0> sigma;
}

transformed parameters {
  vector<lower=0>[N] mu;
  real <lower=0, upper=1> xs[N];
  real e1[numSpp];

  for (i in 1:numSpp) {
    e1[i] = e[i] - d[i];
  }

  for(i in 1:N){
    xs[i]  = (x[i] - d[sppint[i]]) / e1[sppint[i]];
    mu[i] = exp(kimura(xs[i], a[sppint[i]], b[sppint[i]], c[sppint[i]]));
  }

}

model {
  
  a ~ normal(2, 2);

  b ~ normal(2, 2);
  
  c ~ normal(2, 2);
  
  d ~ normal(minx, 2);
  
  e ~ normal(maxx, 2);

  sigma ~ normal(0, 2);

  // likelihood

      y ~ normal(mu, sigma);
}
