data {
  int<lower=0> T; // # of time points
  int<lower=0> V; // # of variants
  int<lower=0> P; // # of positions
  int<lower=0> M; // # of mean group
  array[V] int vMAPp;
  array[V] int vMAPm;
  vector[T] t; // time
  array[V] vector[T] m; // normalized count
}
parameters {
  vector[P] phi; // slope per position
  vector<lower=0>[P] sigma2; 
  vector[V] eta2; // std_normal for beta
  vector<lower=0>[M] epsilon2;
  vector[V] b; // intercept
}
transformed parameters {
  vector[V] beta; // slope per variants
  for (v in 1:V) {
    beta[v] = phi[vMAPp[v]] + eta2[v] * sqrt(sigma2[vMAPp[v]]);
  }
}
model {
  phi ~ normal(0, 1);
  sigma2 ~ inv_gamma(1, 1);
  eta2 ~ normal(0, 1);
  epsilon2 ~ inv_gamma(1, 1);
  b ~ normal(0, 0.25);
  for (v in 1:V) {
    m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
  }
}