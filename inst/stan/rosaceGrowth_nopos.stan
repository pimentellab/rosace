data {
  int<lower=0> T; // # of time points
  int<lower=0> V; // # of variants
  int<lower=0> M; // # of mean group
  array[V] int vMAPm;
  vector[T] t; // time
  array[V] vector[T] m; // normalized count
}
parameters {
  vector[V] beta;
  vector<lower=0>[M] epsilon2;
  vector[V] b; // intercept
}
model {
  beta ~ normal(0, 1);
  epsilon2 ~ inv_gamma(1, 1);
  b ~ normal(0, 0.25);
  for (v in 1:V) {
    m[v] ~ normal(b[v] + beta[v] * t, sqrt(epsilon2[vMAPm[v]]));
  }
}
