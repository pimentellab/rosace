#' @importFrom cmdstanr write_stan_file
#'
NULL

WriteStanModel <- function(type) {

  if (type == "growth_pos") {
    stan_program <- "
    functions{
      vector repeat_vector(vector input, int reps) {
        return(to_vector(rep_matrix(input, reps)));
      }
    }
    data {
      int<lower=0> T; // # of time points
      int<lower=0> V; // # of variants
      int<lower=0> P; // # of positions
      int<lower=0> M; // # of mean group
      array[V] int vMAPp;
      array[V] int vMAPm;
      row_vector[T] t; // time
      matrix[V, T] m; // normalized count
    }
    transformed data {
      vector[V * T] m_vec = to_vector(m); // Vectorized, column-major order
      vector[V * T] t_vec = to_vector(rep_matrix(t, V));
    }
    parameters {
      vector[P] phi; // slope per position
      vector<lower=0>[P] sigma2; 
      vector[V] eta2; // std_normal for beta
      vector<lower=0>[M] epsilon2;
      vector[V] b; // intercept
    }
    transformed parameters {
      vector<lower=0>[V] sigma2_v = sigma2[vMAPp];
      vector[V] phi_v = phi[vMAPp];
      vector[V] beta = phi_v + eta2 .* sqrt(sigma2_v); // slope per variants
      vector<lower=0>[V] epsilon2_v = epsilon2[vMAPm];
    }
    model {
      phi ~ normal(0, 1);
      sigma2 ~ inv_gamma(1, 1);
      eta2 ~ normal(0, 1);
      epsilon2 ~ inv_gamma(1, 1);
      b ~ normal(0, 0.25);
      m_vec ~ normal(repeat_vector(b, T) + repeat_vector(beta, T) .* t_vec, repeat_vector(sqrt(epsilon2_v), T));
    }
    "
  } else if (type == "growth_nopos") {
    stan_program <- "
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
    "
  } else {
    stop("Invalid experiment type. Currently support 'growth' and 'binding'.")
  }

  file_pedantic <- cmdstanr::write_stan_file(code = stan_program,
                       dir = getOption("cmdstanr_write_stan_file_dir", tempdir()))

  return(file_pedantic)
}


