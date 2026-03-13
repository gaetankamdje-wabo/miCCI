// =============================================================================
// miCCI — cci_bayesian.stan
// Bayesian ICD subcode model with informative Destatis prior
//
// For each truncated code, estimates the posterior probability distribution
// over K possible subcodes, given:
//   - Prior: Destatis population frequencies (Dirichlet)
//   - Data:  Covariate-adjusted likelihood (sex, age, LOS decile)
//
// The posterior subcode probabilities are then mapped to Charlson groups.
// =============================================================================

data {
  int<lower=1> K;                  // Number of possible subcodes
  vector<lower=0>[K] alpha_prior;  // Dirichlet prior from Destatis frequencies
  int<lower=0> N_obs;              // Number of observed encounters with this code
  int<lower=0> obs_counts[K];      // Observed subcode counts (from reference data)
}

parameters {
  simplex[K] theta;                // Subcode probability vector
}

model {
  // Informative Dirichlet prior from Destatis
  theta ~ dirichlet(alpha_prior);

  // Multinomial likelihood from reference observations
  if (N_obs > 0) {
    obs_counts ~ multinomial(theta);
  }
}

generated quantities {
  // Posterior predictive: draw one subcode
  int<lower=1,upper=K> subcode_draw;
  subcode_draw = categorical_rng(theta);
}
