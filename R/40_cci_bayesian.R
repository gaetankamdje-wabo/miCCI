# =============================================================================
# miCCI — 40_cci_bayesian.R
# S4: Bayesian CCI — VECTORISED: Dirichlet draws in bulk + vectorised CCI
# Same architecture as S3: draw n_draws datasets, batch-CCI each.
# =============================================================================

#' @export
cci_bayesian <- function(icd_anon, quan_map, cache,
                         age = NULL, n_draws = 25L, alpha_0 = 10, seed = 42L) {
  set.seed(seed)
  codes <- unique(icd3(icd_anon))
  if (length(codes) == 0) return(list(posterior_median = 0))
  posteriors <- lapply(codes, function(code) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) return(NULL)
    p <- pc$child_freqs[prob > 0]
    if (nrow(p) == 0) return(NULL)
    list(alpha = alpha_0 * p$prob + 0.5, names = p$code_nodot)
  })
  posteriors <- Filter(Negate(is.null), posteriors)
  if (length(posteriors) == 0) return(list(posterior_median = 0))
  draws <- integer(n_draws)
  for (d in seq_len(n_draws)) {
    imputed <- character(0)
    for (post in posteriors) {
      theta <- rgamma(length(post$alpha), shape = post$alpha, rate = 1)
      theta <- theta / sum(theta)
      imputed <- c(imputed, add_dot4(post$names[sample.int(length(theta), 1L, prob = theta)]))
    }
    draws[d] <- cci_gold(imputed, quan_map)$cci
  }
  list(posterior_median = round(median(draws), 2))
}

#' VECTORISED S4: Dirichlet bulk draws + vectorised CCI
cci_bayesian_batch <- function(dt, quan_map, cache,
                               n_draws = 25L, alpha_0 = 10, seed = 42L) {
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n <- nrow(dt)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1, 3))
  })

  # Build per-prefix Dirichlet parameters (from cache, ONCE)
  all_prefs <- unique(unlist(code_sets))
  dir_params <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      dir_params[[pref]] <- list(alpha = 1, names = pref)
    } else {
      p <- pc$child_freqs[prob > 0]
      if (nrow(p) == 0) {
        dir_params[[pref]] <- list(alpha = 1, names = pref)
      } else {
        dir_params[[pref]] <- list(alpha = alpha_0 * p$prob + 0.5, names = p$code_nodot)
      }
    }
  }

  message(sprintf("  S4: %d encounters, %d unique prefixes, n_draws=%d", n, length(all_prefs), n_draws))

  # Draw n_draws imputed datasets and compute CCI in bulk
  cci_mat <- matrix(0L, nrow = n, ncol = n_draws)
  for (d in seq_len(n_draws)) {
    imputed_dx <- vapply(seq_len(n), function(i) {
      prefs <- code_sets[[i]]
      drawn <- vapply(prefs, function(pref) {
        par <- dir_params[[pref]]
        theta <- rgamma(length(par$alpha), shape = par$alpha, rate = 1)
        theta <- theta / sum(theta)
        add_dot4(par$names[sample.int(length(theta), 1L, prob = theta)])
      }, character(1))
      paste(drawn, collapse = "|")
    }, character(1))

    imp_dt <- data.table(diagnosen = imputed_dx)
    cci_mat[, d] <- cci_gold_batch(imp_dt, quan_map, pl, dl)
  }

  # Posterior median per encounter
  apply(cci_mat, 1, function(x) round(median(x), 2))
}
