# =============================================================================
# miCCI / 40_cci_bayesian.R
# S4 - Bayesian posterior CCI.
#
# For each truncated prefix p_i with subcode probability vector q_i, place a
# Dirichlet(alpha_0 * q_i + 0.5) prior on the unknown subcode probabilities
# theta_i. Draw n_draws samples; for each draw, sample one full-length code
# from theta_i, recompute the gold CCI on the imputed code set, and return
# the posterior median across draws.
#
# alpha_0 governs prior strength. With alpha_0 large the posterior collapses
# onto q_i and S4 behaves like S3. With alpha_0 small the Dirichlet has high
# variance and S4 is more dispersed than S3. We expose alpha_0 to the
# top-level pipeline so it can be swept.
# =============================================================================

#' Single-encounter S4 (defers to batch on length-1 input).
#' @export
cci_bayesian <- function(icd_anon, quan_map, cache,
                         age = NULL, n_draws = 25L, alpha_0 = 10, seed = 42L,
                         preserve_multiplicity = TRUE,
                         aggregator = c("median", "mean")) {
  aggregator <- match.arg(aggregator)
  codes <- truncate_icd(icd_anon, preserve_multiplicity = preserve_multiplicity)
  one <- data.table(diagnosen = paste(codes, collapse = "|"))
  v <- cci_bayesian_batch(one, quan_map, cache,
                          n_draws = n_draws, alpha_0 = alpha_0, seed = seed,
                          age_idx = if (is.null(age)) NULL else age_to_bin_index(age),
                          preserve_multiplicity = preserve_multiplicity,
                          aggregator = aggregator)
  list(posterior_median = unname(v[1L]))
}

#' Vectorised S4 batch.
#'
#' @inheritParams cci_mi_batch
#' @param n_draws number of Dirichlet draws.
#' @param alpha_0 prior pseudo-count multiplier.
#' @param return_group_count If TRUE, also return per-(idx, gk) trigger
#'   probability across draws (used by the QA mass-conservation check).
#'
#' @return numeric vector of posterior medians, length `nrow(dt)`
#'   (or list with that vector plus a `group_prob` table).
#' @export
cci_bayesian_batch <- function(dt, quan_map, cache,
                               n_draws = 25L, alpha_0 = 10, seed = 42L,
                               return_group_count = FALSE,
                               age_idx = NULL,
                               preserve_multiplicity = TRUE,
                               aggregator = c("median", "mean"),
                               return_draws = FALSE) {
  aggregator <- match.arg(aggregator)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    p <- substr(toupper(gsub("[^A-Z0-9]", "", x)), 1L, 3L)
    p <- p[nchar(p) >= 3L]
    if (preserve_multiplicity) p else unique(p)
  })

  ak_vec <- if (is.null(age_idx)) rep(0L, n) else {
    a <- as.integer(age_idx); a[is.na(a)] <- 0L; a
  }
  key_sets <- lapply(seq_len(n), function(i) {
    if (length(code_sets[[i]]) == 0L) character(0L)
    else paste0(code_sets[[i]], "@", ak_vec[i])
  })

  # The Dirichlet concentration is built from whichever prior the pool
  # resolved to (age-conditioned where available, marginal otherwise), so the
  # smoothing and the prior strength keep the same meaning in both modes.
  all_keys <- unique(unlist(key_sets, use.names = FALSE))
  dir_params <- list()
  n_degenerate <- 0L; n_agefall <- 0L
  for (key in all_keys) {
    parts <- strsplit(key, "@", fixed = TRUE)[[1L]]
    pref  <- parts[1L]; ak <- as.integer(parts[2L])
    pool  <- prefix_pool(pref, cache, age_idx = if (ak == 0L) NULL else ak)
    if (nrow(pool) == 0L) pool <- data.table(code_nodot = pref, prob = 1)
    if (identical(attr(pool, "status", exact = TRUE), "degenerate"))
      n_degenerate <- n_degenerate + 1L
    if (isTRUE(attr(pool, "age_fallback", exact = TRUE)))
      n_agefall <- n_agefall + 1L
    dir_params[[key]] <- list(alpha = alpha_0 * pool$prob + 0.5,
                              names = pool$code_nodot)
  }

  message(sprintf(paste0("S4: n=%d encounters, %d prefix-by-age pools,",
                         " n_draws=%d, alpha_0=%g, agg=%s",
                         " (%d degenerate, %d age fallbacks)"),
                  n, length(all_keys), n_draws, alpha_0, aggregator,
                  n_degenerate, n_agefall))

  cci_mat <- matrix(0, nrow = n, ncol = n_draws)
  group_count_dt <- if (return_group_count)
    data.table(idx = integer(0L), gk = character(0L), n_active = integer(0L))
  else NULL

  .with_local_seed(seed, {
    for (d in seq_len(n_draws)) {
      imputed_dx <- vapply(seq_len(n), function(i) {
        keys <- key_sets[[i]]
        if (length(keys) == 0L) return("")
        drawn <- vapply(keys, function(key) {
          par   <- dir_params[[key]]
          theta <- stats::rgamma(length(par$alpha), shape = par$alpha, rate = 1)
          theta <- theta / sum(theta)
          add_dot4(par$names[sample.int(length(theta), 1L, prob = theta)])
        }, character(1L))
        paste(drawn, collapse = "|")
      }, character(1L))

      imp_dt <- data.table(diagnosen = imputed_dx)
      if (return_group_count) {
        ga <- .gold_active_long_internal(imp_dt, quan_map, pl, dl)
        if (nrow(ga) > 0L) {
          ga[, n_active := 1L]
          group_count_dt <- rbindlist(list(group_count_dt, ga),
                                      use.names = TRUE, fill = TRUE)
        }
      }
      cci_mat[, d] <- cci_gold_batch(imp_dt, quan_map, pl, dl)
    }
  })

  # Draw d consumes the same random numbers whatever n_draws is, so the first k
  # columns of this matrix ARE the k-draw run, and both aggregators can be read
  # from one set of draws. Returned raw for the sensitivity analysis.
  if (return_draws) return(cci_mat)

  # With an odd number of draws over an integer-valued integrand the median is
  # itself an integer for every encounter, so S4 is discrete and sits on the
  # same support as the reference score. The mean aggregator is exposed for the
  # sensitivity analysis and is real-valued.
  out <- if (aggregator == "median")
    apply(cci_mat, 1L, stats::median) else rowMeans(cci_mat)
  if (return_group_count) {
    gp <- group_count_dt[, .(p = sum(n_active) / n_draws), by = .(idx, gk)]
    return(list(posterior_median = out, group_prob = gp))
  }
  out
}
