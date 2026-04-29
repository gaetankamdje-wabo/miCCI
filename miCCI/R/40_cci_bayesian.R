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
                         age = NULL, n_draws = 25L, alpha_0 = 10, seed = 42L) {
  one <- data.table(diagnosen = paste(unique(icd3(icd_anon)), collapse = "|"))
  v <- cci_bayesian_batch(one, quan_map, cache,
                          n_draws = n_draws, alpha_0 = alpha_0, seed = seed)
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
                               return_group_count = FALSE) {
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1L, 3L))
  })
  code_sets <- lapply(code_sets, function(x) x[nchar(x) >= 3L])

  all_prefs <- unique(unlist(code_sets, use.names = FALSE))
  dir_params <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      dir_params[[pref]] <- list(alpha = 1, names = pref)
    } else {
      p <- pc$child_freqs[prob > 0]
      if (nrow(p) == 0L) {
        dir_params[[pref]] <- list(alpha = 1, names = pref)
      } else {
        dir_params[[pref]] <- list(alpha = alpha_0 * p$prob + 0.5,
                                   names = p$code_nodot)
      }
    }
  }

  message(sprintf("S4: n=%d encounters, %d unique prefixes, n_draws=%d, alpha_0=%g",
                  n, length(all_prefs), n_draws, alpha_0))

  cci_mat <- matrix(0, nrow = n, ncol = n_draws)
  group_count_dt <- if (return_group_count)
    data.table(idx = integer(0L), gk = character(0L), n_active = integer(0L))
  else NULL

  .with_local_seed(seed, {
    for (d in seq_len(n_draws)) {
      imputed_dx <- vapply(seq_len(n), function(i) {
        prefs <- code_sets[[i]]
        if (length(prefs) == 0L) return("")
        drawn <- vapply(prefs, function(pref) {
          par   <- dir_params[[pref]]
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

  out <- apply(cci_mat, 1L, stats::median)
  if (return_group_count) {
    gp <- group_count_dt[, .(p = sum(n_active) / n_draws), by = .(idx, gk)]
    return(list(posterior_median = out, group_prob = gp))
  }
  out
}
