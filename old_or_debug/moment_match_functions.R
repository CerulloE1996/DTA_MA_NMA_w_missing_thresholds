# create a named list of draws for use with rstan methods
.rstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function(pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}

# extract original posterior draws
post_draws_stanfit <- function(x, ...) {
  as.matrix(x)
}

# compute a matrix of log-likelihood values for the ith observation
# matrix contains information about the number of MCMC chains
log_lik_i_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[, , i]
}

# transform parameters to the unconstrained space
unconstrain_pars_stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob, object = x,
        adjust_transform = TRUE, gradient = FALSE)
}

# compute log_lik values based on the unconstrained parameters
log_lik_i_upars_stanfit <- function(x, upars, i, parameter_name = "log_lik",
                                    ...) {
  S <- nrow(upars)
  out <- numeric(S)
  for (s in seq_len(S)) {
    out[s] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]][i]
  }
  out
}

