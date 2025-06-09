#' Variational Inference for Multivariate Normal Approximation
#'
#' Estimate the mean and covariance matrix of a multivariate normal
#' distribution that approximates an unknown density using variational inference.
#'
#' @param samples A numeric matrix of samples (rows are observations).
#' @param num_steps Number of optimization steps. Default is 1000.
#' @return A list with elements \code{mu} (estimated mean vector) and \code{Sigma} (estimated covariance matrix).
#' @export
variational_inference <- function(samples, num_steps = 1000) {
  n_samples <- nrow(samples)
  dim <- ncol(samples)
  mu_init <- colMeans(samples) + rnorm(dim)
  L_init <- matrix(0, nrow = dim, ncol = dim)
  diag(L_init) <- log(rep(1, dim))
  L_init[lower.tri(L_init)] <- rnorm(dim * (dim - 1) / 2)
  params_init <- c(mu_init, as.vector(L_init))
  reconstruct_L <- function(params) {
    L <- matrix(0, nrow = dim, ncol = dim)
    idx <- 1
    for (i in 1:dim) {
      for (j in 1:i) {
        L[i, j] <- params[idx]
        idx <- idx + 1
      }
    }
    diag(L) <- exp(diag(L))
    L
  }
  neg_log_likelihood <- function(params) {
    mu <- params[1:dim]
    L_params <- params[-(1:dim)]
    L <- reconstruct_L(L_params)
    Sigma <- L %*% t(L)
    log_det_Sigma <- 2 * sum(log(diag(L)))
    Sigma_inv <- solve(Sigma)
    diffs <- t(t(samples) - mu)
    quad_form <- rowSums((diffs %*% Sigma_inv) * diffs)
    log_density <- -0.5 * (dim * log(2 * pi) + log_det_Sigma + quad_form)
    -mean(log_density)
  }
  opt <- optim(par = params_init, fn = neg_log_likelihood, method = "BFGS",
               control = list(maxit = num_steps))
  final_mu <- opt$par[1:dim]
  final_L_params <- opt$par[-(1:dim)]
  final_L <- reconstruct_L(final_L_params)
  final_Sigma <- final_L %*% t(final_L)
  list(mu = final_mu, Sigma = final_Sigma)
}