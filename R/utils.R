#' Standardize the n -by- p design matrix X to have
#' column means zero and ||X_j||_2^2 = n for all j
#' @param x design matrix
#' @param center should we set column means equal to zero
standardize <- function(x, center = TRUE){
  n <- nrow(x)
  x <- scale(x, center = center, scale = TRUE)
  col_mean <- attr(x, which = "scaled:center")
  col_sd <- attr(x, which = "scaled:scale")
  return(list(x = x, col_mean = col_mean, col_sd = col_sd))
}

#' Generate sparse linear model and random samples
#'
#' Generate design matrix and response following linear models
#' \eqn{y = X \beta + \epsilon}, where
#' \eqn{\epsilon ~ N(0, \sigma^2)}, and \eqn{X ~ N(0, \Sigma)}.
#'
#' @param n the sample size
#' @param p the number of features
#' @param alpha sparsity, i.e., \eqn{n^\alpha} nonzeros in the true regression coefficient.
#' @param rho pairwise correlation among features
#' @param snr signal to noise ratio, defined as \eqn{\beta^T \Sigma \beta / \sigma^2}
#' @param nsim the number of simulations
#'
#' @return A list object containing: \describe{
#' \item{\code{x}: }{The \code{n} by \code{p} design matrix}
#' \item{\code{y}: }{The \code{n} by \code{nsim} matrix of response vector, each column representing one replication of the simulation}
#' \item{\code{beta}: }{The true regression coefficient vector}
#' \item{\code{sigma}: }{The true error standard deviation}}
#'
#' @export
make_sparse_model <- function(n, p, alpha, rho, snr, nsim) {
  # alpha controls the sparsity, n^\alpha = nnz
  # rho is the correlation between columns of the design matrix X
  # snr: signal-to-noise ratio
  # construct the covariance matrix of the design matrix X
  # x_{ij} ~ N(0, 1)
  # columns of x have correlation rho
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(Sigma)
  # standardize x so that it has colmeans 0 and ||x_j||^2 = n - 1
  x <- scale(x, center = TRUE, scale = TRUE)
  # number of nonzeros
  nnz <- ceiling(n^alpha)
  # random indices of nonzero elements in beta
  ind <- sample(p, nnz)
  # nonzero element values are set equals to a random sample from
  # Laplace(1) (center = 0, scale = 1)
  beta <- rep(0, p)
  rr_unif <- stats::runif(n = 1, min = -0.5, max = 0.5)
  beta[ind] <- -sign(rr_unif) * log(1 - 2 * abs(rr_unif))
  # true sigma
  sigma <- sqrt(as.numeric(crossprod(beta, Sigma %*% beta) / snr))
  # true signal
  mu <- as.numeric(x[, ind] %*% beta[ind])
  # ||beta||_1 / sigma^2
  signal_th <- sum(abs(beta)) / sigma

  y <- mu + sigma * matrix(stats::rnorm(nsim * n), n, nsim)
  return(list(x = x, y = y, beta = beta, sigma = sigma))
}

#' Get the two (theoretical) values of lambdas used in scaled lasso
#' @param n number of observations
#' @param p number of features
getLam_slasso <- function(n, p){
  # use lam_1 = \sqrt{2 \log(p) / n}
  lam_1 <- sqrt(2 * log(p) / n)
  # use lam_2 as described in (10) of Reid, et, al (2016)
  L = 0.1
  Lold = 0
  while (abs(L - Lold) > 1e-3) {
    k = (L^4 + 2 * L^2)
    Lold = L
    L = -stats::qnorm(min(k / p, 0.99))
    L = (L + Lold)/2
  }
  if (p == 1) L = 0.5
  lam_2 = sqrt(2 / n) * L
  return(c(lam_1, lam_2))
}

#' Get the two (theoretical) values of lambdas used in the organic lasso
#' @param x design matrix
getLam_olasso <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  # use lam_1 = \log(p) / n
  lam_1 <- log(p) / n
  # and Monte-Carlo estimate (with B replication) of
  # (||X^T e||_\infty / n)^2, where e ~ N(0, I_n)
  B <- 100
  lam_2 <- abs(crossprod(x, matrix(stats::rnorm(B * n), nrow = n, ncol = B)))
  lam_2 <- mean(apply(lam_2, 2, max) / n)^2
  return(c(lam_1, lam_2))
}
