#' Fit a linear model with organic lasso
#'
#' Calculate a solution path of the organic lasso estimate (of error standard deviation) with a list of tuning parameter values. In particular, this function solves the squared-lasso problems and returns the objective function values as estimates of the error variance:
#' \eqn{\tilde{\sigma}^2_{\lambda} = \min_{\beta} ||y - X \beta||_2^2 / n + 2 \lambda ||\beta||_1^2.}
#'
#' This package also includes the outputs of the naive and the degree-of-freedom adjusted estimates, in analogy to \code{\link{nlasso_path}}.
#'
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{flmin}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param flmin The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2}.
#' @param thresh Threshold value for underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.
#' @param intercept Indicator of whether intercept should be fitted. Default to be \code{FALSE}.
#' @return A list object containing: \describe{
#' \item{\code{n} and \code{p}: }{The dimension of the problem.}
#' \item{\code{lambda}: }{The path of tuning parameter used.}
#' \item{\code{a0}: }{Estimate of intercept. A vector of length \code{nlam}.}
#' \item{\code{beta}: }{Matrix of estimates of the regression coefficients, in the original scale. The matrix is of size \code{p} by \code{nlam}. The \code{j}-th column represents the estimate of coefficient corresponding to the \code{j}-th tuning parameter in \code{lambda}.}
#' \item{\code{sig_obj_path}: }{Organic lasso estimates of the error standard deviation. A vector of length \code{nlam}.}
#' \item{\code{sig_naive}: }{Naive estimate of the error standard deviation based on the squared-lasso regression. A vector of length \code{nlam}.}
#' \item{\code{sig_df}: }{Degree-of-freedom adjusted estimate of the error standard deviation, based on the squared-lasso regression. A vector of length \code{nlam}. }
#' \item{\code{type}: }{whether the output is of a natural or an organic lasso.}}
#' @examples
#' set.seed(123)
#' sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
#' ol_path <- olasso_path(x = sim$x, y = sim$y[, 1])
#'
#' @import Matrix
#' @export
#'
#' @useDynLib natural, .registration = TRUE
#' @seealso \code{\link{olasso}, \link{olasso_cv}}
olasso_path <- function(x, y, lambda = NULL,
                   nlam = 100, flmin = 1e-2,
                   thresh = 1e-8, intercept = TRUE) {
  # This solves the lasso problem
  # 0.5 / n ||y - X\beta||^2 + \lambda ||beta||_1^2
  n <- nrow(x)
  p <- ncol(x)
  # first standardize the design matrix to have column means zero
  # and column sds 1
  stan <- standardize(x, center = intercept)
  x <- stan$x
  col_mean <- stan$col_mean
  col_sd <- stan$col_sd

  y_mean <- mean(y)
  if (intercept)
    y <- y - y_mean

  if(!is.null(lambda)){
    # if lambda is given
    nlam <- length(lambda)
  }
  else{
    # max value of lambda for SQRT-lasso as the starting value
    lam_max <- max(abs(crossprod(x, y))) / sqrt(n * sum(y^2))
    lambda <- lam_max* exp(seq(0, log(flmin), length = nlam))
  }
  # since the est of the first lambda value does not make sense (all zero betas)
  # we append to the list of lambda an arbitrary value to start
  # and later we just delete that value
  lambda_c <- c(lambda[1] + 1e-4, lambda)
  # nlam_c is the nlam passed to C sub_routine
  nlam_c <- nlam + 1
  # to compute maxnz, we still use nlam instead of nlam + 1, since the first beta est is always all zeros
  # should be min(n, p), but this could potentially lead to
  # memory problem, so we make maxnz larger
  maxnz <- max(n, p) * nlam

  if (length(y) != n)
    stop("Dimensions of x and y do not match!")
  out = .C("olasso_c",
           # raw data
           as.double(x),
           # response
           as.double(y),
           # dimensions
           as.integer(n),
           as.integer(p),
           # current implementation does not estimate intercept
           # so do not need to standardize the columns
           # lambda path is always provided
           as.integer(nlam_c),
           # lambda sequence
           as.double(lambda_c),
           # beta in compressed column form
           be.i = as.integer(rep(0, maxnz)),
           be.p = as.integer(rep(0, nlam_c + 1)),
           be.x = as.double(rep(0, maxnz)),
           # number of non zeros
           as.integer(maxnz),
           # convergence threshold
           as.double(thresh))

  nz <- out$be.p[length(out$be.p)]
  # get sparse matrix representation of beta
  beta <- sparseMatrix(i=out$be.i[1:nz],
                       p=out$be.p, x=out$be.x[1:nz],
                       dims=c(p, nlam_c),
                       index1 = F)[, -1]
  # beta_est is the beta estimate in the original scale
  beta_est <- beta / col_sd

  if (intercept)
    a0 <- y_mean - as.numeric(crossprod(beta_est, col_mean))
  else
    a0 <- rep(0, p)

  fitted <- x %*% beta
  residual <- matrix(y, nrow = n, ncol = nlam) - fitted
  sse <- colSums(residual^2)

  # compute the degrees of freedom (nnz)
  df <- colSums(abs(beta) >= 1e-8)
  df[df >= n] <- n - 1

  sig_obj_path <- sqrt(sse / n + 2 * lambda * colSums(abs(beta))^2)
  sig_naive_path <- sqrt(sse / n)
  sig_df_path <- sqrt(sse / (n - df))

  out <- list(n = n, p = p, lambda = lambda,
              beta = beta_est,
              a0 = a0,
              sig_obj_path = sig_obj_path,
              sig_naive_path = sig_naive_path,
              sig_df_path = sig_df_path,
              type = "organic")
  class(out) <- "natural.path"
  return(out)
}

#' Solve organic lasso problem with a single value of lambda
#' The lambda values are for slow rates, which could give less satisfying results
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param thresh Threshold value for underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.
olasso_slow <- function(x, y, thresh = 1e-8){
  n <- nrow(x)
  p <- ncol(x)

  # lambda = (lam_1, lam_2)
  lambda <- getLam_slasso(n, p)

  if (lambda[1] < lambda[2]){
    lambda <- rev(lambda)
    ind_1 <- 2
    ind_2 <- 1
  }
  else {
    ind_1 <- 1
    ind_2 <- 2
  }

  # now lambda[ind_1] = lam_1, lambda[ind_2] = lam_2
  fit <- olasso_path(x = x, y = y,
                     lambda = lambda, thresh = thresh)
  return(list(n = n, p = p,
              lam_1 = lambda[ind_1],
              lam_2 = lambda[ind_2],
              a0_1 = fit$a0[ind_1],
              a0_2 = fit$a0[ind_2],
              beta_1 = fit$beta[, ind_1],
              beta_2 = fit$beta[, ind_2],
              sig_obj_1 = fit$sig_obj[ind_1],
              sig_obj_2 = fit$sig_obj[ind_2]))
}

#' Error standard deviation estimation using organic lasso
#'
#' Solve the organic lasso problem
#' \eqn{\tilde{\sigma}^2_{\lambda} = \min_{\beta} ||y - X \beta||_2^2 / n + 2 \lambda ||\beta||_1^2}
#' with two pre-specified values of tuning parameter:
#' \eqn{\lambda_1 = log p / n}, and \eqn{\lambda_2}, which is a Monte-Carlo estimate of \eqn{||X^T e||_\infty^2 / n^2}, where \eqn{e} is n-dimensional standard normal.
#'
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param intercept Indicator of whether intercept should be fitted. Default to be \code{TRUE}.
#' @param thresh Threshold value for underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.
#' @return A list object containing: \describe{
#' \item{\code{n} and \code{p}: }{The dimension of the problem.}
#' \item{\code{lam_1}, \code{lam_2}: }{\eqn{log(p) / n}, and an Monte-Carlo estimate of \eqn{||X^T e||_\infty^2 / n^2}, where \eqn{e} is n-dimensional standard normal.}
#' \item{\code{a0_1}, \code{a0_2}: }{Estimate of intercept, corresponding to \code{lam_1} and \code{lam_2}.}
#' \item{\code{beta_1}, \code{beta_2}: }{Organic lasso estimate of regression coefficients, corresponding to \code{lam_1} and \code{lam_2}.}
#' \item{\code{sig_obj_1}, \code{sig_obj_2}: }{Organic lasso estimate of the error standard deviation, corresponding to \code{lam_1} and \code{lam_2}.}}
#' @examples
#' set.seed(123)
#' sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
#' ol <- olasso(x = sim$x, y = sim$y[, 1])
#' @seealso \code{\link{olasso_path}, \link{olasso_cv}}
#' @export
olasso <- function(x, y, intercept = TRUE, thresh = 1e-8){
  n <- nrow(x)
  p <- ncol(x)

  # lambda = (lam_1, lam_2)
  lambda <- getLam_olasso(x)

  if (lambda[1] < lambda[2]){
    lambda <- rev(lambda)
    ind_1 <- 2
    ind_2 <- 1
  }
  else {
    ind_1 <- 1
    ind_2 <- 2
  }

  # now lambda[ind_1] = lam_1, lambda[ind_2] = lam_2
  fit <- olasso_path(x = x, y = y,
                     lambda = lambda, thresh = thresh,
                     intercept = intercept)
  return(list(n = n, p = p,
              lam_1 = lambda[ind_1],
              lam_2 = lambda[ind_2],
              a0_1 = fit$a0[ind_1],
              a0_2 = fit$a0[ind_2],
              beta_1 = fit$beta[, ind_1],
              beta_2 = fit$beta[, ind_2],
              sig_obj_1 = fit$sig_obj[ind_1],
              sig_obj_2 = fit$sig_obj[ind_2]))
}

#' Cross-validation for organic lasso
#'
#' Provide organic lasso estimate (of the error standard deviation) using cross-validation to select the tuning parameter value
#'
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{flmin}.
#' @param intercept Indicator of whether intercept should be fitted. Default to be \code{TRUE}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param flmin The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2}.
#' @param nfold Number of folds in cross-validation. Default value is 5. If each fold gets too view observation, a warning is thrown and the minimal \code{nfold = 3} is used.
#' @param foldid A vector of length \code{n} representing which fold each observation belongs to. Default to be \code{NULL}, and the program will generate its own randomly.
#' @param thresh Threshold value for underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.
#' @return A list object containing: \describe{
#' \item{\code{n} and \code{p}: }{The dimension of the problem.}
#' \item{\code{lambda}: }{The path of tuning parameter used.}
#' \item{\code{beta}: }{Estimate of the regression coefficients, in the original scale, corresponding to the tuning parameter selected by cross-validation.}
#' \item{\code{a0}: }{Estimate of intercept}
#' \item{\code{mat_mse}: }{The estimated prediction error on the test sets in cross-validation. A matrix of size \code{nlam} by \code{nfold}}
#' \item{\code{cvm}: }{The averaged estimated prediction error on the test sets over K folds.}
#' \item{\code{cvse}: }{The standard error of the estimated prediction error on the test sets over K folds.}
#' \item{\code{ibest}: }{The index in \code{lambda} that attains the minimal mean cross-validated error.}
#' \item{\code{foldid}: }{Fold assignment. A vector of length \code{n}.}
#' \item{\code{nfold}: }{The number of folds used in cross-validation.}
#' \item{\code{sig_obj}: }{Organic lasso estimate of the error standard deviation, selected by cross-validation.}
#' \item{\code{sig_obj_path}: }{Organic lasso estimates of the error standard deviation. A vector of length \code{nlam}.}
#' \item{\code{type}: }{whether the output is of a natural or an organic lasso.}}
# #' \item{\code{sig_naive}: }{Naive estimate of the error standard deviation based on the organic lasso regression, selected by cross-validation.}
# #' \item{\code{sig_naive_path}: }{Naive estimate of the error standard deviation based on the organic lasso regression. A vector of length \code{nlam}.}
# #' \item{\code{sig_df}: }{Degree-of-freedom adjusted estimate of standard deviation of the error, selected by cross-validation. See Reid, et, al (2016).}
# #' \item{\code{sig_df_path}: }{Degree-of-freedom adjusted estimate of standard deviation of the error. A vector of length \code{nlam}.}}
#' @examples
#' set.seed(123)
#' sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
#' ol_cv <- olasso_cv(x = sim$x, y = sim$y[, 1])
#' @seealso \code{\link{olasso_path}, \link{olasso}}
#' @export
olasso_cv <- function(x, y, lambda = NULL,
                      intercept = TRUE,
                      nlam = 100, flmin = 1e-2,
                      nfold = 5, foldid = NULL,
                      thresh = 1e-8){
  n <- nrow(x)
  p <- ncol(x)
  if(!is.null(lambda)){
    # if lambda is given
    nlam <- length(lambda)
  }
  else{
    # standardize the design matrix to have column means zero
    # and column sds 1
    stan <- standardize(x)
    xx <- stan$x

    if (intercept)
      yy <- y - mean(y)

    # use the max lambda for SQRT-lasso
    lam_max <- max(abs(crossprod(xx, yy))) / sqrt(n * sum(yy^2))
    lambda <- lam_max * exp(seq(0, log(flmin), length = nlam))
  }
  if (is.null(foldid)){
    # foldid is a vector of values between 1 and nfold
    # identifying what fold each observation is in.
    # If supplied, nfold can be missing.
    if (n / nfold < 10){
      warning("too few observations, use nfold = 3")
      nfold = 3
    }
    foldid <- sample(rep(seq(nfold), length = n))
    while(all.equal(sort(unique(foldid)), seq(nfold)) != TRUE){
      foldid <- sample(rep(seq(nfold), length = n))
    }
  }

  # use mse as the universal way of cross-validation criterion
  mat_mse <- matrix(NA, nrow = nlam, ncol = nfold)
  for (i in seq(nfold)){
    # train on all but i-th fold
    id_tr <- (foldid != i)
    id_te <- (foldid == i)
    # training/testing data set
    x_tr <- x[id_tr, ]
    x_te <- x[id_te, ]
    y_tr <- y[id_tr]
    y_te <- y[id_te]
    n_te <- sum(id_te)
    # get the fit using olasso on training data
    fit_tr <- olasso_path(x = x_tr, y = y_tr,
                  lambda = lambda, thresh = thresh)
    # and fit/obj on the test data
    mat_mse[, i] <- as.numeric(colMeans((matrix(y_te, nrow = n_te, ncol = nlam) - t(fit_tr$a0 + t(x_te %*% fit_tr$beta)))^2))
  }

  # extract information from CV
  # the mean cross-validated error, a vector of length nlam
  cvm <- rowMeans(mat_mse)
  # the index of best lambda
  ibest <- which.min(cvm)

  cvse <- apply(mat_mse, 1, stats::sd) / sqrt(ncol(mat_mse))

  # now fit the full data
  final <- olasso_path(x = x, y = y,
                  lambda = lambda, thresh = thresh,
                  intercept = intercept)

  out <- list(n = n, p = p, lambda = final$lambda,
              beta = final$beta[, ibest],
              a0 = final$a0[ibest],
              mat_mse = mat_mse,
              cvm = cvm,
              cvse = cvse,
              ibest = ibest,
              foldid = foldid,
              nfold = nfold,
              sig_obj = final$sig_obj[ibest],
              sig_obj_path = final$sig_obj,
              #sig_naive = final$sig_naive[ibest],
              #sig_naive_path = final$sig_naive,
              #sig_df = final$sig_df[ibest],
              #sig_df_path = final$sig_df
              type = "organic"
              )

  class(out) <- "natural.cv"
  return(out)
}
