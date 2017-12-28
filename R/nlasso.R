#' Fit a linear model with natural lasso
#'
#' Calculate a solution path of the natural lasso estimate (of error standard deviation) with a list of tuning parameter values. In particular, this function solves the lasso problems and returns the lasso objective function values as estimates of the error variance:
#' \eqn{\hat{\sigma}^2_{\lambda} = \min_{\beta} ||y - X \beta||_2^2 / n + 2 \lambda ||\beta||_1.}
#' The output also includes a path of naive estimates and a path of degree of freedom adjusted estimates of the error standard deviation.
#'
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{flmin}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param flmin The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2}.
#' @param thresh Threshold value for the underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.
#' @param intercept Indicator of whether intercept should be fitted. Dafault to be \code{TRUE}.
#' @return A list object containing: \describe{
#' \item{\code{n} and \code{p}: }{The dimension of the problem.}
#' \item{\code{lambda}: }{The path of tuning parameters used.}
#' \item{\code{beta}: }{Matrix of estimates of the regression coefficients, in the original scale. The matrix is of size \code{p} by \code{nlam}. The \code{j}-th column represents the estimate of coefficient corresponding to the \code{j}-th tuning parameter in \code{lambda}.}
#' \item{\code{a0}: }{Estimate of intercept. A vector of length \code{nlam}.}
#' \item{\code{sig_obj_path}: }{Natural lasso estimates of the error standard deviation. A vector of length \code{nlam}.}
#' \item{\code{sig_naive_path}: }{Naive estimates of the error standard deviation based on lasso regression, i.e., \eqn{||y - X \hat{\beta}||_2 / \sqrt n}. A vector of length \code{nlam}.}
#' \item{\code{sig_df_path}: }{Degree-of-freedom adjusted estimate of standard deviation of the error. A vector of length \code{nlam}. See Reid, et, al (2016).}
#' \item{\code{type}: }{whether the output is of a natural or an organic lasso.}}
#' @examples
#' set.seed(123)
#' sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
#' nl_path <- nlasso_path(x = sim$x, y = sim$y[, 1])
#' @import glmnet
#' @export
#' @seealso \code{\link{nlasso_cv}}
nlasso_path<- function(x, y, lambda = NULL,
                   nlam = 100, flmin = 1e-2,
                   thresh = 1e-8, intercept = TRUE){
  n <- nrow(x)
  p <- ncol(x)
  # we always standardize the design matrix X so that
  # ||X_j||_2^2 = n, where X_j is the j-th column of X
  stan <- standardize(x, center = intercept)
  x <- stan$x
  col_mean <- stan$col_mean
  col_sd <- stan$col_sd

  if (intercept)
    y <- y - mean(y)

  if (is.null(lambda)){
    lam_max <- max(abs(crossprod(x, y))) / n
    lambda <- lam_max * exp(seq(0, log(flmin), length = nlam))
  }
  else{
    nlam <- length(lambda)
  }

  # fit a lasso
  fit <- glmnet(x = x, y = y, lambda = lambda,
                intercept = FALSE, standardize = FALSE,
                thresh = thresh)

  colnames(fit$beta) <- NULL
  row.names(fit$beta) <- NULL
  # beta_est is the beta estimate in the original scale
  beta_est <- fit$beta / col_sd
  if (intercept)
    a0 <- mean(y) + as.numeric(crossprod(beta_est, col_mean))
  else
    a0 <- rep(0, p)
  # get a prediction
  fitted <- x %*% fit$beta

  residual <- matrix(y, nrow = n, ncol = nlam) - fitted
  sse <- as.numeric(colSums(residual^2))

  # note that it is possible that nnz returned by lasso
  # greater(due to numerical error) or equal to n
  # which makes the denominator less than or equal to 0
  df <- fit$df
  df[df >= n] <- n - 1

  sig_obj_path <- sqrt(sse / n + 2 * lambda * colSums(abs(fit$beta)))
  sig_naive_path <- sqrt(sse / n)
  sig_df_path <- sqrt(sse / (n - df))

  out <- list(n = n, p = p, lambda = lambda,
         beta = beta_est,
         a0 = a0,
         sig_obj_path = sig_obj_path,
         sig_naive_path = sig_naive_path,
         sig_df_path = sig_df_path,
         type = "natural")
  class(out) <- "natural.path"
  return(out)
}

#' Cross-validation for natural lasso
#'
#' Provide natural lasso estimate (of the error standard deviation) using cross-validation to select the tuning parameter value
#' The output also includes the cross-validation result of the naive estimate and the degree of freedom adjusted estimate of the error standard deviation.
#'
#' @param x An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.
#' @param y A response vector of size \code{n}.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{flmin}.
#' @param intercept Indicator of whether intercept should be fitted. Dafault to be \code{TRUE}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param flmin The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2}.
#' @param nfold Number of folds in cross-validaton. Default value is 5. If each fold gets too view observation, a warning is thrown and the minimal \code{nfold = 3} is used.
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
#' \item{\code{sig_obj}: }{Natural lasso estimate of standard deviation of the error, with the optimal tuning parameter selected by cross-validation.}
#' \item{\code{sig_obj_path}: }{Natural lasso estimates of standard deviation of the error. A vector of length \code{nlam}.}
#' \item{\code{sig_naive}: }{Naive estimates of the error standard deviation based on lasso regression, i.e., \eqn{||y - X \hat{\beta}||_2 / \sqrt n}, selected by cross-validation.}
#' \item{\code{sig_naive_path}: }{Naive estimate of standard deviation of the error based on lasso regression. A vector of length \code{nlam}.}
#' \item{\code{sig_df}: }{Degree-of-freedom adjusted estimate of standard deviation of the error, selected by cross-validation. See Reid, et, al (2016).}
#' \item{\code{sig_df_path}: }{Degree-of-freedom adjusted estimate of standard deviation of the error. A vector of length \code{nlam}.}
#' \item{\code{type}: }{whether the output is of a natural or an organic lasso.}}
#' @examples
#' set.seed(123)
#' sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
#' nl_cv <- nlasso_cv(x = sim$x, y = sim$y[, 1])
#' @seealso \code{\link{nlasso_path}}
#' @export
nlasso_cv <- function(x, y, lambda = NULL,
                      intercept = TRUE,
                      nlam = 100, flmin = 1e-2,
                      nfold = 5, foldid = NULL,
                      thresh = 1e-8){
  n <- nrow(x)
  p <- ncol(x)
  # first standardize the design matrix to have column means zero
  # and column sds 1
  stan <- standardize(x, center = intercept)
  x <- stan$x
  col_mean <- stan$col_mean
  col_sd <- stan$col_sd

  if (intercept)
    y <- y - mean(y)

  if(!is.null(lambda)){
    # if lambda is given
    nlam <- length(lambda)
  }
  else{
    # use the max lambda for lasso
    lam_max <- max(abs(crossprod(x, y))) / n
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

  # mse of lasso estimate of beta
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
    fit_tr <- nlasso_path(x = x_tr, y = y_tr,
                  lambda = lambda, thresh = thresh,
                  intercept = intercept)
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
  final <- nlasso_path(x = x, y = y,
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
              sig_obj = final$sig_obj_path[ibest],
              sig_obj_path = final$sig_obj_path,
              sig_naive = final$sig_naive_path[ibest],
              sig_naive_path = final$sig_naive_path,
              sig_df = final$sig_df_path[ibest],
              sig_df_path = final$sig_df_path,
              type = "natural")

  class(out) <- "natural.cv"
  return(out)
}
