% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/olasso_api.R
\name{olasso}
\alias{olasso}
\title{Error standard deviation estimation using organic lasso}
\usage{
olasso(x, y, intercept = TRUE, thresh = 1e-08)
}
\arguments{
\item{x}{An \code{n} by \code{p} design matrix. Each row is an observation of \code{p} features.}

\item{y}{A response vector of size \code{n}.}

\item{intercept}{Indicator of whether intercept should be fitted. Default to be \code{TRUE}.}

\item{thresh}{Threshold value for underlying optimization algorithm to claim convergence. Default to be \code{1e-8}.}
}
\value{
A list object containing: \describe{
\item{\code{n} and \code{p}: }{The dimension of the problem.}
\item{\code{lam_1}, \code{lam_2}: }{\eqn{log(p) / n}, and an Monte-Carlo estimate of \eqn{||X^T e||_\infty^2 / n^2}, where \eqn{e} is n-dimensional standard normal.}
\item{\code{a0_1}, \code{a0_2}: }{Estimate of intercept, corresponding to \code{lam_1} and \code{lam_2}.}
\item{\code{beta_1}, \code{beta_2}: }{Organic lasso estimate of regression coefficients, corresponding to \code{lam_1} and \code{lam_2}.}
\item{\code{sig_obj_1}, \code{sig_obj_2}: }{Organic lasso estimate of the error standard deviation, corresponding to \code{lam_1} and \code{lam_2}.}}
}
\description{
Solve the organic lasso problem
\eqn{\tilde{\sigma}^2_{\lambda} = \min_{\beta} ||y - X \beta||_2^2 / n + 2 \lambda ||\beta||_1^2}
with two pre-specified values of tuning parameter:
\eqn{\lambda_1 = log p / n}, and \eqn{\lambda_2}, which is a Monte-Carlo estimate of \eqn{||X^T e||_\infty^2 / n^2}, where \eqn{e} is n-dimensional standard normal.
}
\examples{
set.seed(123)
sim <- make_sparse_model(n = 50, p = 200, alpha = 0.6, rho = 0.6, snr = 2, nsim = 1)
ol <- olasso(x = sim$x, y = sim$y[, 1])
}
\seealso{
\code{\link{olasso_path}, \link{olasso_cv}}
}
