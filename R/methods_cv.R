#' plot a natural.cv object
#'
#' This function is adapted from the \pkg{ggb} R package.
#'
#' @param x an object of class \code{natural.cv}, as returned by
#'        \code{\link{nlasso_cv}} and \code{\link{olasso_cv}}
#' @param ... additional argument(not used here, only for S3 generic/method consistency)
#' @method plot natural.cv
#' @export
plot.natural.cv <- function(x, ...){
  par(mar = c(5, 5, 5, 1))
  plot(x$lambda, x$cvm,
       main = paste("Cross-validation plot of ", x$type, " lasso", sep = ""),
       xlab = "tuning parameter",
       ylab = "Cross-validation error",
       type="o",
       ylim = c(min(x$cvm - x$cvse), max(x$cvm + x$cvse)),
       pch = 20)

  # 1-se rule
  lines(x$lambda, x$cvm + x$cvse)
  lines(x$lambda, x$cvm - x$cvse)
  abline(v = x$lambda[x$ibest], lty = 2)
}
