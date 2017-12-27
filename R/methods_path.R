#' print a natural.path object
#'
#' This function is adapted from the \pkg{ggb} R package.
#'
#' @param x an object of class \code{natural.path}, as returned by
#'        \code{\link{nlasso_path}} and \code{\link{olasso_path}}
#' @param ... additional argument(not used here, only for S3 generic/method consistency)
#'
#' @method print natural.path
#' @export
print.natural.path <- function(x, ...){
  cat(sprintf("%s lasso path with %s lambda values.", x$type,
              length(x$lambda)), fill = TRUE)
  tab <- data.frame(lambda = x$lambda,
                    sig_objective = x$sig_obj_path,
                    sig_naive = x$sig_naive_path,
                    sig_df = x$sig_df_path)
  print(tab, row.names = FALSE)
}

#' plot a natural.path object
#'
#' This function is adapted from the \pkg{ggb} R package.
#'
#' @param x an object of class \code{natural.path}, as returned by
#'        \code{\link{nlasso_path}} and \code{\link{olasso_path}}
#' @param ... additional argument(not used here, only for S3 generic/method consistency)
#' @method plot natural.path
#' @export
plot.natural.path <- function(x, ...){
  graphics::par(mar = c(5, 5, 5, 1))
  yrange <- range(x$sig_obj_path, x$sig_naive_path, x$sig_df_path)
  graphics::plot(x$lambda, x$sig_naive_path, type = "l", col = "black",
       lwd = 2,
       main = paste("Path plot of ", x$type, " lasso", sep = ""),
       ylim = yrange,
       ylab = "estimated error s.d.", xlab = "lambda")
  graphics::lines(x$lambda, x$sig_df_path, col = "red", lwd = 2)
  graphics::lines(x$lambda, x$sig_obj_path, col = "blue", lwd = 2)
  graphics::legend("bottomright",
                   legend = c("naive", "df", "natural"),
                   lwd = c(2, 2, 2),
                   col = c("black", "red", "blue"))
}
