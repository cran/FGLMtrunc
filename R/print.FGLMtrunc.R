#' @title Print a \code{FGLMtrunc} object
#'
#' @description Print a summary of truncation point of the fitted \code{FGLMtrunc} model.
#' @details
#' Truncation point estimate of \eqn{\delta} is printed.
#' @param x fitted \code{FGLMtrunc} object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @return The fitted object is silently return.
#' @method print FGLMtrunc
#' @export
print.FGLMtrunc <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call, width.cutoff = 100), "\n\n")
  if (x$scalar.pred) {
    printCoefmat(data.frame(t(x$alpha) , row.names = ""))
  }
  cat("\nOptimal truncation point:", x$trunc.point, "\n")
  invisible(x)
}
