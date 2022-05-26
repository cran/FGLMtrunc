#' @title Plot functional parameters \eqn{\beta} from a \code{FGLMtrunc} object
#'
#' @description Plot functional parameters \eqn{\beta} as a function of \eqn{t} for a fitted
#' \code{FGLMtrunc} object.
#'
#' @param x fitted \code{FGLMtrunc} object
#' @param include_smooth If TRUE, smoothing estimate without truncation of \eqn{\beta} is plotted.
#' @param \dots additional plot arguments
#' @return No return value.
#'
#' @method plot FGLMtrunc
#' @importFrom graphics legend lines title
#' @export
plot.FGLMtrunc <- function(x, include_smooth=TRUE, ...){
  plot(x$grid, x$beta.truncated, col="red", type="l",
       main=paste("Truncation point:", x$trunc.point), xlab="", ylab="", ...,
       ylim=c(min(c(x$beta.truncated, x$beta.0)), max(c(x$beta.truncated, x$beta.0)))
       )
  title(xlab=expression(t), ylab=expression(hat(beta)(t)), line=2)
  abline(h=0, col="black", lty=3)
  if (include_smooth) {
    lines(x$grid, x$beta.0, col="blue", lty=2)
    legend("bottomright", legend=c("Truncated", "Smoothing"),
           col=c("red", "blue"), lty=1:2, cex=0.8, title="method")

  }
  NULL
}
