#' @title Make predictions from \code{FGLMtrunc} fitted model
#'
#' @description This function returns truncated estimate of linear predictors, fitted values, and functional parameter \eqn{\beta}
#' for a fitted \code{FGLMtrunc} object.
#'
#'
#' @param object fitted \code{FGLMtrunc} object
#' @param newX.curves Matrix of new values for functional predictors \code{X.curves}.
#' @param newS Matrix of new values for scalar predictors \code{S}.
#' @param type Type of prediction. For logistic regression (\code{family = "binomial"}), \code{type="link"} gives the linear
#' predictors, which is log-odds, and \code{type="response"} gives the predicted probabilities.
#' For linear regression (\code{family = "gaussian"}), both \code{type="link"} and  \code{type="response"} give fitted values.
#' For both linear regression and logistic regression, \code{type="coefficients"} gives truncated estimate of functional parameter \eqn{\beta}.
#' @param \dots additional predict arguments (Not applicable for FGLMtrunc)
#' @return Predictions depends on chosen \code{type}.
#' @seealso \link[stats]{predict.glm}.
#' @method predict FGLMtrunc
#' @export
predict.FGLMtrunc <- function (object, newX.curves, newS=NULL, type = c("link", "response", "coefficients"), ...) {
  type = match.arg(type, type, several.ok = F)

  if(missing(object)) {
    stop("Fitted truncated FGLM model is missing!\n")
  }

  if (type != "coefficients") {
    if(missing(newX.curves)) {
      stop("newX.curves is missing!\n")
    }

    if (object$scalar.pred & is.null(newS)){
      stop("Model is fitted with scalar predictors by newS is missing!\n")
    }

    if (object$family == "gaussian") {
      nbasis = length(object$knots) + object$degree + 1
      res = curves2scalarsVecLinear(X.curves=newX.curves, S=newS, grid=object$grid, nbasis=nbasis,knots=object$knots, degree=object$degree)
      out = res$scalar_mat %*% matrix(c(object$alpha, object$eta.truncated), ncol=1)

    } else if (object$family == "binomial") {
      nbasis = length(object$knots) + object$degree + 1
      res = curves2scalarsVecLogistic(X.curves=newX.curves, S=newS, grid=object$grid, nbasis=nbasis, knots=object$knots, degree=object$degree)
      theta = res$scalar_mat %*% matrix(c(object$alpha, object$eta.truncated), ncol=1)
      if (type == "link") {
        out = theta
      } else if (type == "response") {
        out = 1 / (1 + exp(-theta))
      }
    }
  } else {
    out = object$beta.truncated
  }
  return(out)
}

