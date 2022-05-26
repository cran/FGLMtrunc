#' @title Fit a truncated Functional Generalized Linear Model
#'
#' @description Fit a truncated functional linear or logistic regression model using nested group lasso penalty.
#' The solution path is computed efficiently using active set algorithm with warm start. Optimal tuning parameters (\eqn{\lambda_s, \lambda_t})
#' are chosen by Bayesian information criterion (BIC).
#'
#' @details ## Details on spline estimator
#' For an order \code{q} B-splines (\code{q = degree + 1} since an intercept is used) with \code{k} internal knots 0 < \code{t_1} <...< \code{t_k} < T,
#' the number of B-spline basis equals \code{q + k}. Without truncation (\eqn{\lambda}_t=0), the function returns smoothing estimate that is
#' equivalent to the method of Cardot and Sarda (2005), and optimal smoothing parameter is chosen by Generalized Cross Validation (GCV).
#'
#' ## Details on \code{family}
#' The model can work with Gaussian or Bernoulli responses. If \code{family="gaussian"}, identity link is used. If \code{family="binomial"}, logit link is used.
#'
#' ## Details on scalar predictors
#' \code{FGLMtrunc} allows using scalar predictors together with functional predictors. If scalar predictors are used, their estimated coefficients
#' are included in \code{alpha} form fitted model.
#'
#' @param Y \code{n}-by-\code{1} vector of response.
#' Each row is an observed scalar response, which is continous for family="gaussian" and binary (i.e. 0 and 1) for family="binomal".
#' @param X.curves \code{n}-by-\code{p} matrix of functional predictors.
#' Each row is an observation vector at \code{p} finite points on \code{[0,T]} for some \code{T>0}.
#' @param S (optional) \code{n}-by-\code{s} matrix of scalar predictors. Binary variable should be coded as numeric rather than factor.
#' @param grid A sequence of \code{p} points at which \code{X} is recorded, including both boundaries \code{0} and \code{T}. If not
#' specified, an equally spaced sequence of length p between 0 and 1 will be used.
#' @param family Choice of exponential family for the model. The function then uses corresponding canonical link function to fit model.
#' @param degree Degree of the piecewise polynomial. Default 3 for cubic splines.
#' @param nbasis Number of B-spline basis.
#' If \code{knots} is unspecified, the function choose \code{nbasis - degree - 1} **internal** knots at suitable quantiles of \code{grid}.
#' If \code{knots} is specified, the value of \code{nbasis} will be **ignored**.
#' @param knots \code{k} **internal** breakpoints that define that spline.
#' @param nlambda.s (optional) Length of sequence of smoothing regularization parameters. Default 10.
#' @param lambda.s.seq (optional) Sequence of smoothing regularization parameters.
#' @param precision (optional) Error tolerance of the optimization. Default 1e-5.
#' @param parallel (optional) If TRUE, use parallel \code{foreach} to fit each value of \code{lambda.s.seq}. Must register parallel before hand, such as doMC or others.
#'
#'
#' @references Xi Liu, Afshin A. Divani, and Alexander Petersen. "Truncated estimation in functional generalized linear regression models" (2022). \emph{Computational Statistics & Data Analysis}.
#' @references Hervé Cardot and Pacal Sarda. "Estimation in generalized linear models for functional data via penalized likelihood" (2005). \emph{Journal of Multivariate Analysis}.
#'
#' @return A list with components:
#' \item{grid}{The \code{grid} sequence used.}
#' \item{knots}{The \code{knots} sequence used.}
#' \item{degree}{The degree of the piecewise polynomial used.}
#' \item{eta.0}{Estimate of B-spline coefficients \eqn{\eta} **without** truncation penalty.}
#' \item{beta.0}{Estimate of functional parameter \eqn{\beta} **without** truncation penalty.}
#' \item{eta.truncated}{Estimate of B-spline coefficients \eqn{\eta} **with** truncation penalty.}
#' \item{beta.truncated}{Estimate of functional parameter \eqn{\beta} **with** truncation penalty.}
#' \item{lambda.s0}{Optimal smoothing regularization parameter **without** truncation chosen by GCV.}
#' \item{lambda.s}{Optimal smoothing regularization parameter **with** truncation chosen by BIC.}
#' \item{lambda.t}{Optimal truncation regularization parameter chosen by BIC.}
#' \item{trunc.point}{Truncation point \eqn{\delta} where \eqn{\beta(t)} = 0 for \eqn{t \ge \delta}.}
#' \item{alpha}{Intercept (and coefficients of scalar predictors if used) of truncated model.}
#' \item{scalar.pred}{Logical variable indicating whether any scalar predictor was used.}
#' \item{call}{Function call of fitted model.}
#' \item{family}{Choice of exponential family used.}
#'
#' @seealso \link[splines2]{bSpline} from \link[splines2]{splines2} R package for usage of B-spline basis.
#'
#' @useDynLib FGLMtrunc
#' @import stats
#' @importFrom graphics abline
#' @import foreach
#' @importFrom splines2 bSpline
#'
#' @examples
#' 
#' # Gaussian response
#' data(LinearExample)
#' Y_linear = LinearExample$Y
#' Xcurves_linear = LinearExample$X.curves
#' fit1 = fglm_trunc(Y_linear, Xcurves_linear, nbasis = 20, nlambda.s = 1)
#' print(fit1)
#' plot(fit1)
#' 
#' @export
fglm_trunc <- function(
  Y,
  X.curves,
  S=NULL,
  grid=NULL,
  family = c("gaussian", "binomial"),
  degree = 3,
  nbasis=NULL,
  knots=NULL,
  nlambda.s = 10,
  lambda.s.seq=NULL,
  precision = 1e-5,
  parallel = FALSE){

  this.call=match.call()
  family = match.arg(family, family, several.ok = F)

  if (degree<2) {
    stop("Method requires degree to be at least 2.\n")
  }

  if (is.null(grid)) {
    #cat("grid is not specified! An equally spaced sequence of length p between 0 and 1 will be used.\n")
    grid = seq(0, 1, length.out = dim(X.curves)[2])
  }


  if (is.null(knots) & is.null(nbasis)) {
    stop("Either knots or nbasis needs to be specified.\n")
  } else if (is.null(knots) & !is.null(nbasis)) {
    if (nbasis <= degree) {
      stop("nbasis is too small.\n")
    }
  } else if (!is.null(knots)) {
    k = length(knots)
    nbasis = k + degree + 1
  }

  fit <- switch(family,
                "gaussian" = truncatedLinearModelCpp(Y, X.curves, S, grid, degree, nbasis, knots, lambda.s.seq, nlambda.s, precision, parallel),
                "binomial" = truncatedLogisticModelCpp(Y, X.curves, S, grid, degree, nbasis, knots, lambda.s.seq, nlambda.s, precision, parallel)
                )
  if (!is.null(S)){
    fit$scalar.pred <- TRUE
  } else {
    fit$scalar.pred <- FALSE
  }
  fit$family = family
  fit$call = this.call
  class(fit) = "FGLMtrunc"

  return(fit)
}
