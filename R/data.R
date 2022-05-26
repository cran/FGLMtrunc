#' @title Simulated data for functional linear regression.
#'
#' @description Randomly generated data with Gaussian responses for functional linear regression example follows Case I from Liu et. al. (2022).
#'
#' @name LinearExample
#' @docType data
#' @usage data(LinearExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'  \item{X.curves}{200 by 101 matrix of functional predictors.}
#'  \item{Y}{200 by 1 numeric vector of Gaussian responses.}
#'  \item{beta.true}{The true functional parameter \eqn{\beta}.}
#' }
#' @references Xi Liu, Afshin A. Divani, and Alexander Petersen. "Truncated estimation in functional generalized linear regression models" (2022). \emph{Computational Statistics & Data Analysis}.
NULL



#' Simulated data for functional logistic regression.
#'
#' Randomly generated data with Bernoulli responses for functional logistic regression example follows Case I from Liu et. al. (2022).
#'
#' @name LogisticExample
#' @docType data
#' @usage data(LogisticExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'  \item{X.curves}{200 by 101 matrix of functional predictors.}
#'  \item{Y}{200 by 1 numeric vector of Bernoulli responses.}
#'  \item{beta.true}{The true functional parameter \eqn{\beta}.}
#' }
#' @references Xi Liu, Afshin A. Divani, and Alexander Petersen. "Truncated estimation in functional generalized linear regression models" (2022). \emph{Computational Statistics & Data Analysis}.
NULL
