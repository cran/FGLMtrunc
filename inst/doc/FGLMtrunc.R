## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
hook_output <- knitr::knit_hooks$get("output")
# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
        
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ---- eval=F------------------------------------------------------------------
#  install.packages("FGLMtrunc")

## -----------------------------------------------------------------------------
library(FGLMtrunc)

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
data(LinearExample)
Y_linear = LinearExample$Y
Xcurves_linear = LinearExample$X.curves
timeGrid = seq(0, 1, length.out = 101)
plot(timeGrid, LinearExample$beta.true, type = 'l', 
     main = 'True coefficient function', xlab = "t", ylab=expression(beta(t)))

## -----------------------------------------------------------------------------
fit = fglm_trunc(Y_linear, Xcurves_linear, nbasis = 50)

## ---- eval=F------------------------------------------------------------------
#  library(doMC)
#  registerDoMC(cores = 2)
#  fit = fglm_trunc(Y_linear, Xcurves_linear, nbasis = 50, parallel = TRUE)

## ---- eval=F------------------------------------------------------------------
#  k <- 50 - 3 - 1 #Numbers of knots = nbasis - degree - 1
#  knots_n <- seq(0, 1, length.out = k+2)[-c(1, k+2)] # Remove boundary knots
#  fit2 = fglm_trunc(Y_linear, Xcurves_linear, grid = timeGrid, knots = knots_n)

## -----------------------------------------------------------------------------
print(fit)

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
plot(fit)

## -----------------------------------------------------------------------------
predict(fit, newX.curves = Xcurves_linear[1:5,])

## ----out.lines = 12-----------------------------------------------------------
predict(fit, type = "coefficients")

## -----------------------------------------------------------------------------
data(LogisticExample)
Y_logistic = LogisticExample$Y
Xcurves_logistic = LogisticExample$X.curves

## -----------------------------------------------------------------------------
fit4 = fglm_trunc(Y_logistic, Xcurves_logistic, family="binomial", nbasis = 50)

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
print(fit4)
plot(fit4)

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
logistic_link_pred = predict(fit4, newX.curves = Xcurves_logistic, type="link")
plot(logistic_link_pred, ylab="log-odds")

## ---- fig.align='center', fig.height=4, fig.width=4---------------------------
logistic_response_pred = predict(fit4, newX.curves = Xcurves_logistic, type="response")
plot(logistic_response_pred, ylab="predicted probabilities")

## -----------------------------------------------------------------------------
scalar_coef <- c(1, -1, 0.5) # True coefficients for scalar predictors
set.seed(1234)
S <- cbind(matrix(rnorm(400), nrow=200), rbinom(200, 1, 0.5))  # Randomly generated observations for scalar predictors. Binary coded as 0 and 1.
colnames(S) <- c("s1", "s2", "s3")

## -----------------------------------------------------------------------------
Y_scalar <- Y_linear + (S %*% scalar_coef)

## -----------------------------------------------------------------------------
fit_scalar = fglm_trunc(X.curves=Xcurves_linear, Y=Y_scalar, S=S, nbasis = 50)
fit_scalar

## -----------------------------------------------------------------------------
predict(fit_scalar, newX.curves = Xcurves_linear[1:5,], newS=S[1:5,])

