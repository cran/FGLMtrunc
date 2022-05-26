#' @useDynLib FGLMtrunc
gcvSmoothlinear <- function(Y, scalar_mat, xi_mat, M_aug, denseGrid, nbasis, figure = F) {

  n = length(Y)
  trace_ratio = sum(diag(t(xi_mat) %*% xi_mat))/sum(diag(M_aug))
  gcv_error_vec = NULL
  lambda_s_seq = seq(trace_ratio/10, trace_ratio*100, length.out = 1000)

  for (lambda in lambda_s_seq) {
    gcv_error = linearSmPenalty(Y, scalar_mat = scalar_mat, M_aug = M_aug, lambda_s = lambda)$gcv
    gcv_error_vec = c(gcv_error_vec, gcv_error)
  }
  min_index = which.min(gcv_error_vec)

  if (figure) {
    plot(lambda_s_seq, gcv_error_vec, type = "l", main = "spline FGLM -- BIC(lambda_s)")
    abline(v = lambda_s_seq[min_index], col = 2)
  }
  return(lambda_s_seq[min_index])
}




#' @useDynLib FGLMtrunc
#' @importFrom glmnet glmnet
gcvSmoothLogistic <- function(Y, scalar_mat, xi_mat, M_aug, denseGrid, nbasis, figure = F) {

  n = length(Y)
  trace_ratio = sum(diag(t(xi_mat) %*% xi_mat))/sum(diag(M_aug))
  lambda_s_seq = seq(trace_ratio/10, trace_ratio*100, length.out = 1000)
  warmstart0 = glmnet(scalar_mat[,-1], Y, family = "binomial", alpha = 0, lambda = lambda_s_seq[1])$beta
  warmstart0 = c( 0, as.matrix(warmstart0))

  res = logisticSmPenalty(Y, scalar_mat, M_aug, lambda_s = lambda_s_seq[1], warmstart = warmstart0, precision = 1e-6)
  gcv_error_vec = res$bic

  for (lambda in lambda_s_seq[-1]) {
    gcv_error = logisticSmPenalty(Y, scalar_mat = scalar_mat, M_aug = M_aug, lambda_s = lambda, warmstart = res$beta)$bic
    gcv_error_vec = c(gcv_error_vec, gcv_error)
  }
  min_index = which.min(gcv_error_vec)

  if (figure) {
    plot(lambda_s_seq, gcv_error_vec, type = "l", main = "spline FGLM -- GCV(lambda_s)")
    abline(v = lambda_s_seq[min_index], col = 2)
  }

  return(lambda_s_seq[min_index])

}

