#' @useDynLib FGLMtrunc
truncatedLogisticModelCpp <- function(Y,
                                      X.curves,
                                      S,
                                      grid,
                                      degree,
                                      nbasis,
                                      knots,
                                      lambda_s_seq,
                                      nlambda_s,
                                      precision,
                                      parallel) {
  n = length(Y)
  if (!is.null(S)){
    p.scalar = ncol(S)
  } else {
    p.scalar = 0
  }
  p.funcs = ncol(X.curves)

  # step 1 construct bspline basis function, convert curves to vector and get matrics, M and N's #
  resList = curves2scalarsVecLogistic(X.curves, S, grid, nbasis, knots, degree)
  scalar_mat = resList$scalar_mat # n x (1 + q+k)
  M_aug = resList$M_aug
  xi_mat = resList$xi_mat
  bspline.basis = resList$bspline.basis

  # step 2 estimate weight vector #
  # we denote the estimation and tuning parameter with "0" to indicate they are in the stage of estimating weight vector

  ## step 2.1 GCV to choose lambda_s0
  lambda_s0 = gcvSmoothLogistic(Y, scalar_mat, xi_mat, M_aug, grid, nbasis, figure=F)

  ## step 2.2 estimate \eta with lambda_s0

  warmstart0 = glmnet(scalar_mat[,-1], Y, family = "binomial", alpha = 0, lambda = lambda_s0)$beta
  warmstart0 = c( 0, as.matrix(warmstart0))
  beta0 = logisticSmPenalty(Y, scalar_mat, M_aug, lambda_s = lambda_s0, warmstart = warmstart0, precision = 1e-8)$beta_vec
  eta0 = beta0[(length(beta0) - nbasis + 1):length(beta0), ]

  ## step 2.3 compute weight vector
  weight = 1/sapply(1:(nbasis - degree), function(i) norm(as.matrix(eta0[i:nbasis]), type = "f")^2)

  # step 3 determine optimal lambda_s and lambda_t #
  if (is.null(lambda_s_seq)) {
    lambda_s_seq = seq(lambda_s0/n/20, lambda_s0/n/2, length.out = nlambda_s)
  } else {
    nlambda_s = length(lambda_s_seq)
  }

  esti_list = as.list(seq(nlambda_s))
  if (parallel) {
    esti_list = foreach(i = seq(nlambda_s), .packages = c("FGLMtrunc")) %dopar%
      logisticSolutionPath(Y, p.funcs, scalar_mat, xi_mat, nbasis, M_aug, lambda_s_seq[i], weight, degree, p.scalar, precision)
  } else {
    for (i in seq(nlambda_s)) {
      esti_list[[i]] = logisticSolutionPath(Y, p.funcs, scalar_mat, xi_mat, nbasis, M_aug, lambda_s_seq[i], weight, degree, p.scalar, precision)
    }
  }
  
  lambda_t_sse_vec = sapply(esti_list, get, x="bic")
  s_index = which.min(lambda_t_sse_vec)
  lambda_s_opti = lambda_s_seq[s_index]
  lambda_t_opti = esti_list[[s_index]]$lambda_t

  # step 4 estimate beta using lambda_s_opti and lambda_t_opti #
  beta_vec = esti_list[[s_index]]$beta_path_mat
  eta_truncated = beta_vec[(length(beta0) - nbasis + 1):length(beta0)]
  alpha = beta_vec[-((length(beta0) - nbasis + 1):length(beta0))]
  beta.truncated = c(bspline.basis %*% eta_truncated)

  if (!is.null(S)){
    names(alpha) = c("Intercept", colnames(data.frame(S)))
  }
  # return a list #
  res_list = list(grid = grid,
                  knots=resList$knots,
                  degree=degree,
                  bspline.basis = bspline.basis,

                  eta.0 = eta0,
                  beta.0 = c(bspline.basis %*% eta0),

                  eta.truncated = eta_truncated,
                  beta.truncated = beta.truncated,

                  lambda.s0 = lambda_s0,
                  lambda.s = lambda_s_opti,
                  lambda.t = lambda_t_opti,

                  trunc.point = grid[which(beta.truncated==0)[1]],
                  alpha=alpha
  )
  return(res_list)
}
