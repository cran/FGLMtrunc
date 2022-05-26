#' @useDynLib FGLMtrunc
linearSolutionPath <- function(
  Y,
  scalar_mat,
  xi_mat,
  nbasis,
  M_aug,
  lambda_s,
  weight,
  degree,
  p.scalar,
  precision = 1e-4) {

  lmfit0 = stats::lm(Y ~ 1)
  resi = lmfit0$residual
  beta_A = lmfit0$coefficients
  beta_path_mat = matrix(c(beta_A, rep(0, nbasis+p.scalar)) , ncol = 1)
  n = length(Y)
  d_active_seq = 1

  del_vec2 = (t(xi_mat) %*% (resi)/n)^2
  total_weight = sum(weight)
  pel_vec = (c(cumsum(weight), rep(total_weight, degree)))^2
  lambda_t0_seq = sqrt(cumsum(del_vec2)/cumsum(pel_vec))

  lambda_t_seq = max(lambda_t0_seq)
  d_active = 1 + which.max(lambda_t0_seq) +p.scalar


  d_active_seq = c(d_active_seq, d_active)
  change_index = 1

  while(lambda_t_seq[length(lambda_t_seq)] * n >= precision & d_active <= nbasis + p.scalar - (degree + 1) + 1) {

    lambdaStart = lambda_t_seq[length(lambda_t_seq)]
    d_active = d_active_seq[length(d_active_seq)]
    warmStart = beta_path_mat[, ncol(beta_path_mat)]

    pieceRes <- linearpiecePathCpp(Y, n, scalar_mat, M_aug, warmStart, nbasis, weight, lambda_s, lambdaStart, d_active, p.scalar, degree, precision)

    lambda_t_seq = c(lambda_t_seq, pieceRes$lambda_seq)
    beta_path_mat = cbind(beta_path_mat, pieceRes$betaPath)
    change_index = c(change_index, length(lambda_t_seq))
    d_active = max(which(pieceRes$dL_dPe_ratio[, length(pieceRes$lambda_seq)] > 1)) + 1
    if (d_active > nbasis + p.scalar - (degree + 1) + 1 | lambda_t_seq[length(lambda_t_seq)]  < precision/n) {
      d_active = nbasis + 1
    } else {
      d_active_seq = c(d_active_seq, d_active)
    }
  }

  beta0 = linearSmPenalty(Y, scalar_mat = scalar_mat, M_aug = M_aug, lambda_s = lambda_s)$beta_vec
  lambda_t_seq = c(lambda_t_seq, 0)
  beta_path_mat = cbind(beta_path_mat, beta0)
  change_index = c(change_index, ncol(beta_path_mat))
  d_active_seq = c(d_active_seq, d_active)

  Y_pre = scalar_mat %*% beta_path_mat[, change_index]
  Y_true = matrix(rep(Y, dim(Y_pre)[2]), ncol = dim(Y_pre)[2])
  rss = colSums((Y_true - Y_pre)^2)

  df_lambda_t = compute_df_linear(d_active_seq,
                                  scalar_mat,
                                  M_aug,
                                  lambda_s,
                                  n)
                           
  bic = n*log(rss/n) + df_lambda_t*log(n)
  bic.min.idx = which.min(bic)
  bic.min = bic[bic.min.idx]
  d_active.min = d_active_seq[bic.min.idx]
  change_index.min = change_index[bic.min.idx]

  res_2 <- list(bic=bic.min,
                lambda_t = lambda_t_seq[change_index.min],
                d_active=d_active.min,
                beta_path_mat= beta_path_mat[,change_index.min])
  return(res_2)
}
