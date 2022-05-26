#' @useDynLib FGLMtrunc
#' @importFrom glmnet glmnet
logisticSolutionPath <- function(
  Y,
  p,
  scalar_mat,
  xi_mat,
  nbasis,
  M_aug,
  lambda_s,
  weight,
  degree,
  p.scalar,
  precision = 1e-4) {

  n = length(Y)

  logisticfit0 = stats::glm(Y ~ 1, family = "binomial")$coef
  beta_A = logisticfit0
  beta_path_mat = matrix(c(beta_A, rep(0, nbasis+p.scalar)) , ncol = 1)

  resi = Y - mean(Y)
  del_vec2 = (t(xi_mat) %*% (resi)/n)^2 # the loss function / n
  total_weight = sum(weight)
  pel_vec = (c(cumsum(weight), rep(total_weight, degree)))^2
  lambda_t0_seq = sqrt(cumsum(del_vec2)/cumsum(pel_vec))

  d_active_seq = 1

  lambda_t_seq =  max(lambda_t0_seq)
  d_active = 1 + which.max(lambda_t0_seq) +p.scalar
  d_active_seq = c(d_active_seq, d_active)

  change_index = 1

  while(lambda_t_seq[length(lambda_t_seq)] * n >= precision & d_active <= nbasis + p.scalar - (degree + 1) + 1) {
    lambda.start = lambda_t_seq[length(lambda_t_seq)]
    d_active = d_active_seq[length(d_active_seq)]
    warm.start = beta_path_mat[, ncol(beta_path_mat)]

    pieceRes <- logisticpiecePathCpp(Y, n, scalar_mat, M_aug, warm.start, nbasis, weight, lambda_s, lambda.start, d_active, p.scalar, degree, precision)

    lambda_t_seq = c(lambda_t_seq, pieceRes$lambda_seq)
    beta_path_mat = cbind(beta_path_mat, pieceRes$betaPath)
    change_index = c(change_index, length(lambda_t_seq))
    d_active = max(which(pieceRes$dL_dPe_ratio[, length(pieceRes$lambda_seq)] > 1)) + 1
    if (d_active > nbasis + p.scalar - (degree+1) + 1 |lambda_t_seq[length(lambda_t_seq)]  < precision/n) {
      d_active = nbasis + 1
    } else {
      d_active_seq = c(d_active_seq, d_active)
    }

  }
  beta0_warmstart = glmnet(scalar_mat[,-1], Y, family = "binomial", alpha = 0, lambda = lambda_s*n)$beta
  beta0_warmstart = c( 0, as.matrix(beta0_warmstart))
  beta0 = logisticSmPenalty(Y=Y, scalar_mat=scalar_mat, M_aug=M_aug, lambda_s = lambda_s*n, warmstart=beta0_warmstart)$beta_vec
  lambda_t_seq = c(lambda_t_seq, 0)
  beta_path_mat = cbind(beta_path_mat, beta0)
  change_index = c(change_index, ncol(beta_path_mat))
  d_active_seq = c(d_active_seq, d_active)

  scalar_idx = 1:(1+p.scalar)
  bic = sapply(1:length(d_active_seq), function(j) {
    d_active = d_active_seq[j]
    lambda_t = lambda_t_seq[change_index[j]]

    if (d_active>1 & d_active<p) {
      # approximate truncation penality
      eta = beta_path_mat[, change_index[j]][-scalar_idx]
      index_nonzero = d_active - p.scalar - 1
      norm_eta = sapply(1:index_nonzero, function(i) norm(as.matrix(eta[i:index_nonzero]), type = "f"))
      diag_ele = (lambda_t * sum(norm_eta*weight[1:index_nonzero]))/sum(eta^2)
      diag_vec = c(rep(0, length(scalar_idx)), rep(diag_ele, index_nonzero))
      M_t = diag(diag_vec)
    } else if (d_active == 1) {
      M_t = as.matrix(0)
    } else {
      eta = beta_path_mat[, change_index[j]][-scalar_idx]
      index_nonzero = d_active - p.scalar - 1
      norm_eta = sapply(1:(nbasis - degree), function(i) norm(as.matrix(eta[i:nbasis]), type = "f"))
      diag_ele = (lambda_t * sum(norm_eta*weight))/sum(eta^2)
      diag_vec = c(rep(0, length(scalar_idx)), rep(diag_ele, index_nonzero))
      M_t = diag(diag_vec)
    }

    U_a = scalar_mat[, 1:d_active]
    M_a = M_aug[1:d_active, 1:d_active]
    theta_vec = scalar_mat %*% beta_path_mat[, change_index[j]]
    
    b_prime = exp(theta_vec)/(1 + exp(theta_vec))
    b_2prime = b_prime*(1-b_prime)
    U_aDU_a = t(U_a) %*% diag(as.vector(b_2prime)) %*% U_a
    hessian = (U_aDU_a) + n* lambda_s *M_a + n* M_t
    dim_beta = compute_dim_beta_logistic(U_aDU_a, hessian)
    b_theta = log(1 + exp(theta_vec))
    deviance_logistic = - 2* sum(Y*theta_vec - b_theta)
    bic = deviance_logistic + log(n)*dim_beta

    return(bic)
  })

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
