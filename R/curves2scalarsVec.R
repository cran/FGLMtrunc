#' @useDynLib FGLMtrunc
curves2scalarsVecLinear <- function(X.curves, S, grid, nbasis, knots, degree) {
  diff_t_vec = diff(grid)
  eval_grid = grid[-1] - diff_t_vec/2

  # NO. of basis = order + No. interior knots
  if (is.null(knots)) {
    k <- nbasis - degree - 1 # k = df - order
    knots <- seq(min(grid), max(grid), length.out = k+2)[-c(1, k + 2)] # the k interior knots
  }

  bspline.basis <- bSpline(grid, knots = knots, degree = degree, intercept = T)

  bspline.eval = predict(bspline.basis, newx = eval_grid)
  X.curves_evl = (X.curves[, 1:length(diff_t_vec)] + X.curves[, 2:length(grid)])/2
  xi_mat = X.curves_evl %*% diag(diff_t_vec) %*% bspline.eval

  deri_2 = predict(deriv(bspline.basis, derivs = 2), newx = eval_grid)
  M_mat = t(deri_2) %*% diag(diff_t_vec) %*%  deri_2

  scalarVec = cbind(1, S, xi_mat)

  dim_total = ncol(scalarVec)
  M_aug = diag(rep(0, dim_total))
  M_aug[(dim_total - nbasis + 1 ):dim_total, (dim_total - nbasis +1 ):dim_total] = M_mat

  trace_ratio = sum(diag(t(xi_mat) %*% xi_mat))/sum(diag(M_aug))
  M_aug = M_aug*trace_ratio

  list(M_aug = M_aug,
       scalar_mat = scalarVec,
       xi_mat = xi_mat,
       bspline.basis=bspline.basis,
       knots=knots)
}








#' @useDynLib FGLMtrunc
curves2scalarsVecLogistic <- function(X.curves, S, grid, nbasis, knots, degree) {
  diff_t_vec = diff(grid)
  eval_grid = grid[-1] - diff_t_vec/2

  # NO. of basis = order + No. interior knots
  if (is.null(knots)) {
    k <- nbasis - degree - 1 # k = df - order
    knots <- seq(min(grid), max(grid), length.out = k+2)[-c(1, k + 2)] # the k interior knots
  }

  bspline.basis <- bSpline(grid, knots = knots, degree = degree, intercept = T)

  bspline.eval = predict(bspline.basis, newx = eval_grid)
  X.curves_evl = (X.curves[, 1:length(diff_t_vec)] + X.curves[, 2:length(grid)])/2
  xi_mat = X.curves_evl %*% diag(diff_t_vec) %*% bspline.eval

  deri_2 = predict(deriv(bspline.basis, derivs = 2), newx = eval_grid)
  M_mat = t(deri_2) %*% diag(diff_t_vec) %*%  deri_2


  scalarVec = cbind(1, S, xi_mat)

  dim_total = ncol(scalarVec)
  M_aug = diag(rep(0, dim_total))
  M_aug[(dim_total - nbasis + 1 ):dim_total, (dim_total - nbasis +1 ):dim_total] = M_mat
  M_aug = M_aug/M_aug[dim_total, dim_total]


  list(M_aug = M_aug,
       scalar_mat = scalarVec,
       xi_mat = xi_mat,
       bspline.basis=bspline.basis,
       knots=knots)
}
