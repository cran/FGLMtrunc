#include <iostream>
#include <RcppArmadillo.h>
#include "math.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List linearSmPenalty (const arma::vec& Y,
                            const arma::mat& scalar_mat,
                            const arma::mat& M_aug,
                            double lambda_s) {

  double n = Y.n_rows;
  arma::mat invTemp = arma::inv(scalar_mat.t() * scalar_mat + (n * lambda_s) * M_aug) * scalar_mat.t();
  arma::vec beta_vec = invTemp * Y;
  double df_lambda = arma::trace(scalar_mat * invTemp);

  arma::mat Y_fit = scalar_mat * beta_vec;
  double rss = arma::sum(arma::pow(Y-Y_fit.col(0), 2));
  double gcv = n*log(rss/n) + df_lambda*log(n);

  return Rcpp::List::create(Rcpp::Named("beta_vec") = beta_vec,
                            Rcpp::Named("gcv") = gcv);
}

// accelerated gradient descent combined with smooth dual problem
// the input betaA serves as a warm start
Rcpp::List lipschitzL_C_linear (double n,
                                const arma::mat& X_active,
                                const arma::mat& M_mat_active,
                                double mu,
                                int q_A,
                                double lambda_s,
                                double lambda_t,
                                const arma::vec& weight) {
  arma::mat C(q_A, q_A, arma::fill::zeros);
  C.diag().fill(weight[0]);
  for (int i=1; i<q_A; i++) {
    arma::mat Ci(q_A-i, q_A-i, arma::fill::zeros);
    Ci.diag().fill(weight[i]);
    arma::mat J(q_A-i, q_A, arma::fill::zeros);
    J.tail_cols(q_A-i) = Ci;
    C = join_cols(C, J);
  }
  C = lambda_t * C;
  arma::vec s = svd(C) ;
  double Lmu = pow(s[0], 2)/mu;

  arma::mat temp_mat = (X_active.t() * X_active/n) + (lambda_s*M_mat_active);
  arma::vec eigval = arma::eig_sym(temp_mat);
  double Lc = arma::as_scalar(eigval.tail(1)) + Lmu;
  return Rcpp::List::create(Rcpp::Named("Lc") = Lc,
                            Rcpp::Named("C") = C);
}

arma::vec alpha_star_linear (const arma::vec& etaA,
                             double mu,
                             const arma::vec& weight,
                             int q_A,
                             double lambda_t) {
  double group_norm_A;
  arma::vec alpha_star, temp;
  for (int i = 0; i<q_A; i++) {
    group_norm_A = arma::norm(etaA.tail(q_A-i), "fro");
    if (group_norm_A >  mu/lambda_t/arma::as_scalar(weight[i])) {
      temp = etaA.tail(q_A-i)/group_norm_A;
    } else {
      temp = lambda_t*arma::as_scalar(weight[i])*etaA.tail(q_A-i)/mu;
    }
    alpha_star = arma::join_cols(alpha_star, temp);
  }

  return(alpha_star);
}

void acceleratedGD_linear(const arma::vec& Y,
                          int n,
                          const arma::mat& X_active,
                          const arma::mat& M_mat_active,
                          int d_active,
                          int q_A,
                          int p_scalar,
                          double lambda_s,
                          double lambda_t,
                          const arma::vec& weight,
                          arma::vec& betaA,
                          int maxIter = 50000,
                          double convergePre = 1e-5,
                          double deltaPre = 5e-2) {

  int p = betaA.n_elem;
  int iter = 0;
  double mu = deltaPre/q_A; // when mu is smaller, the smooth begins later, smooth fun is the same with original fun
  arma::vec w_aux = betaA;  //w_0

  Rcpp::List L_C_res = lipschitzL_C_linear(n, X_active, M_mat_active, mu, q_A, lambda_s, lambda_t, weight);
  double Lc = L_C_res[0];
  arma::mat C = L_C_res[1];
  double stepsize = 1./Lc;

  arma::vec alpha_star = alpha_star_linear(w_aux.tail(p-p_scalar-1), mu, weight, q_A, lambda_t);
  arma::vec zero_V(p_scalar+1, arma::fill::zeros);
  arma::vec grad = -X_active.t() * (Y - X_active * betaA) / n + lambda_s * (M_mat_active * betaA) + arma::join_cols(zero_V, C.t() * alpha_star);
  arma::vec beta_old;
  double theta;
  while ((arma::norm(grad, "fro") > convergePre) && (iter < maxIter)) {
    beta_old = betaA; 
    betaA = w_aux - stepsize*grad;

    alpha_star = alpha_star_linear(w_aux.tail(p-p_scalar-1), mu, weight, q_A, lambda_t);
    theta = iter/(iter+3);
    w_aux = betaA + theta*(betaA - beta_old);

    grad = - X_active.t() * (Y - X_active * w_aux) / n + lambda_s * (M_mat_active * w_aux) + arma::join_cols(zero_V, C.t() * alpha_star);
    iter += 1;
  }
  int iter_2 = 1;
  while ((arma::norm(grad, "fro") > convergePre) && (iter_2 < maxIter)) {
    grad = - (X_active.t()) * (Y - X_active * betaA) / n + lambda_s * (M_mat_active * betaA) + arma::join_cols(zero_V, C.t() * alpha_star);
    betaA -= stepsize*grad;
    iter_2 += 1;
  }
  //return (beta_new);
}

arma::vec penalty_d2(const arma::vec& beta_A,
                         int nbasis,
                         int q_A,
                         const arma::vec& weight,
                         int degree) {

  arma::vec beta_rm_scalar = beta_A.tail(q_A);
  arma::vec group_norm_A(q_A);
  arma::vec dPe2_vec_A(q_A);
  for (int i = 0; i<q_A; i++) {
    group_norm_A[i] = arma::norm(beta_rm_scalar.subvec(i,q_A-1), "fro");
  }
  for (int i = 0; i<q_A; i++) {
    dPe2_vec_A[i] = pow(arma::sum(beta_rm_scalar[i] /  group_norm_A.head(i+1) % weight.head(i+1)), 2);
  }

  arma::uvec idx = arma::conv_to<arma::uvec>::from(arma::linspace(q_A, nbasis-degree-1, abs(nbasis-degree-1-q_A)+1));
  double weight_total_sum = arma::sum(weight(idx));
  arma::vec temp(degree);
  temp.fill(weight_total_sum);
  arma::vec dPe2_vec_notA = arma::pow(arma::join_cols(arma::cumsum(weight(idx)),
                                                      temp), 2);
  arma::vec res = arma::cumsum(arma::join_cols(dPe2_vec_A,
                                               dPe2_vec_notA));
  return(res);
}

// a piece of linear solution path -- fixed active set
// [[Rcpp::export]]
Rcpp::List linearpiecePathCpp (const arma::vec& Y,
                               int n,
                               const arma::mat& scalar_mat,
                               const arma::mat& M_aug,
                               const arma::vec& warmStart,
                               int nbasis,
                               const arma::vec& weight,
                               double lambda_s,
                               double lambdaStart,
                               int d_active,
                               int p_scalar,
                               int degree,
                               double precision = 1e-3) {
  int q_A = d_active - p_scalar - 1;
  int iter = 1;
  double delta_lambda = -lambdaStart/10.;
  arma::mat X_active = scalar_mat.head_cols( d_active );
  arma::mat M_mat_active = M_aug( arma::span(0, d_active-1), arma::span(0, d_active-1) );

  double lambda_new = lambdaStart + delta_lambda;
  arma::vec beta_A = warmStart.head(d_active);
  acceleratedGD_linear(Y, n, X_active, M_mat_active, d_active, q_A, p_scalar, lambda_s, lambda_new , weight, beta_A , 100000);
  arma::mat beta_A_mat(beta_A.n_rows, 1);
  beta_A_mat.col(0) = beta_A;

  arma::vec dl_part2 = lambda_s * M_aug * arma::join_cols(beta_A, arma::vec(scalar_mat.n_cols - d_active, arma::fill::zeros));
  arma::vec dl_part1 = scalar_mat.t() * (Y - X_active * beta_A) / n; 
  int p = dl_part1.n_elem;
  arma::vec diff_v = arma::pow(dl_part1 - dl_part2, 2);

  arma::rowvec dL_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(arma::cumsum( diff_v.tail(p-p_scalar-1) )) );
  arma::mat dL_group_mat(1, dL_group.n_elem);
  dL_group_mat.row(0) = dL_group;

  arma::rowvec dPe_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(penalty_d2(beta_A, nbasis, q_A, weight, degree) * (lambda_new * lambda_new)) );
  arma::mat dPe_group_mat(1, dPe_group.n_elem);
  dPe_group_mat.row(0) =  dPe_group;

  Rcpp::NumericVector lambda_seq = {lambda_new}; // start from lambdaStart + delta_lambda

  while (arma::all(dL_group.subvec(d_active-1, nbasis-1) < dPe_group.subvec(d_active-1, nbasis-1)) && (lambda_new >= precision/n))  { //# all block dL < dPe

    while ((lambda_new + delta_lambda) < ( precision/n/2)) {
      delta_lambda = delta_lambda/5.;
    }
    lambda_new = lambda_new + delta_lambda;
    acceleratedGD_linear(Y, n, X_active, M_mat_active, d_active, q_A, p_scalar, lambda_s, lambda_new , weight, beta_A, 100000);

    beta_A_mat = arma::join_rows(beta_A_mat, beta_A);

    dl_part2 = lambda_s * M_aug * arma::join_cols(beta_A, arma::vec(scalar_mat.n_cols - d_active, arma::fill::zeros));
    dl_part1 = scalar_mat.t() * (Y - X_active * beta_A) / n;
    diff_v = arma::pow(dl_part1 - dl_part2, 2);
    dL_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(arma::cumsum( diff_v.tail(p-p_scalar-1) )) );
    dL_group_mat = join_cols(dL_group_mat, dL_group);

    dPe_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(penalty_d2(beta_A, nbasis, q_A, weight, degree) * (lambda_new * lambda_new)) );
    dPe_group_mat = join_cols(dPe_group_mat, dPe_group);
    lambda_seq.push_back(lambda_new);
    iter = iter + 1;
  }

    return Rcpp::List::create(Rcpp::Named("dL_dPe_ratio") = (dL_group_mat/dPe_group_mat).t(),
                              Rcpp::Named("betaPath") = arma::join_cols(beta_A_mat, arma::mat(scalar_mat.n_cols-d_active, lambda_seq.length(), arma::fill::zeros)),
                              Rcpp::Named("lambda_seq") = lambda_seq);
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_df_linear (const Rcpp::NumericVector& d_active_seq,
                                       const arma::mat& scalar_mat,
                                       const arma::mat& M_aug,
                                       double lambda_s,
                                       int n) {
  Rcpp::NumericVector df_v;
  int l = d_active_seq.length();
  arma::mat U_a, M_a;
  double tr;
  int d_active;
  for (int i=0; i<l; i++){
    d_active = d_active_seq[i];
    U_a = scalar_mat.head_cols(d_active);
    M_a = M_aug.submat( 0, 0, d_active-1, d_active-1 );
    tr = arma::trace(U_a * arma::inv(U_a.t() * U_a + n * lambda_s * M_a) * U_a.t());
    df_v.push_back(tr);
  }
  return df_v;
}



// [[Rcpp::export]]
Rcpp::List logisticSmPenalty(const arma::vec& Y,
                             const arma::mat& scalar_mat,
                             const arma::mat& M_aug,
                             double lambda_s,
                             const arma::vec& warmstart,
                             double precision = 1e-6) {
  int n = Y.n_rows;
  arma::vec beta = warmstart;
  arma::vec theta_vec = scalar_mat * beta;
  arma::vec b_prime = arma::exp(theta_vec) / (1. + arma::exp(theta_vec));
  arma::vec grad = - scalar_mat.t() * (Y - b_prime) + lambda_s * (M_aug * beta);

  int iterCount = 1;
  arma::vec b_2prime;
  arma::mat b_2prime_diagmat;
  arma::mat hessianInv;


  // newton raphson
  while((arma::norm(grad, "fro") >  precision) && (iterCount < 1e3)) {
    b_2prime = b_prime % (1. - b_prime);
    b_2prime_diagmat = scalar_mat.t() * arma::diagmat(b_2prime) * scalar_mat;
    hessianInv =  arma::inv(b_2prime_diagmat + lambda_s* M_aug);
    beta -= hessianInv * grad;
    theta_vec = scalar_mat * beta;
    b_prime = arma::exp(theta_vec) / (1. + arma::exp(theta_vec));
    grad = - scalar_mat.t() * (Y - b_prime) + lambda_s * (M_aug * beta);
    iterCount = iterCount + 1;
  }

  //compute the IC value
  //compute dimension of trace(H(lambda_s))
  double dim_beta = arma::trace(hessianInv * b_2prime_diagmat);

  //compute deviance D = -2 * loglikelihood
  arma::vec b_theta = arma::log(1. + arma::exp(theta_vec));
  double deviance_logistic = -2 * arma::sum( (Y%theta_vec) - b_theta);

  //IC = deviance + log(n) * dim_beta => smaller is better
  double bic = deviance_logistic + log(n)*dim_beta;

  return Rcpp::List::create(Rcpp::Named("beta_vec") = beta,
                            Rcpp::Named("bic") = bic,
                            Rcpp::Named("dim_beta") = dim_beta,
                            Rcpp::Named("deviance_logistic") = deviance_logistic);

}



void acceleratedGD_logistic (const arma::vec& Y,
                             int n,
                             const arma::mat& X_active,
                             const arma::mat& M_mat_active,
                             int d_active,
                             int q_A,
                             int p_scalar,
                             double lambda_s,
                             double lambda_t,
                             const arma::vec& weight,
                             arma::vec& betaA,
                             int maxIter = 50000,
                             double convergePre = 1e-5,
                             double deltaPre = 5e-3) {

  int p = betaA.n_elem;
  int iter = 0;
  double mu = deltaPre/q_A; //when mu is smaller, the smooth begins later, smooth fun is the same with original fun
  arma::vec w_aux = betaA;  //w_0

  Rcpp::List L_C_res = lipschitzL_C_linear(n, X_active, M_mat_active, mu, q_A, lambda_s, lambda_t, weight);
  double Lc = L_C_res[0];
  arma::mat C = L_C_res[1];
  double stepsize = 1./Lc;

  arma::vec alpha_star = alpha_star_linear(w_aux.tail(p-p_scalar-1), mu, weight, q_A, lambda_t);
  arma::vec theta_vec = X_active * betaA;
  arma::vec b_prime = arma::exp(theta_vec)/(1 + arma::exp(theta_vec));
  arma::vec zero_V(p_scalar+1, arma::fill::zeros);
  arma::vec grad = - X_active.t() * (Y - b_prime) /n + lambda_s*(M_mat_active * betaA) + arma::join_cols(zero_V, C.t() * alpha_star);
  arma::vec beta_old;
  double theta;
  while ((arma::norm(grad, "fro") > convergePre) && (iter <  maxIter)) {

    beta_old = betaA;
    betaA = w_aux - stepsize*grad;

    alpha_star = alpha_star_linear(w_aux.tail(p-p_scalar-1), mu, weight, q_A, lambda_t);
    theta = iter/(iter+3);
    w_aux = betaA + theta*(betaA - beta_old);

    theta_vec = X_active * w_aux;
    b_prime = arma::exp(theta_vec) / (1. + arma::exp(theta_vec));
    grad = - X_active.t() * (Y - b_prime)/n + lambda_s * (M_mat_active * w_aux) + arma::join_cols(zero_V, C.t() * alpha_star);
    iter = iter + 1;
  }

  int iter_2 = 1;

  while ((arma::norm(grad, "fro") > convergePre) && (iter_2 < 3*maxIter)) {
    alpha_star = alpha_star_linear(betaA.tail(p-p_scalar-1), mu, weight, q_A, lambda_t);
    theta_vec = X_active * betaA;
    b_prime = arma::exp(theta_vec) / (1. + arma::exp(theta_vec));
    grad = - X_active.t() * (Y - b_prime) /n + lambda_s * (M_mat_active * betaA) + arma::join_cols(zero_V, C.t() * alpha_star);
    betaA = betaA - stepsize*grad;
    iter_2 = iter_2 + 1;
  }
}

// [[Rcpp::export]]
Rcpp::List logisticpiecePathCpp (const arma::vec& Y,
                                 int n,
                                 const arma::mat& scalar_mat,
                                 const arma::mat& M_aug,
                                 const arma::vec& warmStart,
                                 int nbasis,
                                 const arma::vec& weight,
                                 double lambda_s,
                                 double lambdaStart,
                                 int d_active,
                                 int p_scalar,
                                 int degree,
                                 double precision = 1e-3) {
  int q_A = d_active - p_scalar - 1;
  int iter = 1;
  double delta_lambda = -lambdaStart/10.;
  arma::mat X_active = scalar_mat.head_cols(d_active);
  arma::mat M_mat_active = M_aug( arma::span(0, d_active-1), arma::span(0, d_active-1) );

  double lambda_new = lambdaStart + delta_lambda;
  arma::vec beta_A = warmStart.head(d_active);
  acceleratedGD_logistic(Y, n, X_active, M_mat_active, d_active, q_A, p_scalar, lambda_s, lambda_new , weight, beta_A, 50000);
  arma::mat beta_A_mat(beta_A.n_rows, 1);
  beta_A_mat.col(0) = beta_A;

  arma::vec dl_part2 = lambda_s * (M_aug * arma::join_cols(beta_A, arma::vec(scalar_mat.n_cols - d_active, arma::fill::zeros))); 
  arma::vec theta_vec = X_active * beta_A;
  arma::vec b_prime = arma::exp(theta_vec)/(1. + arma::exp(theta_vec));

  arma::vec dl_part1 = scalar_mat.t() * (Y - b_prime) /n; 
  int p = dl_part1.n_elem;
  arma::vec diff_v = arma::pow(dl_part1 - dl_part2, 2);
  arma::rowvec dL_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(arma::cumsum( diff_v.tail(p-p_scalar-1) )) );
  arma::mat dL_group_mat(1, dL_group.n_elem);
  dL_group_mat.row(0) = dL_group;

  arma::rowvec dPe_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(penalty_d2(beta_A, nbasis, q_A, weight, degree) * (lambda_new * lambda_new)) );
  arma::mat dPe_group_mat(1, dPe_group.n_elem);
  dPe_group_mat.row(0) =  dPe_group;

  Rcpp::NumericVector lambda_seq = {lambda_new}; // start from lambdaStart + delta_lambda

  while (arma::all(dL_group.subvec(d_active-1, nbasis-1) < dPe_group.subvec(d_active-1, nbasis-1)) && (lambda_new >= precision/n)) {

    while ((lambda_new + delta_lambda) < (precision/n/2)) {
      delta_lambda = delta_lambda/5.;
    }
    lambda_new = lambda_new + delta_lambda;
    acceleratedGD_logistic(Y, n, X_active, M_mat_active, d_active, q_A, p_scalar, lambda_s, lambda_new, weight, beta_A, 50000);
    beta_A_mat = arma::join_rows(beta_A_mat, beta_A);

    dl_part2 = lambda_s * (M_aug * arma::join_cols(beta_A, arma::vec(scalar_mat.n_cols - d_active, arma::fill::zeros)));
    theta_vec = X_active * beta_A;
    b_prime = arma::exp(theta_vec)/(1. + arma::exp(theta_vec));
    dl_part1 = scalar_mat.t() * (Y - b_prime) /n;
    diff_v = arma::pow(dl_part1 - dl_part2, 2);
    dL_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(arma::cumsum( diff_v.tail(p-p_scalar-1) )) );
    dPe_group = arma::conv_to<arma::rowvec>::from( arma::sqrt(penalty_d2(beta_A, nbasis, q_A, weight, degree) * (lambda_new * lambda_new)) );

    dL_group_mat = join_cols(dL_group_mat, dL_group);
    dPe_group_mat = join_cols(dPe_group_mat, dPe_group);
    lambda_seq.push_back(lambda_new);
    iter = iter + 1;
  }

  return Rcpp::List::create(Rcpp::Named("dL_dPe_ratio") = (dL_group_mat/dPe_group_mat).t(),
                            Rcpp::Named("betaPath") = arma::join_cols(beta_A_mat, arma::mat(scalar_mat.n_cols-d_active, lambda_seq.length(), arma::fill::zeros)),
                            Rcpp::Named("lambda_seq") = lambda_seq);
}

// [[Rcpp::export]]
double compute_dim_beta_logistic (const arma::mat& U_aDU_a,
                                  const arma::mat& hessian){
  double dim_beta = arma::trace( arma::inv(hessian) * U_aDU_a );
  return dim_beta;
}





