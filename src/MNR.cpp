#include "genet_types.h"

//' @noRd
// [[Rcpp::export]]
SpMat compute_weighted_interactions(MSpMat H, VectorXi cards, bool dc, Rcpp::NumericVector target){
  SpMat interactions = H;
  for (int k=0; k<interactions.outerSize(); ++k){
    for (SpMat::InnerIterator it(interactions,k); it; ++it)
    {
      if (it.row()==it.col()) {
        if (dc) it.valueRef() = (cards(it.row()) - 1) * it.value();
        else it.valueRef() = (cards(it.row()) - 1) * target(it.row()) * it.value();
      } else {
        if (dc) it.valueRef() = cards(it.row()) * it.value();
        else it.valueRef() = cards(it.row()) * target(it.row()) * it.value();
      }
    }
  }
  return interactions;
}

// [[Rcpp::export]]
VectorXd initialize_root(MSpMat mapped_interactions, VectorXd target) {
  SpMat interactions = mapped_interactions;
  VectorXd a = interactions.diagonal(); // a as in quadratic formula
  VectorXd b = VectorXd::Zero(target.size()); // b as in quadratic formula
  VectorXd upper = interactions.diagonal().cwiseInverse();
  upper = target.cwiseProduct(upper);
  upper = upper.cwiseSqrt();
  // Rcpp::Rcout << "upper bounds:\n" << upper << "\n";
  for (int k=0; k<interactions.outerSize(); ++k){
    for (SpMat::InnerIterator it(interactions,k); it; ++it)
    {
      if (it.row() != it.col()) {
        b(it.col()) += upper(it.row()) * it.value();
      }
    }
  }
  VectorXd lower = b.cwiseAbs2() + 4 * a.cwiseProduct(target);
  lower = - b + lower.cwiseSqrt();
  lower = lower.cwiseProduct(0.5 * a.cwiseInverse());
  // Rcpp::Rcout << "lower bounds:\n" << lower << "\n";
  VectorXd root(target.size());
  for (int l = 0; l < target.size(); l++) {
    root(l) = Rcpp::as<double>(Rcpp::runif(1, lower(l), upper(l)));
  }
  return root;
}

//' @noRd
// [[Rcpp::export]]
void update_J_and_f(SpMat &J, VectorXd &f, MSpMat interactions, VectorXd root, VectorXd target) {
  J = interactions;
  f = VectorXd::Zero(target.size());
  VectorXd J_diag = VectorXd::Zero(root.size());
  for (int k=0; k<J.outerSize(); ++k){
    for (SpMat::InnerIterator it(J,k); it; ++it)
    {
      f(it.col()) += root(it.row()) * it.value();
      if (it.row() == it.col()) {
        it.valueRef() = 2 * root(it.row()) * it.value();
      } else {
        J_diag(it.col()) += root(it.row()) * it.value(); // external contributions to diagonal are saved in J_diag to avoid using expensive .coeffRef()
        it.valueRef() = root(it.col()) * it.value();
      }
    }
  }
  f = f.cwiseProduct(root) - target;
  J.diagonal() += J_diag; // finalize diagonal with contributions saved in J_diag
  J = J.transpose(); // align contributions
}

//' @noRd
// [[Rcpp::export]]
void check_convergence(bool &not_converged, VectorXd diff_root, VectorXd f, double eps) {
  double error = diff_root.norm() + f.norm();
  if (error < eps) not_converged = false;
}

//' @noRd
// [[Rcpp::export]]
VectorXd MNR(MSpMat interactions, VectorXd root_ini, VectorXd target, double eps = 10^(-5), bool verbose = false, int max_ite = 100) {
  VectorXd f, root = root_ini, root_new;
  SpMat J;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  bool not_converged = true;
  int iteration = 0;
  update_J_and_f(J, f, interactions, root, target);
  solver.analyzePattern(J);
  while (not_converged) {
    iteration += 1;
    if (verbose) Rcpp::Rcout << "Iteration "<< iteration << "\n";
    // if (verbose) Rcpp::Rcout << ".. J:\n" << J;
    solver.factorize(J);
    // solver.compute(J);
    f = J * root - f; // constant term
    // if (verbose) Rcpp::Rcout << ".. constant term:\n"<< f << "\n";
    root_new = solver.solve(f);
    // if (verbose) Rcpp::Rcout << ".. current root:\n"<< root_new << "\n";
    update_J_and_f(J, f, interactions, root_new, target);
    if (verbose) Rcpp::Rcout << ".. error of current root is "<< f.norm() << "\n";
    check_convergence(not_converged, root - root_new, f, eps);
    root = root_new;
    if (iteration == max_ite) not_converged = false;
    Rcpp::checkUserInterrupt();
  }
  if (verbose) {
    if (not_converged) Rcpp::Rcout << "MNR did not converge.\n";
    else Rcpp::Rcout << "MNR converged.\n";
  }
  return root;
}
